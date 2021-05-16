# standard libraries
from datetime import datetime
from itertools import combinations
import math
import os
import re
# 3rd party libraries
from Bio import SeqIO
import dendropy as ddp
import numpy as np
import pandas as pd
# my modules
from info_parsing_dicts import file_names_style
from poissonsampling import poisson_mode_and_alpha

"""
Module to make probability-model predictions and validate predictions against 
gold-standard, BEAST-estimated time-scaled phylogenies.

Takes in:
* BEAST trees file
* fasta file(s) containing tree sequences with identifiers containing 
  sampleid_yearfraction

1. make combinations of pairs of samples
2. create index where each level is a sample identifier in a pair of samples
3. add column with difference in time, Dt, between the index sample pair
4. calculate distance between the pair's samples
5. add putative times to common ancestor, Dt_pCA, given in argument
6. calculate cumulative evolutionary time , Dt_CE, 2 * Dt_pCA - Dt
7. calculate upper level of expected substitutions for the range of Dt_CEs
8. make model predictions for each region and pair
9. for each tree, determine whether the BEAST-estimated Most Recent Common 
   Ancestor is earlier or later than the putative CA
10. classify model predictions for each tree into true/false negative/positive
   (TN/TP/FN/FP)
11. append results to a csv file on disk (path created from input trees file) 
"""

# Rates csv columns
# -----------------
RATE = 'Rate'
RATES_REGION = 'Region'
SITES = 'Sites'
RATES_GENOTYPE = 'Genotype'

# Regex file info groups
# ----------------------
REGEX_GENOTYPE = 'Genotype'
REGEX_DATASET = 'Dataset'
REGEX_REGION = 'Region'
REGEX_SETTINGS = 'Settings'

# Data frame headings
# -------------------
# Dataset info columns
GENOTYPE = 'Genotype'
DATASET = 'Dataset'
REGION = 'Region'
TREE_NUMBER = 'Tree'
# sample ID columns
SAMPLE1 = 'Sample1'
SAMPLE2 = 'Sample2'
# sampling time columns
SAMPLE1_T = 't_s1 (weeks)'
SAMPLE2_T = 't_s2 (weeks)'
# time between samples (weeks) column
PAIR_DELTA = 'Dt (weeks)'
# time between most recent sample and Bayesian-estimated MRCA column
T_TO_EST_MRCA = 'Dt_bMRCA (weeks)'
# distance between samples column
DISTANCE = 'Distance'
# Poisson prediction columns
TIME = 'Time'
MODE = 'Mode'
LOWER = 'Lower'
UPPER = 'Upper'
# Modelling variables columns
T_TO_PUTATIVE_CA = 'Dt_pCA (weeks)'
PRED_EVOL_T = 'Dt_CE (weeks)'
# Result classification columns
ESTIMATED = 'Tree'
PREDICTED = 'Model'
RESULT = 'Result'
# Result categories
# -----------------
RELATED = True  # H0 True
UNRELATED = False  # H0 False
TRUE_POS = 'TP'
TRUE_NEG = 'TN'
FALSE_POS = 'FP'
FALSE_NEG = 'FN'
CORRECT = 'correct'
INCORRECT = 'incorrect'


# permissive matching, so long as there is at least one nucleotide in common
# between two bases, they are considered to match
labels = ['A', 'T', 'G', 'C', 'R', 'Y', 'W', 'M', 'S', 'K', 'B', 'V', 'H', 'D',
          'N', '-']

distances = np.array(
    # A  T  G  C  R  Y  W  M  S  K  B  V  H  D  N  -
    [[0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1],  # A
     [1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1],  # T
     [1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1],  # G
     [1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1],  # C
     [0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # R
     [1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # Y
     [0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],  # W
     [0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],  # M
     [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # S
     [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],  # K
     [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # B
     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # V
     [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # H
     [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # D
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # N
     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]   # -
     ],
    dtype='int8')
permissive_matches = pd.DataFrame(distances, index=labels, columns=labels)


def num_diff_seq_pair(seq1, seq2, distance_matrix=None):
    """
    Takes in two aligned sequences (with '-') and returns absolute number of
    differences between them. Counts gaps and difference in length between
    sequences as differences. If a partial_matches_allowed matrix is given, the
    distance between the pair of sequences will be calculated using it.

    Parameters
    ----------
    seq1, seq2 : str
        Two nucleotide sequences to compare.
    distance_matrix : None or pd.DataFrame
        If None, Hamming distance is calculated, otherwise distance is based on
        the matrix given. Data frame should be square, with characters to
        compare as indices and columns.

    Returns
    -------
    int
        Either a Hamming distance or the distance based on the matrix provided
        with difference in sequence length being counted as distance.
    """
    if distance_matrix is not None:
        distance = sum(
            distance_matrix.loc[a, b]
            for a, b in zip(seq1.upper(), seq2.upper()))
    else:
        distance = sum(
            1 for a, b in zip(seq1.upper(), seq2.upper())
            if a != b)
    return distance + abs(len(seq1) - len(seq2))


def pairs_index_from_iter(samples_iter):
    """Produces a pandas index object from an iterator of sample IDs."""
    return pd.Index(combinations(samples_iter, 2))


def pairs_estimated_times_to_mrca(tree, pairs_index=None):
    """
    Given a Bayesian-inferred phylogenetic tree returns a data frame containing
    the time elapsed between the most recent common ancestor estimated for each
    sample pair and the most recent sample in the pair.

    Parameters
    ----------
    tree : ddp.datamodel.treemodel.Tree
        A dendropy tree object.
    pairs_index : pd.Index or None, default None
        Pandas 2-level index of sample pairs (str, str).
        If None, an index is obtained from the leaf node labels in the tree.

    Returns
    -------
    pd.DataFrame
        * 2-level Index with each level being a sample (str) in a sample pair.
        * T_TO_EST_MRCA: weeks between the most recent sample and the Bayesian-
          estimated MRCA for each pair of samples (float).
    """
    # produce tree distance matrix
    distance_matrix = tree.node_distance_matrix()

    pairs_index = (pairs_index_from_iter(tree.taxon_namespace.labels())
                   if pairs_index is None
                   else pairs_index)
    pairs_data = pd.DataFrame(index=pairs_index)

    ids_nodes = {}

    def get_tree_time_to_mrca(row_index):
        """
        Returns estimated time to MRCA for the sample pair in the row index.
        """
        # get the tree nodes for the samples
        nodes_of_interest = []
        for sample in row_index:
            if sample in ids_nodes:
                sample_node = ids_nodes[sample]
            else:
                sample_node = tree.find_node_with_taxon_label(sample)
                ids_nodes[sample] = sample_node
            nodes_of_interest.append(sample_node)

        # get time of the MRCA as estimated by BEAST
        # * find MRCA node
        s1_s2_mrca = distance_matrix.mrca(*nodes_of_interest)
        nodes_of_interest.append(s1_s2_mrca)
        # * calculate the distance of each node from the root (cannot use
        #   distance_from_tip as this calculates time between node and most
        #   recent child of the node (not most recent tip in the tree);
        #   tree distances in years, so converted to weeks
        node_root_dist = {node: node.distance_from_root() * 52
                          for node in nodes_of_interest}
        most_recent_sample_node = max(node_root_dist,
                                      key=lambda node: node_root_dist[node])
        # * number of weeks between the most recent sample and the
        #   Bayesian-estimated time to the pair's MRCA
        delta_t_b_mrca = (node_root_dist[most_recent_sample_node] -
                          node_root_dist[s1_s2_mrca])
        return delta_t_b_mrca

    pairs_data[T_TO_EST_MRCA] = pairs_index.map(get_tree_time_to_mrca)
    return pairs_data


def year_fraction_from_id(identifier):
    """
    Takes sample record id like 'H123456789_2017.1489' and returns year
    fraction (float following '_').
    """
    return float(identifier.split('_')[-1])


def dates_from_fas_identifiers(identifiers):
    """
    Returns a dict of identifiers: dates from an iterable of identifiers with
    date.

    Parameters
    ----------
    identifiers : iter
        Iterable of identifiers like 'H123456789_2017.1489'.

    Returns
    -------
    dict
        identifier (str): date in weeks (float).
    """
    return {identifier: year_fraction_from_id(identifier) * 52  # wk/yr
            for identifier in identifiers}


def calculate_pairs_delta(pairs_index, date_info):
    """
    Given a pandas Index object and a dictionary containing date information
    for each index sample, returns a series with the time difference between
    each pair's samples.

    Parameters
    ----------
    pairs_index : pd.Index
        2-level pandas Index object with each level being a sample identifier
        in a sample pair (str).
    date_info : dict
        Dictionary mapping each sample in the index to year + wk / 52
        (str:float).

    Returns
    -------
    pd.DataFrame
        * 2-level pandas Index object with each level being a sample (str) in a
          sample pair.
        * PAIR_DELTA matching each index pair (float).
    """
    deltas_data = pairs_index.map(
        lambda ind: abs(date_info[ind[0]] - date_info[ind[1]]))
    return pd.DataFrame(data={PAIR_DELTA: deltas_data}, index=pairs_index)


def sample_pairs_distances(seq_records, pairs_data=None):
    """
    Returns a data frame with the pairwise distances for the records in the
    fasta file.

    Parameters
    ----------
    seq_records : dict
        Dictionary mapping a sequence record ID to a SeqIO record object.
    pairs_data : pd.DataFrame or None, default None
        Pandas data frame with sample pairs as index. If None, an Index is
        obtained from the records in the fasta file.
        * 2-level index: pairs' sample ids (str)

    Returns
    -------
    pd.DataFrame
        * 2-level index: pairs' sample ids (str)
        * ...any existing columns in pairs_data...
        * DISTANCE: sample pairs' distance (int)
    """
    pairs_data = (
        pd.DataFrame(index=pairs_index_from_iter(seq_records.keys()))
        if pairs_data is None
        else pairs_data)

    pairs_data[DISTANCE] = pairs_data.index.map(
        lambda ind: num_diff_seq_pair(
            *(seq_records[sample].seq for sample in ind),
            distance_matrix=permissive_matches))

    return pairs_data


def expected_substitutions(rate, sites, alpha, time_range):
    """
    Given a substitution rate, the number of sites, the alpha probability
    interval and the time range to calculate expected values for, returns
    expected, lower and upper values for the expected substitutions for a
    Poisson distribution with lambda=rate*time.

    Parameters
    ----------
    rate : float
        Substitution rate in substitutions/(site*year).
    sites : int
        Number of sites in the genomic region to which the rate applies.
    alpha : float
        Probability interval - fraction of the Poisson distribution to fall
        within interval (0-1). e.g., 0.95: 95% of the Poisson distribution at a
        given time will fall within the lower-upper range of number of
        substitutions.
    time_range : tuple
        Time range for which expected values are calculated (ints or floats),
        (min_time, max_time)

    Returns
    -------
    expected_data : pd.DataFrame
        Pandas data frame containing expected substitution range for each time
        in the interval given.
        TIME: time for which mode, lower and upper are calculated, in weeks
            (int)
        MODE: typical number of substitutions expected at each time (int)
        LOWER: lower # of substitutions at the confidence (alpha) interval
            given (int)
        UPPER: upper # of substitutions for the confidence (alpha) interval
            given (int)
    """
    min_time, max_time = time_range
    expected_data = pd.DataFrame(
        {TIME: range(math.floor(min_time), math.ceil(max_time) + 1)})

    def poisson_values(divergence_time):
        """Expected substitutions for time with rate (poisson lambda)."""
        # rate in subs/(site.year), convert to subs/week; multiply by each time
        # point to obtain expected substitutions at that time
        poisson_rate = (rate / 52) * sites * divergence_time
        mode, lower, upper = poisson_mode_and_alpha(poisson_rate, alpha)
        return {MODE: mode, LOWER: lower, UPPER: upper}

    expected_data[[MODE, LOWER, UPPER]] = expected_data[TIME].apply(
        lambda t: pd.Series(poisson_values(t)))
    return expected_data


def make_predictions(pairs_data, relative_ts_to_ca, substitution_rate,
                     num_sites, incubation_period=2, alpha=0.95):
    """
    Uses the Poisson expected values and the pair distances to predict if the
    samples in each pair may be related relative to a putative common ancestor.

    Parameters
    ----------
    pairs_data : pd.DataFrame
        Pandas data frame containing the model predicted results for
        each file.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * DISTANCE: number of differences between each pair's sequences (int)
    relative_ts_to_ca : iter
        An iterable object containing ints or floats, each representing a time
        between the most recent sample and a presumed common ancestor (Dt_pCA).
    substitution_rate : float
        Substitution rate in substitutions/(site.year)
    num_sites : int
        The number of sites in the MSA.
    incubation_period : int or float, default 2
        Number of weeks between infection and contagion of next patient. 2
        weeks for measles.
    alpha : float, default 0.95
        Fraction of the Poisson distribution to fall within interval. Between
        0 and 1.
        e.g., 0.95: 95% of the Poisson distribution at a given time will fall
        within the lower-upper range of number of substitutions.

    Returns
    -------
    pd.DataFrame
        Pandas data frame containing the observed and predicted results for
        each file.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * T_TO_PUTATIVE_CA: weeks between most recent sample and the putative
            common ancestor (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA in weeks (float)
        * PREDICTED: pair category based on model - related or unrelated (str)
    """
    pairs_to_classify = pd.DataFrame()

    for relative_t_to_ca in relative_ts_to_ca:
        info_for_putative_ca = pairs_data.copy()
        info_for_putative_ca[T_TO_PUTATIVE_CA] = relative_t_to_ca

        # relative_t_to_ca is the number of weeks from most recent sample
        # to the putative common ancestor (Dt_pCA)
        # if the difference in time between samples 1 and 2 is longer than the
        # time from the putative common ancestor (CA cannot be later than a
        # descendant), ignore the pair, e.g.:
        #    /--- Sample 1
        #    \--------- Sample 2
        #          |<->| Dt_pCA    would not make sense to consider these pairs
        # also ignore pairs where Sample 1 occurs less than an incubation
        # period (2 weeks for measles) after the putative common ancestor (CA
        # cannot be contemporaneous to descendant):
        #    / Sample 1
        #    \------ Sample 2
        #      |<--->| Dt_pCA    ignore these pCA/pair relationships too
        info_for_putative_ca = info_for_putative_ca[
            (info_for_putative_ca[PAIR_DELTA] + incubation_period)
            <= relative_t_to_ca]
        # calculate cumulative evolution time (Dt_CE)
        # 2 samples may have accumulated differences over the sum of times
        # between each sample and the putative common ancestor (effectively sum
        # of the number of dashes between samples and pCA)
        info_for_putative_ca[PRED_EVOL_T] = (
                2 * relative_t_to_ca - info_for_putative_ca[PAIR_DELTA])

        pairs_to_classify = pd.concat(
            [pairs_to_classify, info_for_putative_ca])

    # calculate expected substitution range for the full range of Dt_CE
    min_tce = pairs_to_classify[PRED_EVOL_T].min()
    max_tce = pairs_to_classify[PRED_EVOL_T].max()
    # data frame columns: Time, Typical, Lower, Upper
    poisson_values = expected_substitutions(
        substitution_rate, num_sites, alpha, (min_tce, max_tce))
    # classify based on Poisson predictions and model Dt_CE
    # add substitution ranges to df
    pairs_to_classify = pairs_to_classify.join(
        poisson_values.set_index(TIME), on=PRED_EVOL_T)
    pairs_to_classify[PREDICTED] = np.nan
    # distance <= upper end of expected substitution range for Dt_CE: H0 True
    # (related)
    pred_related_pairs = (
            pairs_to_classify[DISTANCE] <= pairs_to_classify[UPPER])
    pairs_to_classify.loc[pred_related_pairs, PREDICTED] = RELATED
    # distance > high end of expected substitution range for Dt_CE: H0 False
    # (unrelated)
    pred_unrelated_pairs = (
            pairs_to_classify[DISTANCE] > pairs_to_classify[UPPER])
    pairs_to_classify.loc[pred_unrelated_pairs, PREDICTED] = UNRELATED
    pairs_to_classify[PREDICTED] = (
        pairs_to_classify[PREDICTED].astype('boolean'))
    # drop unnecessary columns
    return pairs_to_classify.drop(columns=[DISTANCE, MODE, LOWER, UPPER])


def make_estimates(pairs_data):
    """
    Classifies the Bayesian-estimated times to the MRCA against the putative
    times to a common ancestor.

    Parameters
    ----------
    pairs_data : pd.DataFrame
        Pandas data frame containing the observed and predicted results for
        each file.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * T_TO_PUTATIVE_CA: time from putative CA to most recent sample (float)
        * T_TO_EST_MRCA: weeks between most recent sample and MRCA (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA (float)
        * PREDICTED: pair category based on model - related or unrelated (str)
        * REGION: results' genomic region (str)
        * T_TO_EST_MRCA: weeks between the most recent sample and the Bayesian-
          estimated MRCA for each pair of samples (float).

    Returns
    -------
    pd.DataFrame
        Pandas data frame containing the observed and predicted results for
        each file.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * T_TO_PUTATIVE_CA: time from putative CA to most recent sample (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA (float)
        * PREDICTED: pair category based on model - related or unrelated (bool)
        * REGION: genomic region (str)
        * ESTIMATED: pair category based on time tree - related or unrelated
            (bool)
    """
    pairs_data[ESTIMATED] = np.nan
    # Dt_bMRCA <= Dt_pCA: related
    # ___ /--- Sample 1
    #     \--------- Sample 2
    #  |    |
    #  |    Dt_pCA2 --> samples diverged before pCA2, classified as unrelated
    #  |                                              (H0 False)
    #  Dt_pCA1 --> pCA1 before pair's divergence, classified related (H0 True)
    est_related_pairs = (
            pairs_data[T_TO_EST_MRCA] <= pairs_data[T_TO_PUTATIVE_CA])
    pairs_data.loc[est_related_pairs, ESTIMATED] = RELATED
    # Dt_bMRCA > Dt_pCA: unrelated (H0 False)
    est_unrelated_pairs = (
            pairs_data[T_TO_EST_MRCA] > pairs_data[T_TO_PUTATIVE_CA])
    pairs_data.loc[est_unrelated_pairs, ESTIMATED] = UNRELATED

    pairs_data[ESTIMATED] = pairs_data[ESTIMATED].astype('boolean')
    return pairs_data.drop(columns=T_TO_EST_MRCA)


def classify_results(pairs_results):
    """
    Classifies the model-predictions against the Bayesian-estimates.

    Parameters
    ----------
    pairs_results : pd.DataFrame
        Pandas data frame containing the observed and predicted results for
        each file.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * T_TO_PUTATIVE_CA: time from putative CA to most recent sample (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA (float)
        * PREDICTED: pair category based on model - related (True) or unrelated
            (False) (bool)
        * REGION: genomic region (str)
        * ESTIMATED: pair category based on time tree - related (True) or
            unrelated (False) (bool)

    Returns
    -------
    pd.DataFrame
        Pandas data frame containing the pair results classification.
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * REGION: genomic region (str)
        * T_TO_PUTATIVE_CA: time from putative CA to most recent sample (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA (float)
        * RESULT: result classification based on estimated and predicted (str)
    """
    estimated_negative = pairs_results[ESTIMATED] == RELATED
    estimated_positive = pairs_results[ESTIMATED] == UNRELATED
    predicted_negative = pairs_results[PREDICTED] == RELATED
    predicted_positive = pairs_results[PREDICTED] == UNRELATED

    true_negatives = estimated_negative & predicted_negative
    true_positives = estimated_positive & predicted_positive
    false_negatives = estimated_positive & predicted_negative
    false_positives = estimated_negative & predicted_positive

    pairs_results[RESULT] = np.nan
    pairs_results.loc[true_negatives, RESULT] = TRUE_NEG
    pairs_results.loc[true_positives, RESULT] = TRUE_POS
    pairs_results.loc[false_negatives, RESULT] = FALSE_NEG
    pairs_results.loc[false_positives, RESULT] = FALSE_POS

    cat_type = pd.CategoricalDtype(
        categories=[TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS],
        ordered=True)
    pairs_results[RESULT] = pairs_results[RESULT].astype(cat_type)

    pt = pairs_results[[REGION, PAIR_DELTA, T_TO_PUTATIVE_CA, PRED_EVOL_T,
                        RESULT]].pivot_table(
        index=[REGION, PAIR_DELTA, T_TO_PUTATIVE_CA, PRED_EVOL_T],
        columns=RESULT,
        aggfunc=len, fill_value=0)
    if pt.shape[1] < cat_type.categories.shape[0]:
        for cat in cat_type.categories:
            if cat not in pt.columns:
                pt[cat] = 0
    return pt


def results_for_tree(tree, pairs_index, pairs_predictions):
    """
    Given a dendropy tree object, a pandas 2-level index and model predictions
    for each pair, classifies model predictions against the bayesian-estimated
    tree results.

    Parameters
    ----------
    tree : ddp.datamodel.treemodel.Tree
        A dendropy tree object containing a time-scaled Bayesian-estimated
        phylogenetic tree.
    pairs_index : pd.Index
        Pandas 2-level index of sample pair identifiers (str, str).
    pairs_predictions : pd.DataFrame
        * 2-level index with sample identifiers for each sample pair
        * PAIR_DELTA: weeks difference between the pair of samples (int)
        * T_TO_PUTATIVE_CA: weeks between most recent sample and the putative
            common ancestor (float)
        * PRED_EVOL_T: cumulative evolution time (Dt_CE) for the pair since
            the putative CA in weeks (float)
        * PREDICTED: pair category based on model - related or unrelated (str)
        * REGION: genomic region of prediction (str)

    Returns
    -------
    pd.DataFrame
    """
    # get Dt_bMRCA against the pair index
    pairs_tree_data = pairs_estimated_times_to_mrca(tree, pairs_index)
    # classify pairs as related/unrelated
    pairs_tree_data = make_estimates(pairs_predictions.join(pairs_tree_data))
    # compare estimates against predictions (TN/FN/TP/FP)
    return classify_results(pairs_tree_data)


def setup_pairs_data_df(samples_list):
    """
    Given a list of sample identifiers returns a data frame with pairs of
    samples as index and a PAIR_DELTA column containing the time between sample
    collection for each pair in the index.

    Parameters
    ----------
    samples_list : iter
        An iterable object containing a list of sample identifiers in the
        format 'H123456789_2000.5' (str).
    Returns
    -------
    pd.DataFrame
        * 2-level index: pair's sample ids (str)
        * PAIR_DELTA: weeks difference between the pair of samples (int)
    """
    return calculate_pairs_delta(
        pairs_index_from_iter(samples_list),
        dates_from_fas_identifiers(samples_list))


def results_for_dataset(
        trees_path, fasta_paths, ts_to_putative_ca, rates_csv,
        trees_path_regex=file_names_style['style3'],
        fasta_paths_regex=file_names_style['style3'], incubation_period=2,
        probability_interval=0.95, debug=False):
    """
    Makes model predictions of sample pair relatedness for a set of fasta files
    and putative common ancestors and compares them against the results from
    a Bayesian-estimated set of trees. Results for each genotype/dataset are
    saved as csv files with columns (names set as global variables):
        TREE_NUMBER: number of results' tree (int)
        REGION: results' genomic region (str)
        PAIR_DELTA: weeks between the samples in pair, Dt (float)
        T_TO_PUTATIVE_CA: weeks to a putative common ancestor, Dt_pCA (int)
        PRED_EVOL_T: Dt_CE, 2* Dt_pCA - Dt for each row (float)
        TRUE_NEG: number of sample pairs for which H0 is accepted using both
            model and BEAST tree at the given Dt/Dt_pCA combination (int)
        TRUE_POS: number of sample pairs for which H0 is rejected using both
            model and BEAST tree at the given Dt/Dt_pCA combination (int)
        FALSE_NEG: number of sample pairs for which H0 is rejected using the
            BEAST tree and accepted using the model at the given Dt/Dt_pCA
            combination (int)
        FALSE_POS: number of sample pairs for which H0 is accepted using the
            BEAST tree and rejected using the model at the given Dt/Dt_pCA
            combination (int)

    Parameters
    ----------
    trees_path : str
        Path to a file containing a collection of nexus Bayesian-estimated
        trees where branch length is set in years.
    fasta_paths : list or tuple or set
        List containing fasta files with MSAs for the regions used in
        calculating the trees.
    ts_to_putative_ca : list or tuple or set
        A list of putative times (ints) between the most recent sample in each
        pair and a putative common ancestor in weeks.
    rates_csv : str
        The path to a csv file with the rates for the genotype and genomic
        regions of interest.
    trees_path_regex : str
        A regular expression pattern with named groups for genotype, dataset
        and BEAST settings (set in global variables at top of file).
    fasta_paths_regex : str
        A regular expression pattern with named groups for genomic region (set
        in global variables at top of file).
    incubation_period : int or float, default 2
        Number of weeks between infection and contagion of next patient. 2
        weeks for measles.
    probability_interval : float, default 0.95
        Fraction of the Poisson distribution to fall within interval. Between
        0 and 1.
        e.g., 0.95: 95% of the Poisson distribution at a given time will fall
        within the lower-upper range of number of substitutions.
    debug : bool, default False
        When debugging, only process first 5 trees.

    Returns
    -------
    None
    """
    start_time = datetime.now()
    # Setup
    # -----
    # get genotype/dataset info from trees path and load rates
    trees_path_info = re.search(trees_path_regex, trees_path).groupdict()
    genotype = trees_path_info[REGEX_GENOTYPE]
    dataset = trees_path_info[REGEX_DATASET]
    settings = trees_path_info[REGEX_SETTINGS]
    rates = pd.read_csv(rates_csv,
                        index_col=[RATES_GENOTYPE, RATES_REGION],
                        dtype={SITES: int})

    # Process fastas
    # --------------
    # get sample list and pairs index from first fas processed
    pairs_delta = None
    pairs_data = pd.DataFrame()

    for fasta_file in fasta_paths:
        records_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
        pairs_delta = (pairs_delta if pairs_delta is not None
                       else setup_pairs_data_df(records_dict.keys()))
        region_data = pairs_delta.copy()

        # get distances for region
        region_data = sample_pairs_distances(records_dict, region_data)

        # get region info from fasta path
        region = re.search(
            fasta_paths_regex, fasta_file).group(REGEX_REGION).replace(
            'N450', 'N-450').replace('MFNCR', 'MF-NCR')

        # Model
        # -----
        # make predictions for the region
        # rate and number of sites in MSA from rates.csv
        rate = rates.loc[(genotype, region), RATE]
        num_sites = rates.loc[(genotype, region), SITES]
        # predictions against putative CAs
        region_data = make_predictions(
            region_data, ts_to_putative_ca, rate, num_sites,
            incubation_period=incubation_period, alpha=probability_interval)

        region_data[REGION] = region
        # concat to main df
        pairs_data = pd.concat([pairs_data, region_data])

        print(f'Made model predictions for region {region} for file '
              f'{trees_path}')

    region_cats = pd.CategoricalDtype(
        categories=['N-450', 'MF-NCR'], ordered=True)
    pairs_data[REGION] = pairs_data[REGION].astype(region_cats)

    # Process trees
    # -------------
    pairs_index = pairs_delta.index
    print(f'Completed model predictions for {trees_path} in '
          f'{datetime.now() - start_time}. Started processing trees.')
    in_dir, _ = os.path.split(trees_path)
    out_path = os.path.join(
        in_dir.replace('input', 'output'),
        f'{genotype}-{dataset}-{settings}-alpha'
        f'{int(probability_interval * 100)}-trees_summary.csv')
    # load each tree
    for tree_idx, tree in enumerate(ddp.Tree.yield_from_files(
            files=[trees_path], schema='nexus',
            extract_comment_metadata=True, store_tree_weights=True,
            rooting='default-rooted')):
        tree_results = results_for_tree(tree, pairs_index, pairs_data)
        tree_results = pd.concat({tree_idx: tree_results}, names=[TREE_NUMBER])

        write_col_names = False if os.path.exists(out_path) else True
        tree_results.to_csv(out_path, mode='a', header=write_col_names)
        if tree_idx % 100 == 0:
            print(f'Processed {tree_idx} trees from {trees_path} after '
                  f'{datetime.now() - start_time}.')
        if debug:
            if tree_idx > 5:
                break
    print(f'Finished processing {trees_path}. Total time: '
          f'{datetime.now() - start_time}.')
