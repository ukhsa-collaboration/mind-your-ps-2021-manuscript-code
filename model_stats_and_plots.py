# std python libraries
from collections import namedtuple
import math
import os
import re
# 3rd party libraries
from matplotlib import rcParams
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# my modules
from info_parsing_dicts import file_names_style
from model_validation import (
    TREE_NUMBER, REGION, PAIR_DELTA, T_TO_PUTATIVE_CA, PRED_EVOL_T, TRUE_NEG,
    TRUE_POS, FALSE_NEG, FALSE_POS, GENOTYPE, REGEX_GENOTYPE, REGEX_DATASET,
    REGEX_SETTINGS)


# Global vars
# ===========
LEGEND_LABEL = 'Rate'
Y_LABEL = 'rate (%)'
ACCURACY = 'ACC'
FALSE_DETECTION = 'FDR'
FALSE_OMISSION = 'FOR'
NEG_PRED_VAL = 'NPV'
POS_PRED_VAL = 'PPV'
SERIES_COLOURS = {
    FALSE_NEG: '#92c5de', TRUE_NEG: '#0571b0', FALSE_POS: '#f4a582',
    TRUE_POS: '#ca0020', ACCURACY: '#3969AC', FALSE_DETECTION: '#F2B701',
    FALSE_OMISSION: '#E73F74', NEG_PRED_VAL: '#7F3C8D', POS_PRED_VAL: '#11A579'
}

CALCULATE_RATES = [
    TRUE_NEG, FALSE_NEG, TRUE_POS, FALSE_POS, ACCURACY, FALSE_DETECTION,
    FALSE_OMISSION, NEG_PRED_VAL, POS_PRED_VAL]
Rate = namedtuple('Rate', ['numerator', 'denominator'])
tnr = Rate(numerator=(TRUE_NEG,), denominator=(TRUE_NEG, FALSE_POS))
tpr = Rate(numerator=(TRUE_POS,), denominator=(TRUE_POS, FALSE_NEG))
fnr = Rate(numerator=(FALSE_NEG,), denominator=(TRUE_POS, FALSE_NEG))
fpr = Rate(numerator=(FALSE_POS,), denominator=(TRUE_NEG, FALSE_POS))
acc = Rate(numerator=(TRUE_NEG, TRUE_POS),
           denominator=(TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS))
fdr = Rate(numerator=(FALSE_POS,), denominator=(FALSE_POS, TRUE_POS))
fom = Rate(numerator=(FALSE_NEG,), denominator=(FALSE_NEG, TRUE_NEG))
npv = Rate(numerator=(TRUE_NEG,), denominator=(TRUE_NEG, FALSE_NEG))
ppv = Rate(numerator=(TRUE_POS,), denominator=(TRUE_POS, FALSE_POS))
RATES_DEFS = {
    TRUE_NEG: tnr, TRUE_POS: tpr, FALSE_NEG: fnr, FALSE_POS: fpr,
    ACCURACY: acc, FALSE_DETECTION: fdr, FALSE_OMISSION: fom,
    NEG_PRED_VAL: npv, POS_PRED_VAL: ppv
}

X_LABELS = {
    PRED_EVOL_T: r'$\mathbf{\Delta t_{CE}}$ (weeks)',
    PAIR_DELTA: r'$\mathbf{\Delta t}$ (weeks)',
    T_TO_PUTATIVE_CA: r'$\mathbf{\Delta t_{pCA}}$ (weeks)'
}


# for editable text within Inkscape
rcParams['svg.fonttype'] = 'none'


def calculate_rate(data, rate):
    """
    Given a data frame with model vs. estimate results, groups results and
    calculates a rate obtained by dividing the sum of the columns matching the
    numerators by the sum of the columns corresponding to the denominators.

    Parameters
    ----------
    data : pd.DataFrame
        Pandas data frame containing at least the columns needed to calculate
        the rate.
    rate : str
        The name of the rate to be calculated as in RATES_DEFS global dict.

    Returns
    -------
    pd.Series
        Pandas series with the calculated rate.
    """
    if rate not in RATES_DEFS:
        raise KeyError(f'Rate {rate} not defined.')
    numerator = sum(data[col] for col in RATES_DEFS[rate].numerator)
    denominator = sum(data[col] for col in RATES_DEFS[rate].denominator)

    return pd.Series(100 * numerator / denominator)


def prepare_df_for_plotting(data, group_by=PRED_EVOL_T):
    """
    Given a data frame in the format returned by the validation script,
    calculates model statistics and prepares data frame for plotting.

    Parameters
    ----------
    data : pd.DataFrame
        Data frame as produced by the validation script, with columns:
        * TREE_NUMBER, number of the tree for which the estimates were made;
        * REGION, genomic region to which results apply;
        * PAIR_DELTA, time difference between sample pairs classified;
        * T_TO_PUTATIVE_CA, time elapsed between most recent sample in pair and
          a putative common ancestor of the pair;
        * PRED_EVOL_T, predicted time that the samples in the pairs may have
          diverged over since a putative common ancestor;
        * TN/TP/FN/FP, total number of pairs classified as true/false
          negative/positive.
    group_by : str, default PRED_EVOL_T
        Column to group data by. By default, PRED_EVOL_T, but could also be
        PAIR_DELTA or T_TO_PUTATIVE_CA. The remaining time parameter columns
        are dropped and the total of pairs classified in each category is
        calculated for each tree/region/time combination.
        'group_by' column.
    calculate_rates : tuple, default (TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS)
        Tuple containing name of rates to be calculated for dataset.

    Returns
    -------
    pd.DataFrame
        Data frame in long form containing rates for plotting.
    """
    # keep only relevant columns
    keep_columns = [TREE_NUMBER, REGION, group_by,
                    TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS]
    for column in keep_columns:
        if column not in data.columns:
            raise KeyError(f'Column {column} not in the data frame to plot.')
    relevant_data = data[keep_columns].reset_index(drop=True)
    min_aggregate_value = math.floor(relevant_data[group_by].min())
    max_aggregate_value = math.ceil(relevant_data[group_by].max())
    if group_by in [PAIR_DELTA, PRED_EVOL_T]:
        bins = list(
            range(max(0, min_aggregate_value - 4), max_aggregate_value + 4, 4))
        relevant_data[group_by] = pd.cut(
            relevant_data[group_by],
            bins=bins, labels=bins[1:], include_lowest=True)
    # pivot with sum of pairs classified as TN/TP/FN/FP
    group_pivot = relevant_data.pivot_table(
        index=[REGION, TREE_NUMBER, group_by],
        values=[TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS],
        aggfunc=sum)

    # calculate rates
    group_percents = group_pivot.copy()
    for rate in CALCULATE_RATES:
        group_percents[rate] = calculate_rate(group_pivot, rate)

    group_percents = (
        group_percents
            .reset_index()
            .melt(id_vars=[REGION, TREE_NUMBER, group_by],
                  var_name=LEGEND_LABEL, value_name=Y_LABEL))
    rate_cats = pd.CategoricalDtype(categories=CALCULATE_RATES, ordered=True)
    group_percents[LEGEND_LABEL] = (group_percents[LEGEND_LABEL]
                                    .astype(rate_cats))
    return group_percents


def calculate_dataset_rates(
        files_to_process, file_name_style=file_names_style['style4'],
        group_by=PRED_EVOL_T,
        region_order=('N-450', 'MF-NCR', 'N-450+MF-NCR')):
    """
    Process set of files containing tree summaries to produce a table
    containing the statistics for each data set.

    Parameters
    ----------
    files_to_process : tuple
        A tuple containing tuples, each containing paths of csv files with tree
        summaries for each genotype to be processed together.
    file_name_style : regex
        A string containing a regular expression to extract information from
        the csv file name.
    group_by : str
        The column heading to be used to group results by.
    region_order : tuple
        A tuple containing the list of genomic regions in the order to be
        plotted.

    Returns
    -------
    None
    """
    for file_set in files_to_process:
        all_rates = pd.DataFrame()
        out_path = None
        for file in file_set:
            # load results
            rates = pd.read_csv(file)
            region_cats = pd.CategoricalDtype(
                categories=region_order, ordered=True)
            rates[REGION] = rates[REGION].astype(region_cats)
            # get genotype from file name
            in_dir, _ = os.path.split(file)
            name_info = re.search(file_name_style, file).groupdict()
            genotype = name_info[REGEX_GENOTYPE]
            if out_path is None:
                # set output handle
                dataset = name_info[REGEX_DATASET]
                settings = name_info[REGEX_SETTINGS]
                out_path = os.path.join(
                    in_dir.replace('input', 'output'),
                    f'{dataset}-{settings}-{group_by}')
            # calculate rates and set df to long format for plotting
            rates = prepare_df_for_plotting(rates, group_by=group_by)
            rates[GENOTYPE] = genotype
            all_rates = pd.concat([all_rates, rates])
        stats_out_file = f'{out_path}-stats.csv'
        all_rates.to_csv(stats_out_file, index=False)
        files_done = '\n'.join(file_set)
        print(f'Completed processing set of files containing\n{files_done}\n'
              f'Model statistics table saved to {stats_out_file}\n')


def get_axes_numbers(figure):
    """
    Given a figure with one or more subplots, returns a dict of information on
    the position and number of subplots in the figure.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        A matplotlib figure object containing one or multiple axes.

    Returns
    -------
    dict
        * 'x_lims': (x_min, x_max) start and finish of all subplots as
            proportion of the width of the figure
        * 'y_lims': (y_min, y_max) start and finish of all subplots as
            proportion of the height of the figure
        * 'axes_centres': (x_centre, y_centre) midpoint of all subplots as
            proportion of the width/height of the figure
        * 'num_plots': number of subplots in the figure (int)
    """
    # positions of the first and last subplots
    fig_axes = figure.get_axes()
    top_left_subplot_pos = fig_axes[0].get_position()
    bottom_right_subplot_pos = fig_axes[-1].get_position()

    axes_numbers = {}
    # x and y limits of the subplots
    x_min = top_left_subplot_pos.xmin
    x_max = bottom_right_subplot_pos.xmax
    axes_numbers['x_lims'] = (x_min, x_max)
    y_max = top_left_subplot_pos.ymax
    y_min = bottom_right_subplot_pos.ymin
    axes_numbers['y_lims'] = (y_min, y_max)
    # get mid points of all x and y axes
    x_centre = x_min + (x_max - x_min) / 2
    y_centre = y_min + (y_max - y_min) / 2
    axes_numbers['axes_centres'] = (x_centre, y_centre)
    # axes info
    axes_numbers['num_plots'] = len(fig_axes)
    return axes_numbers


def plot_model_stats(
        model_stats_file, x_axis,
        rates_to_plot=(TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS),
        region_order=('N-450', 'MF-NCR', 'N-450+MF-NCR'), sns_style='ticks',
        sns_context='paper', img_formats=('.png', '.svg', '.pdf'),
        add_hline=5):
    """
    Given data frame containing true/false positive/negative rates for each
    genotype, region, and relative t to a CA, produces images with facet grid
    plots in the formats required.

    Parameters
    ----------
    model_stats_file : str
        Path to csv file containing a table with model statistics data.
        Genotype, Dataset, Region: the set of data for which the results were
            calculated (str)
        DtCA: the time between the most recent sample and a presumed common
            ancestor (int)
        Rates: the type of rate that the 'Rate' column value refers to (str)
        Rate: true/false positive/negative rate percentages (float)
    x_axis : str
        Heading of the pandas data frame column where the values for the
        independent variable are found.
    rates_to_plot : tuple, default (TRUE_NEG, TRUE_POS, FALSE_NEG, FALSE_POS)
        A tuple containing the rates to produce plots for.
    region_order : tuple, default ('N-450', 'MF-NCR', 'N-450+MF-NCR')
        Order in which to plot the regions.
    sns_style : str, default 'ticks'
        Argument passed to seaborn set_style method, which adjusts plot
        settings. Any permitted bt sns.set_style.
    sns_context : str or dict, default 'paper'
        Argument passed to seaborn set_context method, which adjusts plot
        settings. Any permitted by sns.set_context.
    img_formats : tuple, default ('.png', '.svg')
        Tuple containing figure formats to save matplotlib plot as. Any file
        type supported by matplotlib.
    add_hline : int or None, default 5
        Add a horizontal line to each plot at y=add_hline. No line added if
        None.

    Returns
    -------
    None
    """
    # load relevant lines of data frame
    model_stats = pd.read_csv(model_stats_file)
    model_stats = (model_stats[model_stats[LEGEND_LABEL]
                   .isin(rates_to_plot)]
                   .reset_index())

    # set up grid args and formatting for legend for different number of
    # regions/genotypes to be plotted
    grid_args = dict(
        row=GENOTYPE, row_order=sorted(model_stats[GENOTYPE].unique()),
        col=REGION,
        col_order=[col 
                   for col in region_order
                   if col in model_stats[REGION].unique()],
        hue=LEGEND_LABEL, hue_order=rates_to_plot,
        palette=SERIES_COLOURS,
        height=4, aspect=1.5)

    num_genotypes = model_stats[GENOTYPE].nunique()
    num_regions = model_stats[REGION].nunique()

    # set up plot
    sns.set_style(sns_style)
    assert isinstance(sns_context, (str, dict)), \
        f'sns_context must be either str or dict, given {type(sns_context)}'
    sns.set_context(**({'context': sns_context} if isinstance(sns_context, str)
                       else sns_context))
    grid = sns.relplot(data=model_stats, x=x_axis, y=Y_LABEL, kind='line',
                       **grid_args)
    if isinstance(add_hline, int):
        grid.map(plt.axhline, y=add_hline, linestyle='--', color='lightgrey')
    grid.set_axis_labels(' ', ' ')  # white space so labels can fit
    [ax_title.set_title('') for ax_title in grid.axes_dict.values()]
    grid.set(xlim=(0, None), ylim=(0, 100))

    # remove side legend
    handles, labels = grid.axes.flat[0].get_legend_handles_labels()
    grid._legend.remove()

    grid.fig.tight_layout()
    axes_numbers = get_axes_numbers(grid.fig)
    x_min, x_max = axes_numbers['x_lims']
    y_min, y_max = axes_numbers['y_lims']
    x_centre, y_centre = axes_numbers['axes_centres']

    legend = plt.legend(
        handles=handles, labels=labels,
        bbox_to_anchor=(x_centre, 0), loc='upper center',
        ncol=len(rates_to_plot), frameon=False,
        bbox_transform=grid.fig.transFigure)

    # label axes
    grid.fig.text(x_centre, 0.2 * y_min,
                  s=X_LABELS[x_axis], weight='bold', ha='center', va='bottom')
    grid.fig.text(0.2 * x_min, y_centre, s=Y_LABEL, weight='bold',
                  rotation='vertical', ha='left', va='center')

    # label rows and columns
    h_gap = (x_max - x_min) / num_regions
    first_col_label_pos = x_min + h_gap / 2
    subplot_labels = []
    for i, region in enumerate(grid.col_names):
        subplot_labels.append(
            grid.fig.text(
                first_col_label_pos + i * h_gap, y_max,
                s=region, ha='center', va='bottom'))
    v_gap = (y_max - y_min) / num_genotypes
    first_row_label_pos = y_max - v_gap / 2
    for i, genotype in enumerate(grid.row_names):
        subplot_labels.append(
            grid.fig.text(
                x_max, first_row_label_pos - i * v_gap,
                s=genotype, rotation=-90, ha='left', va='center'))

    # save plots and close
    out_handle = (
        f'{os.path.splitext(model_stats_file)[0]}-{rates_to_plot[0]}'
            .replace('stats', 'plots')
            .replace('input', 'output'))
    for img_format in img_formats:
        grid.fig.savefig(f'{out_handle}{img_format}',
                         bbox_extra_artists=(*subplot_labels, legend),
                         bbox_inches='tight', dpi=300)
    plt.close()
    print(f'Model statistics plots saved with handle {out_handle}')
