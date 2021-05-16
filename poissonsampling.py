import math
from scipy.stats import poisson


def probability_for_n_changes(n, time, mu):
    """
    Probability that exactly `n` substitutions have been experienced in `time`
    for region with substitution rate `mu` substitutions/time unit.

    (from Bayesian evolutionary analysis with BEAST Drummond and Boukaert,
    4.2 The molecular clock, p. 58; code by DW, adapted by AP)

    Parameters
    ----------
    n : int
        Number of substitutions.
    time : int
        Time t elapsed during which changes can occur.
    mu : float
        Substitutions per time. Multiply substitution rate in
        substitutions/(site*time) by number of sites in the relevant region.

    Returns
    -------
    float
        Poisson probability.
    """
    return (math.exp(- mu * time) * (mu * time) ** n) / math.factorial(n)


def poisson_mode_and_alpha(expected, alpha):
    """
    Returns mode number of expected substitutions and upper & lower limits
    within which falls `alpha` fraction of the occurrences for a poisson
    distribution with lambda=`expected` value.

    Parameters
    ----------
    expected : float
        Expected substitutions.
    alpha : float
        Probability interval containing alpha fraction of the Poisson
        distribution.

    Returns
    -------
    mode : int
        Typical number of substitutions expected.
    lower_limit, upper_limit : int
        Lower number of substitutions between which alpha fraction of the
        Poisson distribution is contained.
    upper_limit : int
        Lower number of substitutions between which alpha fraction of the
        Poisson distribution is contained.
    """
    mean, var = poisson.stats(expected, moments='mv')
    mode = math.floor(mean)  # round down to closest integer
    lower_limit, upper_limit = poisson.interval(alpha, expected, loc=0)

    return mode, lower_limit, upper_limit
