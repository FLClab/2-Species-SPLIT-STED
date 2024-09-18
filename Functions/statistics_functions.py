
import numpy

from matplotlib import pyplot
from scipy import stats
import scikit_posthocs

import ANOVA


def verify_normality(samples, alpha=0.05):
    """
    Verifies the normality of the data
    """
    # print([stats.shapiro(sample)[1] for sample in samples])
    return [stats.shapiro(sample)[1] > alpha for sample in samples]


def get_significance(samples, show_qq=False, force_normal_test=False, show_sr=False, verbose=False):
    """
    Computes a statistic analysis of the given samples.

    :param samples: A list of sample observations
    :param show_qq: Wheter to show the quantiles
    :param force_normal_test: Wheter to force the normal test
    :param show_sr: Wheter to show standard residuals
    :param verbose: Wheter to use verbose
    """
    _print = print if verbose else lambda *args, **kwargs : None # Defines the level of verbosity
    fitted, residuals, standard_res, quantiles = ANOVA.fit_model(samples)
    if show_qq:
        fig, ax = pyplot.subplots()
        for params, fit in quantiles:
            ax.scatter(params[0], params[1], facecolor="white", edgecolor="black")
        lims = [
            numpy.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            numpy.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        ax.plot(lims, lims, color="gray", linestyle="dashed")
        ax.set_title("QQ-Plot")
    if show_sr:
        fig, ax = pyplot.subplots()
        ax.boxplot(standard_res)

    # print([fit[2] for params, fit in quantiles])
    # normality = all([fit[2] > 0.9 for params, fit in quantiles])
    normality = all(verify_normality(standard_res))
    if normality or force_normal_test :
        _print("The standardize residuals are normaly distributed.")
        statistic, pvalue = stats.f_oneway(*samples)
        if pvalue > 0.05:
            _print("The one-way ANOVA test that two or more groups have the same " +
            "population mean cannot be rejected with alpha confidence of 0.05.\n" +
            "\tpvalue : {}".format(pvalue))
            pvalues = numpy.ones((len(samples), len(samples)))
        else:

            _print("The null hypothesis can be rejected from the ANOVA one way test.\n"+
                    "\tpvalue : {}".format(pvalue))
            _print("Posthoc ttest is computed on the data")
            pvalues = scikit_posthocs.posthoc_ttest(samples)
            _print("Posthoc t-test results\n" +
                    "\tpvalues : {}".format(pvalues))
            pvalues=pvalues.values
            
    else :
        _print("The standardize residuals are not normaly distributed.")
        statistics, pvalue = stats.kruskal(*samples)
        if pvalue > 0.05:
            _print("The null hypothesis cannot be rejected from the Kruskal-Wallis H-test.")
            pvalues = numpy.ones((len(samples), len(samples)))
        else:
            _print("The null hypothesis can be rejected from the Kruskal-Wallis H-test.\n"+
                    "\tpvalue : {}".format(pvalue))
            _print("Posthoc Dunn is computed on the data.")
            pvalues = scikit_posthocs.posthoc_dunn(samples)
            _print("Posthoc Dunn results\n" +
                    "\tpvalues\n : {}".format(pvalues))
            pvalues=pvalues.values

    combinations = [(i, j) for i in range(len(samples)) for j in range(i + 1, len(samples))]
    significance = []
    for (i, j) in combinations:
       
        pvalue = pvalues[i,j]
    

        if pvalue <= 0.001:
            significance.append((pvalue, 3, i, j))
        elif pvalue <= 0.01:
            significance.append((pvalue, 2, i, j))
        elif pvalue <= 0.05:
            significance.append((pvalue, 1, i, j))
        else:
            significance.append((pvalue, 0, i, j))
    return significance
