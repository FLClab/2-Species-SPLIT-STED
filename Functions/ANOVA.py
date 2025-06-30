
import numpy
from scipy import stats


def fit_model(data):
    """
    This function implements the linear model for the one-way ANOVA.

    :param data: A list of list of each observations in each treatments

    :returns : The fitted model, the residuals, the standardize residuals and
               the quantiles for probability plot

    Ref : https://onlinelibrary.wiley.com/doi/full/10.1111/jac.12220
    """
    fitted, residuals, standard_residuals, quantiles = [],[],[],[]
    for treatment in data:
        treatment = numpy.array(treatment)
        mean = numpy.mean(treatment)
        fitted.append(mean)
        residuals.append(treatment - mean)
        sr = (treatment - mean) / numpy.std(treatment)
        standard_residuals.append(sr)
        quantiles.append(stats.probplot(sr, dist="norm", fit=True))
    return fitted, residuals, standard_residuals, quantiles
