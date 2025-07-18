import numpy 

def ExpFun_bi_MLE_tau_and_alpha(p0,args):  
    """Compute the MLE error of bi-exponential over taus and alphas (lifetimes and amplitudes).

    Parameters
    ----------
    p0 : Initial Guess
        [tau1, tau2, alpha1, alpha2]
    args : Argument
        [Temporal array, data]

    Returns
    -------
    Float
        MLE error
    """
    
    t = args[0]
    y = args[1]
    p = [p0[0], p0[1]]
    alphas = [p0[2], p0[3]]
    ind = y > 0
    n = numpy.asarray(t).size
    dt = (t[-1] - t[1]) / n

    z = numpy.zeros_like(t)
    for tau, alpha in zip(p, alphas):
        z = alpha * numpy.exp(- t / tau) / tau + z # Norm analytically
    z = z * dt

    if numpy.all(ind):
        err = sum(numpy.multiply(- y,numpy.log(z + 1e-30)) + z)
    else:
        err = sum(numpy.multiply(- y[ind],numpy.log(z[ind] + 1e-30))) + sum(z)


    return err

