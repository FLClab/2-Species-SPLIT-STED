import numpy as np
""" 
Convertion out of the box du script Matlab
def ExpFunLSQ(p = None,t = None,y = None,mle = None): 
    #[err, z] = ExpFunLSQ(p, t, y, mle), minimised, faster version of ExpFun
# Does a nonnegative least-square fit of y with the lifetimes p.
# p   - parameters (tau1 tau2 ...)
# t   - time (x) axis, does not have to start at zero but has to be unifrom
# y   - function values (y-axis)
# mle - if true the MLE (negative log likelihood) is returned, otherwise chi^2
    
    # err - error: negative log likelihood OR chi^2
# c   - compontens amplitudes (bg amp1 amp2 ...)
# zz  - compontens function value
# z   - total function value
    
    # Christoph Thiele, 2020
    
    t = t
    p = np.transpose(p)
    y = y
    ind = y > 0
    n = np.asarray(t).size
    dt = (t(end()) - t(1)) / n
    
    zz = np.array([np.ones((y.shape,y.shape)) / n,np.exp(- t / p) / p * dt])
    c = lsqnonneg(zz,y)
    z = zz * c
    if mle:
        # MLE
        if np.all(ind):
            err = sum(np.multiply(- y,np.log(z + 1e-30)) + z)
        else:
            err = sum(np.multiply(- y(ind),np.log(z(ind) + 1e-30))) + sum(z)
    else:
        # chi2
        err = sum((y - z) ** 2.0 / np.abs(z)) / np.asarray(t).size
    
    return err,c,zz,z
""" 
# Traduit et adapté du code Matlab
def ExpFun(p0 = None,alphas = None,t = None,y = None,mle = None): 
    """Compute the MLE & Chi2 error of multi-exponential

    Parameters
    ----------
    p0 : Numpy array
        Array of lifetimes Tau
    alphas : Numpy array
        Array of concentration for each liftimes
    t : Numpy array
        Temporal array
    y : Numpy array
        Array of measured data
    mle : Bool
        MLE or Chi2

    Returns
    -------
    Float
        Error
    From Christoph Thiele, 2020
    """
    t = t
    #p = np.transpose(p)

    y = y
    ind = y > 0
    n = np.asarray(t).size
    dt = (t[-1] - t[1]) / n

    z = np.zeros_like(t)
    for tau, alpha in zip(p0, alphas):
        z = alpha * np.exp(- t / tau) / tau + z # Norm analytically
    z = z * dt

    if mle:
        # MLE
        if np.all(ind):
            err = sum(np.multiply(- y,np.log(z + 1e-30)) + z)
        else:
            err = sum(np.multiply(- y[ind],np.log(z[ind] + 1e-30))) + sum(z)
    else:
        # chi2
        err = sum((y - z) ** 2.0 / np.abs(z)) / np.asarray(t).size
    return err


# Fonction un peu modifié pour fonctionner avec minize de SCYPI


def ExpFun_bi_MLE_tau_and_alpha(p0,args):  
    """Compute the MLE error of bi-exponential over taus and alphas

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
    n = np.asarray(t).size
    dt = (t[-1] - t[1]) / n

    z = np.zeros_like(t)
    for tau, alpha in zip(p, alphas):
        z = alpha * np.exp(- t / tau) / tau + z # Norm analytically
    z = z * dt

    if np.all(ind):
        err = sum(np.multiply(- y,np.log(z + 1e-30)) + z)
    else:
        err = sum(np.multiply(- y[ind],np.log(z[ind] + 1e-30))) + sum(z)


    return err


def ExpFun_bi_MLE_tau(p0,args):  
    """Compute the MLE error of bi-exponential over taus and alphas

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
    p = [args[0], args[1]]
    t = args[2]
    y = args[3]
    alphas = [p0[0], 1-p0[0]]
    ind = y > 0
    n = np.asarray(t).size
    dt = (t[-1] - t[1]) / n

    z = np.zeros_like(t)
    for tau, alpha in zip(p, alphas):
        z = alpha * np.exp(- t / tau) / tau + z # Norm analytically
    z = z * dt

    if np.all(ind):
        err = sum(np.multiply(- y,np.log(z + 1e-30)) + z)
    else:
        err = sum(np.multiply(- y[ind],np.log(z[ind] + 1e-30))) + sum(z)


    return err