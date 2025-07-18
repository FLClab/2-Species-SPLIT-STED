from ast import arg
import numpy as np
    

def ExpFunMono_MLE(p,args): 
#[err, z] = ExpFunMono(p, t, y, mle), minimised, faster version of ExpFunFull
# p   - parameters (bg, , tau1 amp1 tau2 amp2 ...)
# args: t   - time (x) axis, does not have to start at zero but has to be unifrom
# args: y   - function values (y-axis)
    
# err - error: negative log likelihood

    
    # Christoph Thiele, 2020
    t = args[0]
    p = p
    y = args[1]
    n = np.asarray(t).size
    dt = (t[-1] - t[0]) / n
    
    z = np.exp(- (t) / p[1])
    # # z = p(3)*z/sum(z) + p(0)/n;     # Norm numerically
    z = p[2] / p[1] * dt * z + p[0] / n # Norm analytically
    
    # MLE
    ind = y > 0
    if np.all(ind):
        err = sum(np.multiply(- y,np.log(z + 1e-30)) + z)
    else:
        err = sum(np.multiply(- y[ind],np.log(z[ind] + 1e-30))) + sum(z)

    return err

