import skimage.io as skio
import numpy as np
import pandas as pd

from inspect import currentframe, getframeinfo
import warnings




def decorr_res(imname=None, image=None):
    #Input image and analysis parameter
    Nr = 50
    Ng = 10
    w = 20
    find_max_tol = 0.0005 #In the suppl info of the article it was 0.001...
    if image is None:
        image = skio.imread(imname)

    # Edge apodization
    edge_apod = (np.sin(np.linspace(-np.pi/2, np.pi/2, w))+1)/2
#    Wx, Wy = np.meshgrid(*(np.concatenate([edge_apod, np.ones(l-2*w), np.flip(edge_apod)]) for l in image.shape)) #BEFORE
    Wx, Wy = np.meshgrid(*(np.concatenate([edge_apod, np.ones(l-2*w), np.flip(edge_apod)]) for l in np.flip(image.shape))) #AFTER
    W = Wx*Wy
    mean = image.mean()
    image = W*(image-mean) + mean


    # Do those two has any use ?
    image = image.astype('float32')
    image = image[:image.shape[0]-1+image.shape[0]%2,:image.shape[1]-1+image.shape[1]%2] #Odd array sizes

#    X, Y = np.meshgrid(*(np.linspace(-1, 1, l) for l in image.shape)) #J'ai un ordre différent de dans le code matlab #BEFORE
    X, Y = np.meshgrid(*(np.linspace(-1, 1, l) for l in np.flip(image.shape))) #AFTER
    R = np.sqrt(X**2 + Y**2)


    Ik = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(image))) #Pourquoi utiliser fftshift 2 fois???
    Ik = Ik*(R<1) # ??? That was not in the "manual" ...
    
    with np.errstate(divide='ignore', invalid='ignore'):
        Ikn =  Ik/np.abs(Ik)
    Ikn[np.isnan(Ikn)] = 0 #Nécessaire?
    Ikn[np.isinf(Ikn)] = 0 #Nécessaire?


    def decorr(Ik, Ikn, R, r, highpass_sigma=None):

        Nr = len(r)
        if highpass_sigma:
           Ik =  Ik*(1 - np.exp(-2*highpass_sigma**2*R**2))

        d=[]
        denom2 = np.sum(np.abs(Ik)**2)
        Ikn_conj = np.conj(Ikn)
        Ikn_abs_exp2 = np.abs(Ikn)**2
        
        for i in range(Nr):
    #         Mr = R<r[i]
    #         I1 = Ik[Mr]
    #         I2_conj = Ikn_conj[Mr]
    #         I2_abs_exp2 = Ikn_abs_exp2[Mr]
    #         d.append(
    #         np.real(I1*I2_conj).sum()/np.sqrt(np.sum(I2_abs_exp2)*denom2)
    #         )
            if i==0:
                Mr = R<r[i]
                I1 = Ik[Mr]
                I2_conj = Ikn_conj[Mr]
                I2_abs_exp2 = Ikn_abs_exp2[Mr]
                nom = np.real(I1*I2_conj).sum()
                denom1 = np.sum(I2_abs_exp2)
            else:
                Msk = np.logical_and(R<r[i], R>=r[i-1])
                I1 = Ik[Msk]
                I2_conj = Ikn_conj[Msk]
                I2_abs_exp2 = Ikn_abs_exp2[Msk]
                nom = nom + np.real(I1*I2_conj).sum()
                denom1 = denom1 + np.sum(I2_abs_exp2)
            with np.errstate(divide='ignore', invalid='ignore'):
                d.append(
                nom/np.sqrt(denom1*denom2)
                )
            
        
        
        d=np.array(d)
        d = np.floor(1000*d)/1000 # WHY ???
        d[np.isnan(d)] = 0
        return d


    r = np.linspace(0,1,Nr)
    d = decorr(Ik, Ikn, R, r)


    def decorr_peak(d, find_max_tol):
        idx = np.argmax(d)
        A = d[idx]
        while len(d) > 1:
            if ((A-np.min(d[idx:])) >= find_max_tol) or (idx==0):
                break
            else:
                d = d[:-1]
                idx = np.argmax(d)
                A = d[idx]
        return idx, A

    
    idx_0, A0 = decorr_peak(d, find_max_tol)
    r0 = r[idx_0]

    

    sigmas = np.exp(np.arange(Ng+1)/Ng*(np.log(2/r0)-np.log(0.15)) + np.log(0.15)) #Not like in the code, Ng+1 values???
    gMax = 2/r0

    if gMax==np.inf: gMax = max(image.shape[0],image.shape[1])/2
    sigmas = np.array([image.shape[0]/4] + list(np.exp(np.linspace(np.log(gMax),np.log(0.15),Ng))))
    

    idxs = np.array([])
#    As = np.array([A0]) #BEFORE
#    rs = np.array([r0]) #BEFORE
    As = np.array([A0]) #AFTER
    rs = np.array([r0]) #AFTER

    for i, sig in enumerate(sigmas):
        d = decorr(Ik, Ikn, R, r, highpass_sigma=sig)
        idx, A = decorr_peak(d, find_max_tol) #Note that in the matlab code, the first element in the array is SNR0, so the length of the array A0 is 12 instead of 11 here
        idxs = np.append(idxs, idx)
        As = np.append(As, A)
        rs = np.append(rs, r[idx])
    
#    print("line =", getframeinfo(currentframe()).lineno, "rs=", rs)
    max_freq_peak = rs.max()
    max_freq_peak_idx = np.where(rs == max_freq_peak)[0][-1]
#    if r0>max_freq_peak: #BEFORE
    if max_freq_peak_idx==0: #AFTER
        ind1 = 0
    elif max_freq_peak_idx>=(len(sigmas)-1): #AFTER: == changed to >=
        ind1 = max_freq_peak_idx-2 #AFTER: -1 changed to -2
    else:
        ind1 = max_freq_peak_idx-1
    ind2 = ind1+1
    
    
    r1 = rs[max_freq_peak_idx] - (r[1]-r[0])
    r2 = rs[max_freq_peak_idx] + 0.4
    r_finetune = np.linspace(r1, min(r2,r[-1]), Nr)

#    import pdb; pdb.set_trace()
#    print("line =", getframeinfo(currentframe()).lineno, "r_finetune=", r_finetune)
    sigmas_fine_tune = np.exp(np.linspace(np.log(sigmas[ind1]), np.log(sigmas[ind2]) ,Ng))
    
    As=As[:-1] #AFTER (weird...)
    rs=rs[:-1] #AFTER (weird...)
    for i, sig in enumerate(sigmas_fine_tune):
        d = decorr(Ik, Ikn, R, r_finetune, highpass_sigma=sig)
        idx, A = decorr_peak(d, find_max_tol) #Note that in the matlab code, the first element in the array is SNR0, so the length of the array A0 is 12 instead of 11 here
        idxs = np.append(idxs, idx)
        As = np.append(As, A)
        rs = np.append(rs, r_finetune[idx])

#    print("line =", getframeinfo(currentframe()).lineno, "rs=", rs)
    rs = np.append(rs, r0) #AFTER (nothing before)
    As = np.append(As, A0) #AFTER (nothing before)
    
    rs[As<0.05] = 0 # IN the matlab code, Ks set to zero too
    


    
    res = 2/np.max(rs)
    
    return res






    #TODO:
    #Compare my results with the matlab code
    #Do some final review/tests
