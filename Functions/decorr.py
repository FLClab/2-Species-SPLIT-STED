
import numpy

def linmap(val, mapMin, mapMax):

    valMin = numpy.min(val)
    valMax = numpy.max(val)

    tempVal = (val - valMin) / (valMax - valMin)

    map0 = tempVal < 0
    map1 = tempVal > 1

    tempVal[map0] = 0
    tempVal[map1] = 1

    return tempVal * (mapMax - mapMin) + mapMin

def apodImRect(im, N):

    Ny, Nx = im.shape

    x = numpy.abs(numpy.linspace( -Nx / 2, Nx / 2, Nx))
    y = numpy.abs(numpy.linspace( -Ny / 2, Ny / 2, Ny))

    mapx = x > Nx / 2 - N
    mapy = y > Ny / 2 - N

    val = numpy.mean(im)

    d = (-1 * numpy.abs(x) - numpy.mean(-1 * numpy.abs(x[mapx]))) * mapx
    d = linmap(d, -numpy.pi/2, numpy.pi/2)
    d[~mapx] = numpy.pi / 2
    maskx = (numpy.sin(d) + 1) / 2

    d = (-1 * numpy.abs(y) - numpy.mean(-1 * numpy.abs(y[mapy]))) * mapy
    d = linmap(d, -numpy.pi/2, numpy.pi/2)
    d[~mapy] = numpy.pi / 2
    masky = (numpy.sin(d) + 1) / 2
    masky = masky[:, numpy.newaxis]

    mask = numpy.tile(masky, (1, Nx)) * numpy.tile(maskx, (Ny, 1))

    return (im - val) * mask + val

def getCorrcoef(I1, I2, c1):

    c2 = numpy.sqrt(numpy.sum(numpy.abs(I2) ** 2))
    if c1 * c2 == 0:
        return 0

    cc = numpy.sum(numpy.real(I1 * numpy.conj(I2))) / (c1 * c2)
    return numpy.floor(cc * 1000) / 1000

def getDCorrLocalMax(d):

    Nr = len(d)
    if Nr < 2:
        ind = 0
        A = d[ind]
    else:
        t = d.copy()
        ind = numpy.argmax(t)
        A = t[ind]
        while len(t) > 1:
            if ind == len(t):
                t = t[:-1]
                ind = numpy.argmax(t)
                A = t[ind]
            elif ind == 0:
                break
            elif (A - min(t[ind : ])) >= 0.0005:
                break
            else:
                t = t[:-1]
                ind = numpy.argmax(t)
                A = t[ind]
    return ind, A

def getDCorr(im, r, Ng=10):

    im = im.astype(numpy.float32)

    im = im[:im.shape[0]-1+im.shape[0]%2,:im.shape[1]-1+im.shape[1]%2]
    X, Y = numpy.meshgrid(*(numpy.linspace(-1, 1, l) for l in im.shape[::-1]))
    R = numpy.sqrt(X ** 2 + Y ** 2)

    Nr = len(r)

    fourier = numpy.fft.fftshift(numpy.fft.fftn(numpy.fft.fftshift(im)))
    In = fourier / numpy.abs(fourier)
    In = numpy.nan_to_num(In, nan=0., posinf=0., neginf=0.)

    mask0 = R ** 2 < 1
    In = mask0 * In

    Ik = mask0 * fourier

    c = numpy.sqrt(numpy.sum(Ik * numpy.conj(Ik)))
    c = numpy.real(c) # This ensures a real number

    r0 = numpy.linspace(r[0], r[-1], Nr)

    d0 = []
    for value in r0:
        cc = getCorrcoef(Ik, (R ** 2 < value ** 2) * In, c)
        if numpy.isnan(cc):
            cc = 0
        d0.append(cc)
    d0 = numpy.array(d0)
    ind0, snr0 = getDCorrLocalMax(d0)
    k0 = r[ind0]

    if r0[ind0] == 0.:
        gMax = max(im.shape)/2
    else:
        gMax = 2 / r0[ind0]

    # Search highest frequency peak
    g = [im.shape[0] / 4, *numpy.exp(numpy.linspace(numpy.log(gMax), numpy.log(0.15), Ng))]
    d = numpy.zeros((Nr, 2 * Ng + 1))

    kc = numpy.zeros(len(g) + Ng + 1)
    kc[0] = k0
    SNR = numpy.zeros(len(g) + Ng + 1)
    SNR[0] = snr0

    ind0 = 1
    for refin in range(2):

        for h in range(len(g)):
            Ir = Ik * (1 - numpy.exp(-2 * g[h] * g[h] * R ** 2))
            c = numpy.sqrt(numpy.sum(numpy.abs(Ir) ** 2))

            for k in reversed(range(ind0 - 1, len(r))):

                mask = (R ** 2 < r[k] ** 2)
                cc = getCorrcoef(Ir[mask], In[mask], c)
                if numpy.isnan(cc):
                    cc = 0
                d[k, h + Ng * refin] = cc

            ind, snr = getDCorrLocalMax(d[k :, h + Ng * refin])
            ind = ind + k

            kc[h + Ng * refin + 1] = r[ind]
            SNR[h + Ng * refin + 1] = snr

        if refin == 0:
            # In original matlab code, kc is built on the fly hence min(indmax, len(g)) is not required
            indmax = numpy.argwhere(kc == max(kc)).ravel()
            ind1 = min(indmax[-1], len(g))

            if ind1 == 0: # Peaks only without highpass
                ind1 = 0
                ind2 = 1
                g1 = im.shape[0]
                g2 = g[0]
            elif ind1 >= len(g):
                ind2 = ind1 - 1
                ind1 = ind1 - 2
                g1 = g[ind1]
                g2 = g[ind2]
            else:
                ind2 = ind1
                ind1 = ind1 - 1
                g1 = g[ind1]
                g2 = g[ind2]
            g = numpy.exp(numpy.linspace(numpy.log(g1), numpy.log(g2), Ng))

            r1 = kc[indmax[-1]] - (r[1] - r[0])
            r2 = kc[indmax[-1]] + 0.4

            if r1 < 0:
                r1 = 0
            if r2 > 1:
                r2 = 1
            r = numpy.linspace(r1, min(r2, r[-1]), Nr)
            ind0 = 1
            r2 = r

    kc = numpy.append(kc, k0)
    SNR = numpy.append(SNR, snr0)

    kc[SNR < 0.05] = 0
    SNR[SNR < 0.05] = 0

    snr = SNR

    ind = numpy.argmax(kc)
    kcMax = kc[ind]
    AMax = SNR[ind]
    A0 = snr0

    # print(2 / kcMax * 20, A0, AMax)

    return kcMax, A0

def calculate(image, N=20, Nr=50, Ng=10):
    pps=5
    r = numpy.linspace(0, 1, Nr)

    image = apodImRect(image, N)
    KcMax, A0 = getDCorr(image, r, Ng)

    return 2 / KcMax

if __name__ == "__main__":

    import tifffile
    im = tifffile.imread("test_image.tif")
    calculate(im)

