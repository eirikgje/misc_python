from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import fplcmod
import pyfits
import scipy.stats
import os
import subprocess
import shlex
from scipy import interpolate

def analytic(fname='taulikes_an.dat', lmax=30, path=None, lmax_comm=90):
    anpi = fplcmod.Pluginfo('ANALYTIC', lmin=2, lmax=lmax)
    pis = (anpi,)
    if path is None:
        fplcmod.update_files(pis, lmax=lmax_comm)
    else:
        fplcmod.update_files(pis, path=path, lmax=lmax_comm)
    fplcmod.run_program(app='tauspec')
    fplcmod.savefile(fname)

def brcg_new(fname='taulikes_saved.dat', lmax_like=30, firstchain=None, 
            lastchain=None,firstsamp=None, lastsamp=None, suffix='.fits', 
            path=None, lmax=90, comm_lmax=None, l_thresh=9, logfile=None):

    if suffix is None:
        sigmafile = None
        clfile = None
    else:
        sigmafile = 'sigma' + suffix
        clfile = 'cls' + suffix
    pi = fplcmod.Pluginfo('BRCG', lmin=2, lmax=lmax, lmax_like=lmax_like, 
            firstchain=firstchain, lastchain=lastchain, firstsamp=firstsamp, 
            lastsamp=lastsamp, lranges=(2, l_thresh), numranges=1, 
            rangetypes=('TT_TE_EE',), sigmafile=sigmafile, 
            clfile=clfile, spectra='TT_TE_EE', l_thresh=l_thresh,
            comm_lmax=comm_lmax)
    if path is None:
        fplcmod.update_files((pi,), lmax=lmax)
    else:
        fplcmod.update_files((pi,), path=path, lmax=lmax)
    fplcmod.run_program(app='tauspec', logfile=logfile)
    fplcmod.savefile(fname)

def brcg(fname='taulikes_onechain_v1.dat', lmax=30, firstchain=None, 
        lastchain=None, firstsamp=None, lastsamp=None, 
        suffix='_ctp3_real10_hke_onechain_v2.fits', path=None, lmax_br=9, 
        lmax_comm=90):

    if lmax > lmax_br:
        brpi = fplcmod.Pluginfo('BR_MOD', lmin=2, lmax=lmax_br, numranges=1,
                        lranges=(2, lmax_br), rangetypes=('TT_TE_EE',),
                        sigmafile = 'sigma' + suffix,
                        firstchain=firstchain, lastchain=lastchain,
                        firstsamp=firstsamp, lastsamp=lastsamp)
        cgpi = fplcmod.Pluginfo('COMM_GAUSS', lmin=lmax_br + 1, lmax=lmax, 
                                spectra='TT_TE_EE',
                                clfile='cls' + suffix,
                                firstchain=firstchain, lastchain=lastchain,
                                firstsamp=firstsamp, lastsamp=lastsamp)
        pis = (brpi, cgpi)
    else:
        brpi = fplcmod.Pluginfo('BR_MOD', lmin=2, lmax=lmax, numranges=1,
                        lranges=(2, lmax), rangetypes=('TT_TE_EE', ), 
                        sigmafile='sigma' + suffix, firstchain=firstchain,
                        lastchain=lastchain, firstsamp=firstsamp, 
                        lastsamp=lastsamp)
        pis = (brpi,)
    if path is None:
        fplcmod.update_files(pis, lmax=lmax_comm)
    else:
        fplcmod.update_files(pis, path=path, lmax=lmax_comm)
    fplcmod.run_program(app='tauspec')
    fplcmod.savefile(fname)


def run_allreal(suffix, fnamemidfix='allreal', numreal=10, burnin=0, lmax=30, 
                lmax_like=30, logprefix=None):
    sigmafile = 'datasets/brcg/sigma' + suffix
    clfile =  'datasets/brcg/cls' + suffix
    inithdu = pyfits.open(sigmafile)
    maxnumsamp = inithdu[0].header['NUMSAMP']
    numchains = inithdu[0].header['NUMCHAIN']
    numchains_per_real = numchains // numreal
    currfirstchain = 1
    for i in range(numreal):
        print 'Python: Real', i + 1
        if logprefix is None:
            logfile = None
        else:
            logfile = logprefix + '_%02d.log' % (i + 1)
        numsamp = int(inithdu[0].data[0, i, 0, 0])
        brcg_new(fname='taulikes_' + fnamemidfix + '_real%02d.dat'%(i + 1), 
                lmax=lmax, firstchain=currfirstchain, 
                lastchain=currfirstchain + numchains_per_real - 1, 
                firstsamp=1 + burnin, lastsamp=numsamp, suffix=suffix, 
                lmax_like=lmax_like, logfile=logfile)
        currfirstchain += numchains_per_real
    inithdu.close()

def run_allreal_convergence(suffix, numreal=10, burnin=0):
    sigmafile = 'sigma' + suffix
    clfile =  'cls' + suffix
    inithdu = pyfits.open(sigmafile)
    maxnumsamp = inithdu[0].header['NUMSAMP']
    numchains = inithdu[0].header['NUMCHAIN']
    inithdu.close()
    numchains_per_real = numchains // numreal
    currfirstchain = 1
    for i in range(numreal):
        brcg(fname='taulikes_allreal_pt1_real%02d.dat' % (i + 1), lmax=30, 
                firstchain=currfirstchain, 
                lastchain=currfirstchain + numchains_per_real - 1, 
                firstsamp=1 + burnin, lastsamp=(maxnumsamp-burnin)//2 + burnin,
                suffix=suffix)
        brcg(fname='taulikes_allreal_pt2_real%02d.dat' % (i + 1), lmax=30, 
                firstchain=currfirstchain, 
                lastchain=currfirstchain + numchains_per_real - 1, 
                firstsamp=(maxnumsamp-burnin)//2 + burnin + 1, 
                lastsamp=maxnumsamp, suffix=suffix)
        currfirstchain += numchains_per_real

def calculate_ml_and_asymm_errorbars_from_slices(x, lnL, area_fraction=0.68):
    dx = x[1] - x[0]
    lnL_norm = lnL / (np.sum(lnL) * dx)
    ml_ind = np.argmax(lnL_norm)
    ml = x[ml_ind]
    mean = np.sum(x * lnL_norm) * dx
    currpoint_high = ml_ind
    currpoint_low = ml_ind
    if ml_ind == 0:
        currpoint_high += 1
    elif ml_ind == len(x) - 1:
        currpoint_low -= 1
    elif lnL_norm[ml_ind + 1] > lnL_norm[ml_ind - 1]:
        currpoint_high += 1
    else:
        currpoint_low -= 1
    area = np.sum(lnL_norm[currpoint_low:currpoint_high + 1] * dx)
    while area < area_fraction:
        if currpoint_high == len(x) - 1 and currpoint_low == 0:
            raise ValueError("Cannot find bounds")
        if currpoint_high == len(x) - 1:
            currpoint_low -= 1
        elif currpoint_low == 0:
            currpoint_high += 1
        else:
            if lnL_norm[currpoint_high + 1] > lnL_norm[currpoint_low - 1]:
                currpoint_high += 1
            else:
                currpoint_low -= 1
        area = sum(lnL_norm[currpoint_low:currpoint_high + 1] * dx)

    upper = x[currpoint_high]
    lower = x[currpoint_low]

    return ml, upper, lower, mean

def calculate_ml_and_asymm_errorbars_from_samples(samples, x=None, area_fraction=0.68, smooth=True):
    #Assume samples is a fits file. Pretty simple now.
    if isinstance(samples, str):
        samples = pyfits.open(samples)[0].data.flatten()
    if smooth:
        kernel = scipy.stats.kde.gaussian_kde(samples)
        ml_ind = np.argmax(kernel(x))
        ml = x[ml_ind]
        dx = x[1] - x[0]
        currpoint_high = ml_ind
        currpoint_low = ml_ind
        y = kernel(x)
        if y[ml_ind + 1] > y[ml_ind - 1]:
            currpoint_high += 1
        else:
            currpoint_low -= 1
        area = sum(y[currpoint_low:currpoint_high + 1]*dx)
        while area < area_fraction:
            if (currpoint_high == len(x) - 1 and currpoint_low == 0):
                raise ValueError("Cannot find bounds")
            if currpoint_high == len(x) - 1:
                currpoint_low -= 1
            elif currpoint_low == 0:
                currpoint_high += 1
            else:
                if y[currpoint_high + 1] > y[currpoint_low - 1]:
                    currpoint_high += 1
                else:
                    currpoint_low -= 1
            area = sum(y[currpoint_low:currpoint_high + 1] * dx)
    
        upper = x[currpoint_high]
        lower = x[currpoint_low]
    else:
        sortedSamps = np.sort(samples)
        if x is None:
            x = np.linspace(sortedSamps[0], sortedSamps[-1], 10000)
        hist = scipy.stats.histogram(sortedSamps, numbins=len(x), defaultlimits=(x[0], x[-1]))
        ml = x[np.argmax(hist[0])]
        mlarg = np.searchsorted(sortedSamps, np.array(ml))
        numsamp_thresh = int(round(area_fraction * len(samples)))
        if mlarg == 0:
            lower = sortedSamps[0]
            upper = sortedSamps[numsamp_thresh]
            return ml, upper, lower
        if mlarg < numsamp_thresh:
            currminbracket = sortedSamps[numsamp_thresh] - sortedSamps[0]
            currminind = 0
            start = 1
            if mlarg + numsamp_thresh >= len(samples):
                stop = len(samples) - numsamp_thresh
            else:
                stop = mlarg + 1
        elif mlarg + numsamp_thresh < len(samples):
            currminbracket = sortedSamps[mlarg] - sortedSamps[mlarg - numsamp_thresh]
            currminind = mlarg - numsamp_thresh
            start = mlarg - numsamp_thresh + 1
            stop = mlarg + 1
        elif mlarg + numsamp_thresh >= len(samples):
            currminbracket = sortedSamps[mlarg] - sortedSamps[mlarg - numsamp_thresh]
            currminind = mlarg - numsamp_thresh
            start = mlarg - numsamp_thresh + 1
            stop = len(samples) - numsamp_thresh
        for i in range(start, stop):
            currbrack = sortedSamps[numsamp_thresh + i] - sortedSamps[i]
            if currbrack < currminbracket:
                currminbracket = currbrack
                currminind = i
        upper = sortedSamps[currminind + numsamp_thresh]
        lower = sortedSamps[currminind]

#        if mlarg - numsamp_thresh / 2 < 0:
#            lower = sortedSamps[0]
#            upper = sortedSamps[numsamp_thresh]
#        elif mlarg + numsamp_thresh / 2 > len(sortedSamps):
#            upper = sortedSamps[-1]
#            lower = sortedSamps[len(sortedSamps) - numsamp_thresh]
#        else:
#            upper = sortedSamps[mlarg + numsamp_thresh / 2]
#            lower = sortedSamps[mlarg - numsamp_thresh / 2]
#        currhigh = mlarg
#        if mlarg == 0:
#            currlow = 0
#        else:
#            currlow = mlarg - 1
#            num = 2
#        while num < numsamp_thresh:
#            if currlow == 0:
#                currhigh += 1
#            elif currhigh == len(samples):
#                currlow -= 1
#            else:
#                if (sortedSamps[currhigh + 1] - ml) > (ml - sortedSamps[currlow - 1]):
#                    currlow -= 1
#                else:
#                    currhigh += 1
#            num += 1
#        upper = sortedSamps[currhigh]
#        lower = sortedSamps[currlow]

    return (ml, upper, lower)
    
def collect_statistics(samples, x, area_fraction=0.68, smooth=True):
    if isinstance(samples, str):
        samples = pyfits.open(samples)[0].data.flatten()

    mean = np.mean(samples)
    std = np.std(samples)
    (ml, upper, lower) = calculate_ml_and_asymm_errorbars_from_samples(samples, x, area_fraction, smooth=smooth)
    return mean, std, ml, upper, lower

def collect_statistics_from_slice(x, y, area_fraction=0.68):
    dx = x[1] - x[0]
    ynorm = y / (np.sum(y) *dx)
    mean = np.sum(ynorm * x) *dx
    std = np.sqrt(np.sum((x-mean) ** 2 * ynorm) * dx)
    ml, upper, lower, dum = calculate_ml_and_asymm_errorbars_from_slices(x, y, area_fraction)
    return mean, std, ml, upper, lower

def replace_line_in_file(fname, startswith, replacement):
    """Will overwrite tempfile.txt."""
    subprocess.call(shlex.split("cp %s tempfile.txt" % fname))
    infile = open('tempfile.txt')
    outfile = open(fname, 'w')
    i = 0
    for line in infile:
        i += 1
        written = False
        for entry, repl in zip(startswith, replacement):
            if isinstance(entry, int):
                if i == entry:
                    outfile.write(repl)
                    written = True
                    break
            elif isinstance(entry, str):
                if line.startswith(entry):
                    outfile.write(repl)
                    written = True
                    break
        if not written:
            outfile.write(line)
    subprocess.call(shlex.split("rm tempfile.txt"))
    infile.close()
    outfile.close()

def cp(fname1, fname2):
    subprocess.call(shlex.split('cp ' + fname1 + ' ' + fname2))

def rm(fname):
    subprocess.call(shlex.split('rm ' + fname))

def flatten_commander_chain(chain, burnin=0):
    if isinstance(chain, str):
        chain = pyfits.open(chain)[0].data
    numsamps = 0
    for i in range(chain.shape[1]):
        numsamps += chain[0, i, 0, 0] - burnin
    shape = (numsamps, chain.shape[2], chain.shape[3])
    flatchain = np.zeros(shape)
    k = 0
    for i in range(chain.shape[1]):
        for j in range(burnin, chain[0, i, 0, 0]):
            flatchain[k] = chain[j + 1, i]
            k += 1
    return flatchain

def get_noise_powerspec(rms, lmax=47, nside=16, beam_fwhm=None, pol=False):
    "FWHM in arcmin"
    l = np.arange(0, lmax + 1)
    npix = 12 * nside ** 2
    spec = 4 * np.pi / npix * rms ** 2 * l * (l + 1) / (2 * np.pi)
    if beam_fwhm is not None:
        sigma = beam_fwhm/(60.0 * 180.0) * np.pi / np.sqrt(8.0*np.log(2.0))
        sig2 = sigma ** 2
        g = np.exp(-0.5*l*(l+1.0)*sig2)
        if pol:
            factor_pol = np.exp([0.0, 2.0*sig2, 2.0*sig2])
            gout = g * np.exp(2.0 *sig2)
        else:
            gout = g
        spec = spec / g**2
    return spec

def cambcls_to_comminit(fnamein, fnameout, mode='TT_EE_TE'):
    if mode == 'TT_EE_TE':
        a = np.loadtxt(fnamein)
        b = np.zeros((len(a), 13))
        b[:, 0] = a[:, 0]
        b[:, 1] = a[:, 1]
        b[:, 2] = a[:, 1] / 10.0
        b[:, 3] = a[:, 3]
        b[:, 4] = a[:, 3] / 10.0
        b[:, 6] = a[:, 3] / 100.0
        b[:, 7] = a[:, 2]
        b[:, 8] = a[:, 2] / 10.0
        b[:, 12] = a[:, 2] / 20.0
    np.savetxt(fnameout, b)

def normalize_2d_probdist(dist):
    dx1 = dist[1, 0, 0] - dist[0, 0, 0]
    dx2 = dist[0, 1, 1] - dist[0, 0, 1]
    norm = sum(sum(dist))[2] * dx1 * dx2
    dist[:, :, 2] = dist[:, :, 2] / norm
    return dist

def calc_normprobdiff(dist1, dist2):
    numbins = int(np.sqrt(dist1.shape[0]))
    ndist1 = np.reshape(dist1, (numbins, numbins, 3))
    ndist2 = np.reshape(dist2, (numbins, numbins, 3))
    ndist1 = normalize_2d_probdist(ndist1)
    ndist2 = normalize_2d_probdist(ndist2)
    dx1 = ndist1[1, 0, 0] - ndist1[0, 0, 0]
    dx2 = ndist1[0, 1, 1] - ndist1[0, 0, 1]
    absdiff = sum(sum(abs(ndist1[:, :, 2] - ndist2[:, :, 2]))) * dx1 * dx2
    return absdiff

def bestline(x, y):
    #Returns the best fit linear coefficients that describe a line 
    #through the points x_i, y_i, assuming error only in the y direction

    if len(x) != len(y):
        raise ValueError('x and y must have same length')
    xi = np.zeros((2, 2))
    xy = np.zeros(2)
    xi[0, 0] = len(x)
    xi[0, 1] = np.sum(x)
    xi[1, 0] = np.sum(x)
    xi[1, 1] = np.sum(x**2)
    xy[0] = np.sum(y)
    xy[1] = np.sum(x * y)
    xi = np.matrix(xi)
    xy = np.matrix(xy)
    return np.array(xi.I * xy.T)

def marginalize_2d_dist(dist, skiprows=1):
    if isinstance(dist, str):
        dist = np.loadtxt(dist, skiprows=skiprows)
    nbins = int(np.sqrt(np.shape(dist)[0]))
    dist = np.reshape(dist, (nbins, nbins, 3))
#    dist = normalize_2d_probdist(dist)
    par = np.zeros((4, nbins))
    for i in range(nbins):
        par[0, i] = dist[i, 0, 0]
        par[1, i] = np.sum(dist[i, :, 2])
        par[2, i] = dist[0, i, 1]
        par[3, i] = np.sum(dist[:, i, 2])

    par[1, :] = par[1, :] / np.sum(par[1, :] * (par[0, 1] - par[0, 0]))
    par[3, :] = par[3, :] / np.sum(par[3, :] * (par[2, 1] - par[2, 0]))

    return par

def calc_alm_chisq(alm, cl):
    #Calculates the l-by-l chisquared of the alms relative to the cls
    #Assumes alms are lmax, lmax, 2 and cls are lmax

    chisq = np.zeros(len(cl))
    for i in range(len(cl)):
        l = i + 2
        for m in range(l):
            if m == 0:
                chisq[i] += alm[i, m, 0] ** 2 / cl[i]
            else:
                for j in range(2):
                    chisq[i] += 2 * alm[i, m, j] ** 2 / cl[i]
        chisq[i] = chisq[i] / (2*l + 1) * (l * (l + 1)) / (2 * np.pi)
    return chisq

def calc_alm_chisq_fromfits(almfile, clfile):
    alm = pyfits.open(almfile)[0].data
    cls = pyfits.open(clfile)[0].data
    numiter = alm.shape[0]
    numchain = alm.shape[1]
    alm = np.reshape(alm, (numiter*numchain,alm.shape[2], alm.shape[3], alm.shape[4], alm.shape[5]))
    alm = alm[:, :, :, 2:, :]
    cls = cls[1:]
    cls = np.reshape(cls, (numiter * numchain, cls.shape[2], cls.shape[3]))
    if alm.shape[1] == 3:
        cls = np.concatenate((cls[:, 0:1, :], cls[:, 3:4, :], cls[:, 5:6, :]), 1)
    elif alm.shape[1] == 1:
        cls = cls[:, 0:1, :]
    cls = cls[:, :, 2:]
    cls = np.transpose(cls).copy()
    alm = np.transpose(alm).copy()
    chisq = np.zeros(cls.shape)
    for i in range(cls.shape[0]):
        l = i + 2
        for m in range(l):
            if m == 0:
                chisq[i, :, :] += alm[0, i, m, :, :] ** 2
            else:
                chisq[i, :, :] += np.sum(2 * alm[:, i, m, :, :] ** 2, 0)
        chisq[i, :, :] = chisq[i, :, :] / cls[i, :, :] / (2 * l + 1) * (l * (l + 1)) / (2 * np.pi)
    return chisq

def reduce_alm_chisq(chisq):
    totdf = 0
    redchi = np.zeros((chisq.shape[1], chisq.shape[2]))
    for i in range(len(chisq)):
        l = i + 2
        redchi[:, :] += chisq[i, :, :] * (2 *l + 1)
        totdf += 2 * l + 1
    redchi /= totdf
    print totdf
    return redchi

def separate_contours_data(fname, islog=False):
    file = open(fname, 'r')

    numpoint = file.readline()

    numpoint = np.fromstring(numpoint, sep=' ', dtype='int')

    file.close()
    
    Q_args = np.zeros(numpoint[0])
    N_args = np.zeros(numpoint[1])
    lnL = np.zeros((numpoint[1], numpoint[0]))
    
    data = np.loadtxt(fname, skiprows=1)
    
    for i in range(numpoint[0]):
        for j in range(numpoint[1]):
            Q_args[i] = data[i*numpoint[1] + j][0]
            N_args[j] = data[i*numpoint[1] + j][1]
            if islog:
                lnL[j, i] = data[i*numpoint[1] + j][2]
            else:
                lnL[j, i] = np.log(data[i*numpoint[1] + j][2])
    return lnL, Q_args, N_args

def calc_1d_mean_and_symm_sigmas(x, like):
    dx = x[1] - x[0]
    mean = np.sum(like * x) * dx
    sigma = np.sqrt(np.sum((x - mean) ** 2 * like) * dx)

    return mean, sigma

def get_splined_2d_dist(fname, islog=False, num_spline_points=100):
    data = np.loadtxt(fname, skiprows=1)
    if islog:
        lnL = data[:, 2]
    else:
        lnL = np.log(data[:, 2])
    lnL = -2* (lnL - np.max(lnL))
    tck = interpolate.bisplrep(data[:, 0], data[:, 1], lnL, s=1)
    num_spline_points = complex(0, num_spline_points)
    Q_args, N_args = np.mgrid[data[0, 0]:data[-1, 0]:num_spline_points, data[0, 1]:data[-1, 1]:num_spline_points]
    lnL_splined = interpolate.bisplev(Q_args[:, 0], N_args[0, :], tck)

    return Q_args, N_args, lnL_splined

def spline_and_marginalize_2d_dist(fname, skiprows=1):
#    if isinstance(dist, str):
#        dist = np.loadtxt(dist, skiprows=skiprows)
    Q_args, N_args, lnL_splined = get_splined_2d_dist(fname)
#    nbins = int(np.sqrt(np.shape(dist)[0]))
#    dist = np.reshape(dist, (nbins, nbins, 3))
#    dist = normalize_2d_probdist(dist)
    like = np.exp(-0.5*lnL_splined)
    par = np.zeros((4, len(Q_args)))
    for i in range(len(Q_args)):
        par[0, i] = Q_args[i, 0]
        par[1, i] = np.sum(like[i, :])
        par[2, i] = N_args[0, i]
        par[3, i] = np.sum(like[:, i])

    par[1, :] = par[1, :] / np.sum(par[1, :] * (par[0, 1] - par[0, 0]))
    par[3, :] = par[3, :] / np.sum(par[3, :] * (par[2, 1] - par[2, 0]))

    return par

def collect_statistics_from_sigma_samples(sigma, burnin=0, smooth=True, area_fraction=0.68, numpoints=1000):
#    sigma = sigma[1:, :, :, :]
    sigma = flatten_commander_chain(sigma, burnin)
    lmax = sigma.shape[2] - 1
    means = []
    stds = []
    mls = []
    uppers = []
    lowers = []
    for l in range(2, lmax + 1):
        print l
        samps = sigma[:, 0, l].flatten()
        x = np.linspace(np.min(samps), np.max(samps), numpoints)
        print np.min(samps), np.max(samps)
        mean, std, ml, upper, lower = collect_statistics(samps, x, area_fraction=area_fraction, smooth=smooth)
        means.append(mean)
        stds.append(std)
        mls.append(ml)
        uppers.append(upper)
        lowers.append(lower)
    return np.array(means), np.array(stds), np.array(mls), np.array(uppers), np.array(lowers)

def collect_statistics_from_sigma_bins(sigma, bins_start, bins_end, burnin=0, smooth=True, area_fraction=0.68, numpoints=1000):
    sigma = flatten_commander_chain(sigma, burnin)
    lmax = sigma.shape[2] - 1
    means = []
    stds = []
    mls = []
    uppers = []
    lowers = []
    for lstart, lend in zip(bins_start, bins_end):
        vars = []
        sigmas = []
        for l in range(lstart, lend+1):
            print l
            vars.append(np.var(sigma[:, 0, l]))
            print vars[-1]
            sigmas.append(sigma[:, 0, l] / vars[-1])
        vars = np.array(vars)
        sigmas = np.array(sigmas)
        if lstart == lend:
            samps = sigmas
            samps = samps * vars
        else:
            samps = np.sum(sigmas, axis=0)
            samps = samps / np.sum(1 / vars)
        print np.min(samps)
        print np.max(samps)
        x = np.linspace(np.min(samps), np.max(samps), numpoints)
        mean, std, ml, upper, lower = collect_statistics(samps, x, area_fraction=area_fraction, smooth=smooth)
        means.append(mean)
        stds.append(std)
        mls.append(ml)
        uppers.append(upper)
        lowers.append(lower)
    return means, stds, mls, uppers, lowers

def corr(chain):
    c = []
    mean = np.mean(chain)
    var = np.var(chain)
    for i in range(len(chain) - 1):
        c.append(np.mean((chain[:len(chain)-(i + 1)] - mean) * \
                (chain[i+1:] - mean))/var)
    return np.array(c)

def calc_multipole_correlation_length_from_fits_files(file, lmax=None, spec=0, burnin=0, corrlength_thresh=0.1):
    data = flatten_commander_chain(file, burnin)
    if lmax is None:
        lmax = data.shape[2] - 1 
    corrlengths = np.empty(0)
    for l in xrange(2, lmax + 1):
        currcorr = np.where(corr(data[:, spec, l]) < corrlength_thresh)
        if np.size(currcorr) == 0:
            corrlengths = np.append(corrlengths, -1)
        else:
            corrlengths = np.append(corrlengths, currcorr[0])
    return np.arange(2, lmax+1), corrlengths

def calc_corr_raw(l=2, ps='cls', lmax=90, spec=1, chain=None, burnin=None):
    if burnin is None:
        burnin = 0
    dat = []
    if ps == 'cls':
        ind = l - 2
    elif ps == 'sigma':
        ind = l - 1 + lmax - 2
    else:
        raise ValueError("ps has unknown value")
    if chain is None:
        chain = 1
    filenum = 1
    cont = True
    filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
    cont = os.path.exists(filename)
    while cont:
        temp = np.loadtxt(filename)
        dat.append(temp[ind, spec])
        filenum += 1
        filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
        cont = os.path.exists(filename)
    return corr(dat)
