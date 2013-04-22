#Module for various plot routines used all over the place
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import gelman_rubin
import corrmod
import pyfits
import scipy.integrate as integrate
import gen_utils
from scipy import interpolate

def plot_raw_chains(mode='multipole', l=None, ps='cls', lmax=90, spec=1,
                    chains=None, burnin=0, lmin=2, numsamples=None, clf=True):
    if ps == 'cls':
        if mode == 'multipole':
            ind = l - 2
        elif mode == 'spectrum':
            ind = np.arange(lmin-2, lmax-1)
    elif ps == 'sigma':
        if mode == 'multipole':
            ind = l - 1 + lmax - 2
        elif mode == 'spectrum':
            ind = np.arange(lmax - 1 + lmin-2, 2 * lmax - 2)
    else:
        raise ValueError("ps has unknown value")
    if chains is None:
        chains = (1,)
    dat = []
    for chain in chains:
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
            if numsamples is not None:
                cont = cont and filenum <= numsamples + burnin
    if clf:
        plt.clf()
    if mode == 'multipole':
        plt.plot(dat)
    elif mode == 'spectrum':
        l = np.arange(lmin, lmax + 1)
        for i in range(len(dat)):
            plt.plot(l, dat[i])


def plot_raw_all_real(l=None, ps='cls', lmax=90, spec=1,
                        chain_start=1, chain_end=100, chain_intervals=10, 
                        burnin=None, separate_plots=False):
    if burnin is None:
        burnin = 0
    if ps == 'cls':
        ind = l - 2
    elif ps == 'sigma':
        ind = l -2 + lmax - 1
    else:
        raise ValueError("ps has unknown value")
    minnum = 1000000
    
    tdat = []
    for num in range(chain_start, chain_end, chain_intervals):
        dat = []
        chains = range(num, num + chain_intervals)
        filenum = 0
        for chain in chains:
            currfilenum = 1
            cont = True
            filename = 'cl_c%04d_k%05d.dat' % (chain, currfilenum + burnin)
            cont = os.path.exists(filename)
            while cont:
                temp = np.loadtxt(filename)
                dat.append(temp[ind, spec])
                currfilenum += 1
                filename = 'cl_c%04d_k%05d.dat' % (chain, currfilenum + burnin)
                cont = os.path.exists(filename)
            filenum += currfilenum - 1
        tdat.append(dat)
        if filenum < minnum:
            minnum = filenum
    ndat = np.zeros((len(tdat), minnum))
    for i in range(len(tdat)):
        for j in range(minnum):
            ndat[i, j] = tdat[i][j]
    for i in range(ndat.shape[0]):
        if separate_plots:
            plt.figure()
        plt.plot(ndat[i])
    
def plot_gr(lmin=2, lmax=90, ps='cls', spec=1, burnin=None, chain_start=1, chain_end=100, chain_intervals=10, separate_plots=False):
    totgelm = gelman_rubin.calc_gr_all_multipoles_all_chains_raw(lmax=lmax, lmin=lmin, ps=ps, spec=spec, burnin=burnin, chain_start=chain_start, chain_end=chain_end, chain_intervals=chain_intervals)
    for i in range(len(totgelm)):
        if separate_plots:
            plt.figure()
        plt.plot(totgelm[i])

def plot_raw_correlations(l, ps='cls', chains=None, lmax=90, spec=1, burnin=0, lmin=2):
    if ps == 'cls':
        ind = np.array(l) - 2
    elif ps == 'sigma':
        ind = np.array(l) - 1 + lmax - 2 
    else:
        raise ValueError("ps has unknown value")
    dat = []
    dat2= []
    for chain in chains:
        filenum = 1
        cont = True
        filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
        cont = os.path.exists(filename)
        while cont:
            temp = np.loadtxt(filename)
            dat.append(temp[ind[0], spec])
            dat2.append(temp[ind[1], spec])
            filenum += 1
            filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
            cont = os.path.exists(filename)
    plt.scatter(dat, dat2)

def plot_raw_corrcoff(l, ps='cls', chain=None, lmax=90, spec=1, burnin=0):
    dat = corrmod.calc_corr_raw(l=l, ps=ps, chain=chain, lmax=lmax, 
                                spec=spec, burnin=burnin)
    plt.clf()
    plt.plot(dat[:len(dat)//2])

def plot_raw_histogram(l, ps='cls', chains=None, lmax=90, spec=1, burnin=0, bins=30, normed=False, xrange=None):
    if ps == 'cls':
        ind = l - 2
    elif ps == 'sigma':
        ind = l - 1 + lmax - 2 
    else:
        raise ValueError("ps has unknown value")
    if chains is None:
        chains = (1,)
    dat = []
    for chain in chains:
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
    plt.clf()
    plt.hist(dat, bins=bins, range=xrange, normed=normed)

def plot_fits_chains(fname, mode='multipole', l=None, lmax=90, spec=1, 
                    chain_range=[0, 1], start_sample=1, num_samples=None, 
                    lmin=2, clf=True, xlim=None, ylim=None, thinning=None, 
                    linewidth=1):
    if mode == 'multipole':
        ind = l 
    elif mode == 'spectrum':
        ind = np.arange(lmin, lmax + 1)
    else:
        raise ValueError("mode has unknown value")
    hdulist = pyfits.open(fname)
    data = hdulist[0].data
    if thinning is None:
        thinning = 1
    if data.shape[1] == 1:
        if num_samples is None:
            last_sample = data.shape[0]
        else:
            last_sample = start_sample + num_samples * thinning
        if mode == 'multipole':
            dat = np.reshape(data[start_sample:last_sample:thinning, 0, spec - 1, ind], ((num_samples)))
        elif mode == 'spectrum':
            dat = np.reshape(data[start_sample:last_sample:thinning, 0, spec - 1, ind], (num_samples, len(ind)))
    elif data.shape[1] > 1:
        #Means there are several chains and we must read number of samples form each chain
        if num_samples is None:
            last_sample = 1000000
        else:
            last_sample = start_sample + num_samples * thinning
        dat = []
        for chain in range(chain_range[0], chain_range[1]):
            last_samp = min(last_sample, int(data[0, chain, 0, 0]))
            for samp in range(start_sample, last_samp, thinning):
                dat.append(data[samp, chain, spec-1, ind])
    hdulist.close()
    if clf:
        plt.clf()

    if mode == 'multipole':
        plt.plot(dat)
    elif mode == 'spectrum':
        l = np.arange(lmin, lmax + 1)
        for i in range(len(dat)):
            plt.plot(l, dat[i], linewidth=linewidth)
    if xlim is not None:
        plt.xlim(xlim)

def plot_fits_histogram(fname, l=None, spec=1, bins=30, chain_range=[0, 1], 
                        start_sample=1, last_sample=None, clf=True, xlim=None, 
                        histtype='bar', normed=False, color=None, ylim=None,
                        xrange=None, label=None):
    if not isinstance(l, int):
        raise ValueError("l must be an integer")
    ind = l 
    hdulist = pyfits.open(fname)
    data = hdulist[0].data
    if data.shape[1] == 1:
        if last_sample is None:
            last_sample = data.shape[0] - 1
        dat = np.reshape(data[start_sample:last_sample + 1, 0, spec - 1, ind], ((last_sample - start_sample + 1)))
    elif data.shape[1] > 1:
        #Means there are several chains and we must read number of samples form each chain
        dat = []
        ls = last_sample is None
        for chain in range(chain_range[0], chain_range[1]):
            if ls:
                last_sample = int(data[0, chain, 0, 0])
            for samp in range(start_sample, last_sample + 1):
                dat.append(data[samp, chain, spec-1, ind])

    hdulist.close()
    if clf:
        plt.clf()

    dat = plt.hist(dat, bins=bins, histtype=histtype, normed=normed, color=color, range=xrange, label=label)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    return dat

def plot_lorisspectra(numrange=None, spec=1, lmax=30):
    prefix = '/mn/svati/d1/eirikgje/data/spline_project/slicer_data/ctp3_trieste/loris_powerspec/'
    if numrange is None:
        numrange = [0, 100]
    taudiff = 150
    taustart = 4150
    for i in range(numrange[0], numrange[1]):
        file = prefix + 'cls_ctp3_fiducial_tau.%05d_scalCls.dat' % (taustart + i * taudiff)
        curr = np.loadtxt(file)
        plt.plot(curr[:lmax-1, 0], curr[:lmax-1, spec], linewidth=2)

def integrand(x, C, D, l):
    return np.exp((2*l-3)/2*np.log(x) - 1/x - (2*l-1)/2*np.log((x-C)**2 + D**2))
#    return x**((2*l-3)/2)*np.exp(-1/x)/(((x-C)**2 + D**2)**((2*l-1)/2))
def I2(l, C, D):
    result = np.zeros(len(C))
    for i in range(len(C)):
        result[i], err = integrate.quad(integrand, 0, integrate.Inf, (C[i], D[i], l))
    return result

def plot_cl_marginals(l, sigma, spec, cl):
    #Assumes a specific shape of sigma
    if spec == 'TT':
        sig = sigma[l-2, 1]
    elif spec == 'EE':
        sig = sigma[l-2, 2]
    elif spec == 'TE':
        sig = [sigma[l-2, 1], sigma[l-2, 2], sigma[l-2, 3]]
    if spec in ('TT', 'EE'):
        dist = sig**((2*l -3)/2)/cl**((2*l-1)/2)*np.exp(-(2*l + 1)/2*sig/cl)
    elif spec == 'TE':
        det = sig[0]*sig[1] - sig[2]**2
        C = sig[2]*cl/(sig[0]*sig[1])
        D = np.sqrt(det)*cl/((2*l+1)/2*sig[0]*sig[1])
        dist = det**((2*l-2)/2)/((sig[0]*sig[1])**((2*l-1)/2)) * I2(l, C, D)
    #normalize
    norm = sum(dist*(cl[1]-cl[0]))
    dist = dist / norm
    plt.plot(cl, dist, linewidth=3)
    return dist

def plot_cl_single_marginal(l, sig, spec, cl):
    ##TT, EE only
    dist = sig**((2*l -3)/2)/cl**((2*l-1)/2)*np.exp(-(2*l + 1)/2*sig/cl)
    #normalize
    norm = sum(dist*(cl[1]-cl[0]))
    dist = dist / norm
    plt.plot(cl, dist, linewidth=3)
    return dist


def plot_hists_and_marginals(histfname, sigma, l, chain_range=[0, 1], start_sample=1, last_sample=None, clf=True, xrange=None):
    for spec, num in zip(('TT', 'TE', 'EE'), (1, 2, 4)):
        plt.figure()
        dat = plot_fits_histogram(histfname, l, spec=num, chain_range=chain_range, start_sample=start_sample, last_sample=last_sample, xrange=xrange, normed=True)
        cl = np.linspace(dat[1][0], dat[1][-1], 500)
        plot_cl_marginals(l, sigma, spec, cl)

def plot_diag_noise_powspec(rms, lmax=47, nside=16, beam_fwhm=None, pol=False):
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

    plt.plot(l, spec, linewidth=3)

def plot_lorisslices(num, label='Loris', nside=32, color=None, linewidth=1, normalize=False):
    if nside == 32:
        lorisdata = np.loadtxt(
        '/mn/svati/d1/eirikgje/data/spline_project/slicer_data/ctp3_trieste/loris_slices/tau_slice_%04d.dat' % num)
        if normalize:
            dx = lorisdata[1, 0] - lorisdata[0, 0]
            area = sum(lorisdata[:, 1] * dx)
            lorisdata[:, 1] = lorisdata[:, 1] / area
        if color is None:
            p = plt.plot(lorisdata[:, 0], lorisdata[:, 1], label=label, linewidth=linewidth)
        else:
            p = plt.plot(lorisdata[:, 0], lorisdata[:, 1], label=label, color=color, linewidth=linewidth)
        return p

def plot_tauspec_output(fname='taulikes.dat', label='taulikes', color=None, linewidth=1, linestyle='-', normalize=False):
    if fname.endswith('npy'):
        data = np.load(fname)
    else:
        data = np.loadtxt(fname)
    data[:, 1] = np.exp(data[:, 1] - data[np.argmax(data[:, 1]), 1])
    if normalize:
        dx = data[1, 0] - data[0, 0]
        area = sum(data[:, 1] * dx)
        data[:, 1] = data[:, 1] / area
    if color is None:
        p = plt.plot(data[:, 0], data[:, 1], label=label, linewidth=linewidth, linestyle=linestyle)
    else:
        p = plt.plot(data[:, 0], data[:, 1], label=label, color=color, linewidth=linewidth, linestyle='-')
    return p

def plot_allreal(prefix='taulikes_allreal', numreal=10, suffix='.dat', nside=32, label='BR+G', color=None, linewidth=1, linestyle='-', normalize=False):
    if numreal == 10:
        for i in range(numreal):
            if nside == 32:
                fname = prefix + '_real%02d' % (i + 1) + suffix
            elif nside == 16:
                fname = prefix + '_real%02d' % (i + 1) + suffix
            plt.subplot(2, 5, i + 1)
            plot_output(fname, label=label, color=color, linewidth=linewidth, linestyle=linestyle, normalize=normalize)
            plt.xlim((0.05, 0.12))

def plot_allloris(nside=32, label='Loris', color=None, linewidth=1, normalize=False):
    for i in range(10):
        print i
        plt.subplot(2, 5, i + 1)
        plot_loris(i, color=color, label='Loris', linewidth=linewidth, normalize=normalize)

def plot_allreal_convergence(prefix='taulikes_allreal', numreal=10, 
                             suffix='.dat'):
    if numreal == 10:
        for i in range(numreal):
            plt.subplot(2, 5, i + 1)
            fname = prefix + '_pt1_real%02d' % (i + 1) + suffix
            plot_output(fname, label='BR+G_pt1')
            fname = prefix + '_pt2_real%02d' % (i + 1) + suffix
            plot_output(fname, label='BR+G_pt2')
            plot_loris(i)
            plt.xlim((0.05, 0.12))

def plot_tau_sampler_histogram(fname, normed=True, firstval=0.0100, interval=0.0015, nvals=194, chains=None, label=None, histtype='bar'):
    hdulist = pyfits.open(fname)
    if chains is None:
        dat = hdulist[0].data.flatten()
    else:
        dat = hdulist[0].data[:, chains].flatten()
    hdulist.close()
    start = firstval - interval / 2
    bins = np.zeros(nvals + 1)
    for i in range(nvals + 1):
        bins[i] = start + interval * i
    plt.hist(dat, bins=bins, normed=normed, label=label, histtype=histtype)

def plot_all_tausamples_vs_pixel_based(fname):
    tau= 505
    dtau = 105
    for i in range(10):
        plt.subplot(2, 5, i + 1)
        plot_tau_sampler_histogram(fname, chains=(i,))
#        plot_tauspec_output(fname='/mn/svati/d1/eirikgje/data/vault/tau_likelihood_data/tau_slices/taulikes_dx7_ns16_tau%04d_quietbf_lmax15_EEonly_taucls_WMAP7_pol_mask.dat' % tau, linewidth=2, normalize=True)
        plot_tauspec_output(fname='/mn/svati/d1/eirikgje/data/vault/tau_likelihood_data/tau_slices/taulikes_ctp3_ns16_real%02d_quietbf_lmax15_EEonly.dat' % (i + 1), linewidth=2, normalize=True)
        plt.xlim((0.05, 0.12))
        tau += dtau

def plot_contours(fname, label=None, colors=None, islog=False, linestyle='-'):
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
    lnL = -2*(lnL - np.max(lnL))
   
    my_levels = np.array([0.1, 2.3, 6.17, 11.8])
    
    a = plt.contour(Q_args, N_args, lnL, my_levels, colors=colors, linestyles=linestyle)

    if label is not None:
        plt.plot(Q_args[numpoint[0]/2], N_args[numpoint[1]/2], linestyle, color=colors, label=label)

    return a

def plot_splined_contours(fname, label=None, colors=None, islog=False, linestyle='-', num_spline_points=100):
    data = np.loadtxt(fname, skiprows=1)
    if islog:
        lnL = data[:, 2]
    else:
        lnL = np.log(data[:, 2])
    lnL = -2* (lnL - np.max(lnL))
    #tck = interpolate.bisplrep(data[:, 0], data[:, 1], lnL, s=0)
    tck = interpolate.bisplrep(data[:, 0], data[:, 1], lnL, s=1)
    num_spline_points = complex(0, num_spline_points)
    Q_args, N_args = np.mgrid[data[0, 0]:data[-1, 0]:num_spline_points, data[0, 1]:data[-1, 1]:num_spline_points]
    lnL_splined = interpolate.bisplev(Q_args[:, 0], N_args[0, :], tck)
    
    my_levels = np.array([0.1, 2.3, 6.17, 11.8])
    
    a = plt.contour(Q_args, N_args, lnL_splined, my_levels, colors=colors, linestyles=linestyle)

    if label is not None:
        plt.plot(Q_args[int(abs(num_spline_points))/2, 0], N_args[0, int(abs(num_spline_points))/2], linestyle, color=colors, label=label)

    return a

def plot_errbars_from_sigma_sample_marginals(ls, sigmas, burnin=0, spec=0, lmax=20, sample_fraction=0.68, label=None, color=None, shift=0):
    if isinstance(sigmas, str):
        sigmas = pyfits.open(sigmas)[0].data
    numsamps = np.sum(sigmas[0, :, 0, 0]) - len(sigmas[0, :, 0, 0]) * burnin
    nsigmas = np.zeros(np.append(numsamps, sigmas.shape[2:]))
    j = 0
    for i in range(len(sigmas[0, :, 0, 0])):
        jprev = j
        j += sigmas[0, i, 0, 0] - burnin
        nsigmas[jprev:j] = sigmas[burnin + 1:sigmas[0, i, 0, 0] + 1, i]
#    sigmas = np.reshape(sigmas[1:], ((sigmas.shape[0] - 1) * sigmas.shape[1], 1, sigmas.shape[2], sigmas.shape[3]))
    res = np.zeros((len(ls), 3))
    i = 0
    print sigmas.shape
    for l in ls:
        samps = nsigmas[:, spec, l]
        x = np.linspace(min(samps), max(samps), 100)
        ml, upper, lower = gen_utils.calculate_ml_and_asymm_errorbars_from_samples(samps, x, sample_fraction, smooth=False)
        res[i, 0] = ml
        res[i, 1] = ml - lower
        res[i, 2] = upper - ml
        i += 1

    ls = ls + shift
    plt.scatter(ls, res[:, 0], color=color, label=label)
    plt.errorbar(ls, res[:, 0], res[:, 1:].T, ecolor=color, fmt=None)

def plot_ml_powspec_with_band(sigmas, lmax=50, sample_fraction=0.68, label=None, color=None, spec=0, burnin=0):
    if isinstance(sigmas, str):
        sigmas = pyfits.open(sigmas)[0].data
    numsamps = np.sum(sigmas[0, :, 0, 0]) - len(sigmas[0, :, 0, 0]) * burnin
    nsigmas = np.zeros(np.append(numsamps, sigmas.shape[2:]))
    j = 0
    for i in range(len(sigmas[0, :, 0, 0])):
        jprev = j
        j += sigmas[0, i, 0, 0] - burnin
        nsigmas[jprev:j] = sigmas[burnin + 1:sigmas[0, i, 0, 0] + 1, i]
#    sigmas = np.reshape(sigmas[1:], ((sigmas.shape[0] - 1) * sigmas.shape[1], 1, sigmas.shape[2], sigmas.shape[3]))
    ls = np.arange(2, lmax + 1)
    res = np.zeros((len(ls), 3))
    i = 0
    for l in ls:
        samps = nsigmas[:, spec, l]
        x = np.linspace(min(samps), max(samps), 100)
        ml, upper, lower = gen_utils.calculate_ml_and_asymm_errorbars_from_samples(samps, x, sample_fraction, smooth=False)
        res[i, 0] = ml
#        res[i, 1] = ml - lower
#        res[i, 2] = upper - ml
        res[i, 1] = lower
        res[i, 2] = upper
        i += 1
#    plt.plot(ls, res[:, 0], color=color, label=label)
    plt.plot(ls, res[:, 1], color=color, label=label)
    plt.plot(ls, res[:, 2], color=color)
