from __future__ import division
import numpy as np
import os
import pyfits

def gr(chains):
    chain_means = np.zeros(len(chains))
    chain_vars = np.zeros(len(chains))
    b = 0
    w = 0
    for i in range(len(chains)):
        chain_means[i] = np.mean(chains[i])
        chain_vars[i] = np.sum((chains[i] - chain_means[i])**2)
    chain_vars = chain_vars/(len(chains[0]) - 1)
    w = np.mean(chain_vars)
    b = np.sum((chain_means - np.mean(chain_means))**2)
    b = b * len(chains[0]) / (len(chains) - 1)
    var = (len(chains[0]) - 1) / len(chains[0]) * w + 1 / len(chains[0]) * b
    return np.sqrt(var/w)

def calc_gr_all_fits(chain, numchains, lmax=1024, spec=1, burnin=None):
    if isinstance(chain, str):
        chain = pyfits.open(chain)[0].data
    chain = chain[:, :, spec-1]
    ret = np.zeros(lmax)
    for i in range(lmax):
        ret[i] = gr(chain[:, :, i].T)

    return ret

def calc_gr_raw(l=2, ps='cls', lmax=90, spec=1, chains=None, burnin=None):
    if burnin is None:
        burnin = 0
    dat = []
    if ps == 'cls':
        ind = l - 2
    elif ps == 'sigma':
        ind = l - 1 + lmax - 2
    else:
        raise ValueError("ps has unknown value")
    minnum = 1000000
    for chain in chains:
        filenum = 1
        cont = True
        filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
        cont = os.path.exists(filename)
        tdat = []
        while cont:
            temp = np.loadtxt(filename)
            tdat.append(temp[ind, spec])
            filenum += 1
            filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
            cont = os.path.exists(filename)
        if filenum - 1 < minnum:
            minnum = filenum - 1
        dat.append(tdat)
    ndat = np.zeros((len(dat), minnum))
    for i in range(len(dat)):
        for j in range(minnum):
            ndat[i,j] = dat[i][j]
    return gr(ndat)

def calc_gr_all_multipoles_raw(lmax, lmin=2, ps='cls', spec=1, chains=None, burnin =None):
    if burnin is None:
        burnin = 0
    if ps == 'cls':
        inds = np.arange(lmin-2, lmax-1)
    elif ps == 'sigma':
        inds = np.arange(lmin - 2 + lmax - 1, 2 * lmax - 1)
    else:
        raise ValueError("ps has unknown value")
    dat = []
    minnum = 1000000
    for chain in chains:
        filenum = 1
        cont = True
        filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
        cont = os.path.exists(filename)
        tdat = []
        while cont:
            ttdat = np.zeros(len(inds))
            temp = np.loadtxt(filename)
            ttdat[:] = (temp[inds, spec])
            tdat.append(ttdat)
            filenum += 1
            filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
            cont = os.path.exists(filename)
        if filenum - 1 < minnum:
            minnum = filenum - 1
        dat.append(tdat)
    ndat = np.zeros((len(dat), minnum, len(inds)))
    for i in range(len(dat)):
        for j in range(minnum):
            ndat[i, j] = dat[i][j]

    ndat = np.rollaxis(ndat, 2)
    gelm = []
    for l in range(len(inds)):
        gelm.append(gr(ndat[l]))
    return gelm

def calc_gr_all_multipoles_all_chains_raw(lmax=90, lmin=2, ps='cls', spec=1, burnin=None, chain_start=1, chain_end=100, chain_intervals=10):
    if burnin is None:
        burnin = 0
    if ps == 'cls':
        inds = np.arange(lmin-2, lmax-1)
    elif ps == 'sigma':
        inds = np.arange(lmin - 2 + lmax - 1, 2 * lmax - 2)
    else:
        raise ValueError("ps has unknown value")
    minnum = 1000000
    
    totgelm = []
    for num in range(chain_start, chain_end, chain_intervals):
        dat = []
        chains = range(num, num + chain_intervals)
        for chain in chains:
            filenum = 1
            cont = True
            filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
            cont = os.path.exists(filename)
            tdat = []
            while cont:
                ttdat = np.zeros(len(inds))
                temp = np.loadtxt(filename)
                ttdat[:] = (temp[inds, spec])
                tdat.append(ttdat)
                filenum += 1
                filename = 'cl_c%04d_k%05d.dat' % (chain, filenum + burnin)
                cont = os.path.exists(filename)
            if filenum - 1 < minnum:
                minnum = filenum - 1
            print minnum
            dat.append(tdat)
        ndat = np.zeros((len(dat), minnum, len(inds)))
        for i in range(len(dat)):
            for j in range(minnum):
                ndat[i, j] = dat[i][j]
    
        ndat = np.rollaxis(ndat, 2)
        gelm = []
        print ndat.shape
        for l in range(len(inds)):
            gelm.append(gr(ndat[l]))
        totgelm.append(gelm)
    return totgelm

