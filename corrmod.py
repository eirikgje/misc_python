from __future__ import division
import numpy as np
import os

def corr(chain):
    c = []
    mean = np.mean(chain)
    var = np.var(chain)
    for i in range(len(chain) - 1):
        c.append(np.mean((chain[:len(chain)-(i + 1)] - mean) * \
                (chain[i+1:] - mean))/var)
    return c

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
