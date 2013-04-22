#Just a module with functions that easily perform tasks needed to run the slicer. Will be updated ad-hoc.

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import shlex
import pyfits
import sys

def cp(fname1, fname2):
    subprocess.call(shlex.split('cp ' + fname1 + ' ' + fname2))

def rm(fname):
    subprocess.call(shlex.split('rm ' + fname))

class Pluginfo(object):
#    def __init__(self, lmin=2, lmax=30, numranges=2, lranges=(2, 9, 30),
#                 name='BR_MOD', sigmafile=None, clfile=None,
#                 relpath=None, rangetypes = ('TT_TE_EE', 'TT'),
#                 firstsamp=1, lastsamp=10000, spectras='TE_EE'):
    def __init__(self, name, lmin=2, lmax=30, numranges=None, lranges=None,
                 sigmafile=None, clfile=None,
                 relpath=None, rangetypes = None,
                 firstsamp=None, lastsamp=None, spectra=None, chatfile=None,
                 sigma_t_uniform=None, firstchain=None, lastchain=None,
                 lmax_like=None, lmin_like=None, comm_lmin=None, comm_lmax=None,
                 l_thresh=None, cmbfile=None):
        if relpath is None:
            self.relpath = 'datasets/' + name.lower() + '/'
        else:
            self.relpath = relpath
        self.lmin = lmin
        self.lmax = lmax
        if lmin_like is None:
            self.lmin_like=lmin
        else:
            self.lmin_like=lmin_like
        if lmax_like is None:
            self.lmax_like=lmax
        else:
            self.lmax_like=lmax_like
        self.comm_lmin=comm_lmin
        self.comm_lmax=comm_lmax
        self.l_thresh = l_thresh
        self.numranges = numranges
        self.lranges = lranges
        self.name = name
        self.sigmafile = sigmafile
        self.clfile = clfile
        self.rangetypes = rangetypes
        self.firstsamp = firstsamp
        self.spectra = spectra
        self.chatfile = chatfile
        self.sigma_t_uniform = sigma_t_uniform
        self.firstchain = firstchain
        self.lastchain = lastchain
        self.cmbfile = cmbfile
        if lastsamp is None:
            if sigmafile is not None:
                self.lastsamp = pyfits.open(self.relpath + sigmafile)[0].header['NUMSAMP']
            elif clfile is not None:
                self.lastsamp = pyfits.open(self.relpath + clfile)[0].header['NUMSAMP']
            else:
                self.lastsamp = None
        else:
            self.lastsamp = lastsamp

def update_files(pluginlist, app='tauspec', lmax=90, 
             path = '/mn/svati/d1/eirikgje/data/spline_project/slicer_data/ctp3_trieste/'):
    if app == 'tauspec':
        fname1 = path + 'param_tauspec.txt'
        fname2 = path + 'param_tauspec.txt.bak'
        cp(fname1, fname2)
        rm(fname1)
        rf = open(fname2, 'r')
        wf = open(fname1, 'w')
        for line in rf:
            if line.startswith('LMAX'):
                wf.write('LMAX = ' + str(lmax) + '\n')
            else:
                wf.write(line)
        rf.close()
        wf.close()

    fname1 = path + 'fplc_config.txt'
    rm(fname1)
    wf = open(fname1, 'w')
    for plugin in pluginlist:
        wf.write(plugin.name + '   ' + str(plugin.lmin_like) + '   ' + 
                str(plugin.lmax_like) + '   T    "' + plugin.relpath[:-1] + '"\n')
    wf.close()
    for plugin in pluginlist:
        fname1 = path + plugin.relpath + 'par.txt'
        fname2 = path + plugin.relpath + 'par.txt.bak'
        cp(fname1, fname2)
        rm(fname1)
        rf = open(fname2, 'r')
        wf = open(fname1, 'w')
        for line in rf:
            if line.startswith('COMMANDER_SIGMAFILE'):
                if plugin.sigmafile is not None:
                    wf.write("COMMANDER_SIGMAFILE = '" + plugin.sigmafile + 
                             "'\n")
                    continue
            elif line.startswith('COMMANDER_CLFILE'):
                if plugin.clfile is not None:
                    wf.write("COMMANDER_CLFILE = '" + plugin.clfile + "'\n")
                    continue
            elif line.startswith('NUM_L_RANGES'):
                if plugin.numranges is not None:
                    wf.write('NUM_L_RANGES = ' + str(plugin.numranges) + '\n')
                    continue
            if plugin.numranges is not None:
                written = False
                for i in range(plugin.numranges + 1):
                    if line.startswith('L%02d' % i):
                        wf.write('L%02d = ' % i + str(plugin.lranges[i]) + '\n')
                        written = True
                    if i > 0:
                        if line.startswith('RANGE_TYPE%02d' % i):
                            wf.write("RANGE_TYPE%02d = '" % i  + 
                                     plugin.rangetypes[i-1] + "'\n")
                            written = True
                if written:
                    continue
            if line.startswith('COMMANDER_FIRST_SAMPLE'):
                if plugin.firstsamp is not None:
                    wf.write('COMMANDER_FIRST_SAMPLE = ' + 
                             str(plugin.firstsamp) + '\n')
                    continue
            elif line.startswith('COMMANDER_LAST_SAMPLE'):
                if plugin.lastsamp is not None:
                    wf.write('COMMANDER_LAST_SAMPLE = ' + 
                             str(plugin.lastsamp) + '\n')
                    continue
            elif line.startswith('COMMANDER_FIRST_CHAIN'):
                if plugin.firstchain is not None:
                    wf.write('COMMANDER_FIRST_CHAIN = ' + 
                             str(plugin.firstchain) + '\n')
                    continue
            elif line.startswith('COMMANDER_LAST_CHAIN'):
                if plugin.lastchain is not None:
                    wf.write('COMMANDER_LAST_CHAIN = ' + 
                             str(plugin.lastchain) + '\n')
                    continue
            if line.startswith('L_THRESH'):
                if plugin.l_thresh is not None:
                    wf.write('L_THRESH = ' + 
                             str(plugin.l_thresh) + '\n')
                    continue
            elif line.startswith('LMIN'):
                if plugin.lmin is not None:
                    wf.write('LMIN = ' + str(plugin.lmin) + '\n')
                    continue
            elif line.startswith('LMAX'):
                if plugin.lmax is not None:
                    wf.write('LMAX = ' + str(plugin.lmax) + '\n')
                    continue
            elif line.startswith('COMMANDER_LMIN'):
                if plugin.comm_lmin is not None:
                    wf.write('COMMANDER_LMIN = ' + str(plugin.comm_lmin) + '\n')
                    continue
            elif line.startswith('COMMANDER_LMAX'):
                if plugin.comm_lmax is not None:
                    wf.write('COMMANDER_LMAX = ' + str(plugin.comm_lmax) + '\n')
                    continue
            elif line.startswith('SPECTRA'):
                if plugin.spectra is not None:
                    wf.write("SPECTRA = '" + plugin.spectra + "'\n")
                    continue
            elif line.startswith('CHAT_FILE'):
                if plugin.chatfile is not None:
                    wf.write("CHAT_FILE = '" + plugin.chatfile + "'\n")
                    continue
            elif line.startswith('SIGMA_T_UNIFORM'):
                if plugin.sigma_t_uniform is not None:
                    wf.write("SIGMA_T_UNIFORM = " + str(plugin.sigma_t_uniform)
                             + "\n")
                    continue
            elif line.startswith('CMBFILE'):
                if plugin.cmbfile is not None:
                    wf.write("CMBFILE = '" + plugin.cmbfile + "'"
                             + "\n")
                    continue
            wf.write(line)
        wf.close()
        rf.close()

def run_program(app='tauspec', logfile=None):
    if app == 'tauspec':
        print logfile
        if logfile is None:
            subprocess.call(shlex.split('/mn/svati/u1/eirikgje/src/uio_svn/src/flikelihood/branches/exp_eirik/src/applications/tauspec_eval/tauspec_eval param_tauspec.txt'))
        else:
            f = open(logfile, 'w')
            subprocess.call(shlex.split('/mn/svati/u1/eirikgje/src/uio_svn/src/flikelihood/branches/exp_eirik/src/applications/tauspec_eval/tauspec_eval param_tauspec.txt'), stdout=f)
            f.close()


def savefile(fname):
    if fname.endswith('npy'):
        np.save(fname, np.loadtxt('taulikes.dat'))
    else:
        cp('taulikes.dat', fname)
