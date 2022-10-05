#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
#
# Originally contributed by Mawson Sammons
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Constraints for optimizers to evaluate shark models against observations
"""

import os
import sys

import common
import numpy as np
#import smf
import re

#sys.path.insert(0, '../output')
import routines as r

GyrToYr = 1e9

#######################
# Binning configuration
mlow = 5
mupp = 14
dm = 0.2
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

mlow2 = 5
mupp2 = 14
dm2 = 0.3
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

ssfrlow = -6
ssfrupp = 4
dssfr = 0.2
ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)


# These are two easily create variables of these different shapes without
# actually storing a reference ourselves; we don't need it
zeros1 = lambda: np.zeros(shape=(1, 3, len(xmf)))
zeros2 = lambda: np.zeros(shape=(1, 3, len(xmf2)))
zeros3 = lambda: np.zeros(shape=(1, len(mbins)))
zeros4 = lambda: np.empty(shape=(1), dtype=np.bool_)
zeros5 = lambda: np.zeros(shape=(1, len(ssfrbins)))

class Constraint(object):
    """Base classes for constraint objects"""

    def __init__(self):
        self.redshift_table = None
        self.weight = 1
        self.rel_weight = 1

    def _load_model_data(self, modeldir, subvols):

        if  len(subvols) > 1:
            subvols = ["multiple_batches"]

        # Histograms we are interested in
        hist_smf = zeros3()
        hist_HImf = zeros3()

#        fields = {
#            'galaxies': (
#                'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
#                'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
#                'matom_bulge', 'mmol_bulge', 'mgas_bulge',
#                'mgas_metals_disk', 'mgas_metals_bulge',
#                'mstars_metals_disk', 'mstars_metals_bulge', 'type',
#                'mvir_hosthalo', 'rstar_bulge')
#        }

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
        # hard coding stuff here that should be generalised
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
        fields = ['StellarMass', 'DiscHI', 'LenMax', 'DiscStars', 'MergerBulgeMass', 'InstabilityBulgeMass', 'IntraClusterStars', 'LocalIGS']
        files = [2,3]
        Nage = 30
        G = r.darksage_snap(modeldir+'model_z0.000', files, Nannuli=30, Nage=Nage, fields=fields)
#        h0 = 0.6751
#        vol = (75.0/h0)**3
        h0 = 0.6774
        Omega0 = 0.3089
        vol = (205.0/h0)**3 * (1.0*len(files)/128.)
        age_alist_file = '/Users/adam/Illustris/alist_TNG.txt'
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
        

#        for index, z in enumerate(self.z):
#            hdf5_data = common.read_data(modeldir, self.redshift_table[z], fields, subvols)
#            h0 = hdf5_data[0]
#            smf.prepare_data(hdf5_data, index, hist_smf, zeros3(), zeros3(),
#                             zeros3(), zeros3(), hist_HImf, zeros3(), zeros3(),
#                             zeros3(), zeros3(), zeros3(), zeros1(), zeros1(),
#                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(),
#                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(),
#                             zeros1(), zeros4(), zeros4(), zeros2(), zeros5(),
#                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(),
#                             zeros1())

        logSM = np.log10(G['StellarMass']*1e10/h0)
        logSM[~np.isfinite(logSM)] = -20
        hist_smf, _ = np.histogram(logSM, bins=mbins)
        hist_smf = hist_smf / (dm * vol)
        hist_smf = hist_smf[np.newaxis]

        logHIM = np.log10(np.sum(G['DiscHI'], axis=1)*1e10/h0)
        logHIM[~np.isfinite(logHIM)] = -20
        hist_HImf, _ = np.histogram(logHIM, bins=mbins)
        hist_HImf = hist_HImf / (dm * vol)
        hist_HImf = hist_HImf[np.newaxis]
        
        # sum mass from age bins of all components
        StarsByAge = np.zeros(Nage)
        for k in range(Nage):
            StarsByAge[k] += np.sum(G['DiscStars'][:,:,k])
            StarsByAge[k] += np.sum(G['MergerBulgeMass'][:,k])
            StarsByAge[k] += np.sum(G['InstabilityBulgeMass'][:,k])
            StarsByAge[k] += np.sum(G['IntraClusterStars'][:,k])
            StarsByAge[k] += np.sum(G['LocalIGS'][:,k])
            
        # get the edges of the age bins
        alist = np.loadtxt(age_alist_file)
        if Nage>=len(alist)-1:
            alist[::-1]
            RedshiftBinEdge = 1./ alist - 1.
        else:
            indices_float = np.arange(Nage+1) * (len(alist)-1.0) / Nage
            indices = indices_float.astype(np.int32)
            alist = alist[indices][::-1]
            RedshiftBinEdge = 1./ alist - 1.
        TimeBinEdge = np.array([r.z2tL(redshift, h0, Omega0, 1.0-Omega0) for redshift in RedshiftBinEdge]) # look-back time [Gyr]
        dT = np.diff(TimeBinEdge) # time step for each bin
        m, lifetime, returned_mass_fraction_integrated, ncum_SN = r.return_fraction_and_SN_ChabrierIMF()
        eff_recycle = np.interp(TimeBinEdge[:-1] + 0.5*dT, lifetime[::-1], returned_mass_fraction_integrated[::-1])
        SFRbyAge = StarsByAge*1e10/h / (dT*1e9) / (1.-eff_recycle)
        SFRD_Age = np.log10(np.append(SFRbyAge[0],SFRbyAge)/vol)


        #########################
        # take logs
        ind = np.where(hist_smf > 0.)
        hist_smf[ind] = np.log10(hist_smf[ind])
        ind = np.where(hist_HImf > 0.)
        hist_HImf[ind] = np.log10(hist_HImf[ind])

        return h0, Omega0, hist_smf, hist_HImf, TimeBinEdge, SFRD_Age


    def load_observation(self, *args, **kwargs):
        obsdir = os.path.normpath(os.path.abspath(os.path.join(__file__, '..')))#, '..', 'data')))
#        obsdir = os.path.normpath(os.path.abspath(__file__))
        return common.load_observation(obsdir, *args, **kwargs)

    def _get_raw_data(self, modeldir, subvols):
        """Gets the model and observational data for further analysis.
        The model data is interpolated to match the observation's X values."""

        h0, Omega0, hist_smf, hist_HImf, TimeBinEdge, SFRD_Age = self._load_model_data(modeldir, subvols)
        x_obs, y_obs, y_dn, y_up = self.get_obs_x_y_err(h0, Omega0)
        x_mod, y_mod = self.get_model_x_y(hist_smf, hist_HImf, TimeBinEdge, SFRD_Age)
        return x_obs, y_obs, y_dn, y_up, x_mod, y_mod

    def get_data(self, modeldir, subvols):

        x_obs, y_obs, y_dn, y_up, x_mod, y_mod = self._get_raw_data(modeldir, subvols)

        # Linearly interpolate model Y values respect to the observations'
        # X values, and only take those within the domain.
        # We also consider the biggest relative error as "the" error, in case
        # they are different
        y_mod = np.interp(x_obs, x_mod, y_mod)
        ind = np.where((x_obs >= self.domain[0]) & (x_obs <= self.domain[1]))
        err = np.maximum(np.abs(y_dn[ind]), np.abs(y_up[ind]))
        print('in get_data:')
        print('obs x:', x_obs[ind])
        print('obs y:', y_obs[ind])
        print('mod y:', y_mod[ind])
        return y_obs[ind], y_mod[ind], err

    def __str__(self):
        s = '%s, low=%.1f, up=%.1f, weight=%.2f, rel_weight=%.2f'
        args = self.__class__.__name__, self.domain[0], self.domain[1], self.weight, self.rel_weight
        return s % args


class HIMF(Constraint):
    """The HI Mass Function constraint"""

    domain = (7, 12)
    z = [0]

    def get_obs_x_y_err(self, h0, _):
        # Load Jones18 data and correct data for their choice of cosmology
        hobs = 0.7
        log_mHI, phiHI, delta_phiHI = self.load_observation('Jones2018_HIMF.dat', cols=[0,1,2])
        x_obs = log_mHI + 2.0 * np.log10(hobs/h0)
        y_obs = np.log10(phiHI) + 3.0 * np.log10(h0/hobs)
        y_dn = np.log10(phiHI - delta_phiHI) + 3.0 * np.log10(h0/hobs)
        y_up = np.log10(phiHI + delta_phiHI) + 3.0 * np.log10(h0/hobs)
        return x_obs, y_obs, y_dn, y_up

    def get_model_x_y(self, _, hist_HImf, _2, _3):
        y = hist_HImf[0]
        ind = np.where(y < 0.)
        return xmf[ind], y[ind]

class SMF(Constraint):
    """Common logic for SMF constraints"""

    domain = (8, 13)

    def get_model_x_y(self, hist_smf, _, _2, _3):
        y = hist_smf[0,:]
        ind = np.where(y < 0.)
        return xmf[ind], y[ind]

class SMF_z0(SMF):
    """The SMF constraint at z=0"""

    z = [0]

    def get_obs_x_y_err(self, h0, Omega0):
        # Load data from Driver et al. (2022)
        logm, logphi, dlogphi = self.load_observation('GAMA_SMF.dat', cols=[0,1,2])
        cosmology_correction = np.log10( r.comoving_distance(0.079, 100*h0, 0, Omega0, 1.0-Omega0) / r.comoving_distance(0.079, 70.0, 0, 0.3, 0.7) )
        x_obs = logm + 2.0 * cosmology_correction 
        y_obs = logphi + 3.0 * cosmology_correction + 0.0769 # last factor accounts for average under-density of GAMA
        y_dn = y_obs - dlogphi
        y_up = y_obs + dlogphi
        return x_obs, y_obs, y_dn, y_up


class SMF_z1(SMF):
    """The SMF constraint at z=1"""

    z = [1]

    def get_obs_x_y_err(self, _, _2):

        # Wright et al. (2018, several reshifts). Assumes Chabrier IMF.
        zD17, lmD17, pD17, dp_dn_D17, dp_up_D17 = self.load_observation('mf/SMF/Wright18_CombinedSMF.dat', cols=[0,1,2,3,4])
        hobs = 0.7
        pD17 = pD17 - 3.0 * np.log10(hobs)
        lmD17 = lmD17 - np.log10(hobs)
        in_redshift = np.where(zD17 == 1)
        x_obs = lmD17[in_redshift]
        y_obs = pD17[in_redshift]
        y_dn = dp_dn_D17[in_redshift]
        y_up = dp_up_D17[in_redshift]

        return x_obs, y_obs, y_dn, y_up
                              
class CSFRDH(Constraint):

    z = [0]
    domain = (0, 14) # look-back time in Gyr
    
    def get_obs_x_y_err(self, h0, Omega0):
        zmin, zmax, logSFRD, err1, err2, err3 = self.load_observation('Driver_SFRD.dat', cols=[1,2,3,5,6,7])
        Np = len(logSFRD)
        x_obs = r.z2tL(0.5*(zmin+zmax), 100*h0, Omega0,  1.0-Omega0)
        y_obs = np.zeros(Np)
        for i in range(Np):
            y_obs[i] = logSFRD[i] + np.log10(pow(r.comoving_distance(zmax[i], 70.0, 0, 0.3, 0.7), 3.0) - pow(r.comoving_distance(zmin[i], 70.0, 0, 0.3, 0.7), 3.0)) - np.log10(pow(r.comoving_distance(zmax[i], 100*h0, 0, Omega0,  1.0-Omega0), 3.0) - pow(r.comoving_distance(zmin[i], 100*h0, 0, Omega0,  1.0-Omega0), 3.0))
            
        err_total = err1 + err2 + err3
        y_dn = y_obs - err_total
        y_up = y_obs + err_total
        
        return x_obs, y_obs, y_dn, d_up
        
    def get_model_x_y(self, _, _2, TimeBinEdge, SFRD_Age):
        return 0.5*(TimeBinEdge[1:]+TimeBinEdge[:-1]), SFRD_Age
        
    

_constraint_re = re.compile((r'([0-9_a-zA-Z]+)' # name
                              '(?:\(([0-9\.]+)-([0-9\.]+)\))?' # domain boundaries
                              '(?:\*([0-9\.]+))?')) # weight
def parse(spec):
    """Parses a comma-separated string of constraint names into a list of
    Constraint objects. Specific domain values can be specified in `spec`"""

    _constraints = {
        'HIMF': HIMF,
        'SMF_z0': SMF_z0,
        'SMF_z1': SMF_z1,
        'CSFRDH': CSFRDH
    }

    def _parse(s):
        m = _constraint_re.match(s)
        if not m or m.group(1) not in _constraints:
            raise ValueError('Constraint does not specify a valid constraint: %s' % s)
        c = _constraints[m.group(1)]()
        if m.group(2):
            dn, up = float(m.group(2)), float(m.group(3))
            if dn < c.domain[0]:
                raise ValueError('Constraint low boundary is lower than lowest value possible (%f < %f)' % (dn, c.domain[0]))
            if up > c.domain[1]:
                raise ValueError('Constraint up boundary is higher than lowest value possible (%f > %f)' % (up, c.domain[1]))
            c.domain = (dn, up)
        if m.group(4):
            c.weight = float(m.group(4))
        return c

    constraints = [_parse(s) for s in spec.split(',')]
    total_weight = sum([c.weight for c in constraints])
    for c in constraints:
        c.rel_weight = c.weight / total_weight
    return constraints
