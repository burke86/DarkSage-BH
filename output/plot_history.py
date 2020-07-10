# Plots that require multiple snapshots of Dark Sage. Potentially usable for calibration. First version just does the Madau plot.

from pylab import *
import os
import routines as r
import random

# Warnings are annoying
import warnings
warnings.filterwarnings("ignore")


###### USER NEEDS TO SET THESE THINGS ######
indir = '/Users/adam/DarkSage_runs/571j/' # directory where the Dark Sage data are
sim = 0 # which simulation Dark Sage has been run on -- if it's new, you will need to set its defaults below.
#   0 = Mini Millennium, 1 = Full Millennium, 2 = SMDPL

Nannuli = 30 # number of annuli used for discs in Dark Sage
Nage = 20 # number of age bins used for stars -- not advised to use a run with this for calibration
RecycleFraction = 0.43 # Only needed for comparing stellar-age based SFR with raw SFR
###### ============================== ######



##### SET PLOTTING DEFAULTS #####
fsize = 26
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize-4, 'xtick.top': True, 'ytick.right': True})

NpartMed = 100 # minimum number of particles for finding relevant medians for minima on plots

outdir = indir+'plots/' # where the plots will be saved
if not os.path.exists(outdir): os.makedirs(outdir)
######  =================== #####



##### SEARCH DIRECTORY FOR SNAPSHOTS AND FILENUMBERS PRESENT #####
filenames = os.listdir(indir)
redshiftstr = np.array([], dtype=str)
filenumbers = np.array([], dtype=np.int32)
fpre = None
for f in filenames:
    w1 = f.find('_z')
    if w1==-1: continue
    w2 = f.find('_', w1+3)
    zstr = f[w1+2:w2]
    fno = int(f[w2+1:])
    if zstr not in redshiftstr: redshiftstr = np.append(redshiftstr, zstr)
    if fno not in filenumbers: filenumbers = np.append(filenumbers, fno)
    if fpre is None: fpre = f[:w1+2]
redshifts = redshiftstr.astype(np.float32)
args = np.argsort(redshifts) # want redshifts to be listed in order
redshiftstr = redshiftstr[args]
redshifts = redshifts[args]
##### ====================================================== #####



##### SIMULATION DEFAULTS #####
if sim==0:
    h = 0.73
    vol = (62.5/h)**3 * len(filenumbers)/8. # comoving volume of the (part of the) simulation
    Omega_M, Omega_L = 0.25, 0.75
elif sim==1:
    h = 0.73
    vol = (500.0/h)**3 * len(filenumbers)/512.
elif sim==2:
    h = 0.6777
    vol = (400.0/h)**3 * len(filenumbers)/1000.
# add here 'elif sim==3:' etc for a new simulation
else:
    print('Please specify a valid simulation.  You may need to add its defaults to this code.')
    quit()
######  ================= #####



##### READ DARK SAGE DATA AND BUILD RELEVANT ARRAYS #####
Nsnap = len(redshifts)
SFRD = np.zeros(Nsnap)
SFRD_resolved = np.zeros(Nsnap)

for i in range(Nsnap):
    G = r.darksage_snap(indir+fpre+redshiftstr[i], filenumbers, Nannuli=Nannuli, Nage=Nage)
    if i==0 and Nage>1: G0 = G # save the lowest-z snap to compare history reconstruction from age bins
    SFRD[i] = np.sum(G['SfrFromH2']+G['SfrInstab']+G['SfrMergeBurst']) / vol
    SFRD_resolved[i] = np.sum((G['SfrFromH2']+G['SfrInstab']+G['SfrMergeBurst'])[G['LenMax']>=100]) / vol
##### ============================================= #####




##### PLOT 1: MADAU-LILLY DIAGRAM (UNIVERSAL STAR FORMATION RATE DENSITY HISTORY) #####
try:
    fig, ax  = plt.subplots(1, 1)
    plt.plot(1+redshifts, np.log10(SFRD), 'k--', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p}\!\geq\!20$')
    plt.plot(1+redshifts, np.log10(SFRD_resolved), 'k.-', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p,max}\!\geq\!100$')
    
    if Nage>1: # check consistency from stellar-age bins
#        print 'All Close', np.allclose(np.sum(G0['DiscStars'],axis=(1,2)) + np.sum(G0['InstabilityBulgeMass'],axis=1) + np.sum(G0['MergerBulgeMass'],axis=1), G0['StellarMass'])
    
        RedshiftBinEdge = np.append(np.append(0, [0.007*1.47**i for i in range(Nage-1)]), 1000.) # defined within Dark Sage hard code
        TimeBinEdge = np.array([r.z2tL(z,h,Omega_M,Omega_L) for z in RedshiftBinEdge]) # look-back time [Gyr]
        dT = np.diff(TimeBinEdge)*1e9 # time step for each bin [yr]
        
        # sum mass from age bins of all components
        StarsByAge = np.zeros(Nage)
        for k in range(Nage):
            StarsByAge[k] += np.sum(G0['DiscStars'][:,:,k])
            StarsByAge[k] += np.sum(G0['MergerBulgeMass'][:,k])
            StarsByAge[k] += np.sum(G0['InstabilityBulgeMass'][:,k])
            StarsByAge[k] += np.sum(G0['IntraClusterStars'][:,k])
        SFRbyAge = StarsByAge*1e10/h/dT/(1.-RecycleFraction)
        SFRD_Age = np.log10(np.append(SFRbyAge[0],SFRbyAge)/vol)
        plt.step(1+RedshiftBinEdge, SFRD_Age, color='silver', label=r'{\sc Dark Sage} age recon.')

    
    r.SFRD_obs(h, plus=1)
    plt.xlabel(r'Redshift')
    plt.ylabel(r'$\log_{10}\left( \bar{\rho}_{\rm SFR}~[{\rm M}_{\odot}\, {\rm yr}^{-1}\, {\rm cMpc}^{-3}] \right)$')
    plt.xscale('log')
    plt.axis([1,9,-2.8,0.5])
    plt.minorticks_off()
    plt.xticks(range(1,10), (str(i) for i in range(9)))
    plt.legend(loc='best', frameon=False, ncol=2)
    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'H1-SFRDH', xsize=768, ysize=400)
except Exception as excptn:
    print('Unexpected issue with plot H1: {0}'.format(excptn))
##### =========================================================================== #####



##### PLOT 2: SFRDH and SMH BREAKDOWN BY z=0 GALAXY MASS #####
if Nage>1: # currently just built for stellar-age bins -- can definitely be extended to use other snapshots (and should be as a consistency check!)
    Mstar_bins = 10**np.array([8.5, 9.5, 10.5]) * h * 1e-10
    c = ['y', 'g', 'c', 'b']
    Nbins = len(Mstar_bins)+1
    
    SFRDH = SFRbyAge/vol
    SMDH = np.cumsum(StarsByAge[::-1])[::-1]*1e10/h / vol
    TimeBinCentre = 0.5*(TimeBinEdge[1:] + TimeBinEdge[:-1])
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False)

    ax1.plot(TimeBinCentre, SFRDH, '-', color='grey', lw=3)
    ax2.plot(TimeBinEdge, np.append(SMDH,0), '-', color='grey', lw=3)
    
    t_LB = np.array([r.z2tL(z) for z in redshifts])
    ax1.plot(t_LB, SFRD, '.-', color='k', lw=2)

    
    for i in xrange(Nbins):
        if i==0:
            f = (G0['StellarMass']<Mstar_bins[0])
        elif i==Nbins-1:
            f = (G0['StellarMass']>=Mstar_bins[-1])
        else:
            f = (G0['StellarMass']>=Mstar_bins[i-1]) * (G0['StellarMass']<Mstar_bins[i])
            
        StarsByAge = np.zeros(Nage)
        for k in range(Nage):
            StarsByAge[k] += np.sum(G0['DiscStars'][:,:,k][f])
            StarsByAge[k] += np.sum(G0['MergerBulgeMass'][:,k][f])
            StarsByAge[k] += np.sum(G0['InstabilityBulgeMass'][:,k][f])
            StarsByAge[k] += np.sum(G0['IntraClusterStars'][:,k][f])
        SFRDH = StarsByAge*1e10/h/dT/(1.-RecycleFraction) / vol
        SMDH = np.cumsum(StarsByAge[::-1])[::-1]*1e10/h / vol
        
        ax1.plot(TimeBinCentre, SFRDH, '-', color=c[i], lw=2)
        ax2.plot(TimeBinEdge, np.append(SMDH,0), '-', color=c[i], lw=2)

    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.axis([0,14,2e-4,1e-1])
    ax2.set_ylim(5e4,5e8)
    
    ax1.set_ylabel(r'SFRD [M$_\odot$\,yr$^{-1}$\,cMpc$^{-3}$]')
    ax2.set_ylabel(r'SMD [M$_\odot$\,cMpc$^{-3}$]')
    ax2.set_xlabel(r'Lookback time [Gyr]')
    
    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'H2-SFRDH+SMDH', xsize=768, ysize=800)
##### ================================================== #####
