# Routines used for reading and plotting data (e.g. used by plot_constraints_z0.py).  Many of these are copied from the arhstevens/Dirty-AstroPy GitHub repository.

from pylab import *


def galdtype_darksage(Nannuli=30):
    floattype = np.float32
    Galdesc_full = [
                    ('Type'                         , np.int32),
                    ('GalaxyIndex'                  , np.int64),
                    ('HaloIndex'                    , np.int32),
                    ('SimulationHaloIndex'          , np.int32),
                    ('TreeIndex'                    , np.int32),
                    ('SnapNum'                      , np.int32),
                    ('CentralGalaxyIndex'           , np.int64),
                    ('CentralMvir'                  , floattype),
                    ('mergeType'                    , np.int32),
                    ('mergeIntoID'                  , np.int32),
                    ('mergeIntoSnapNum'             , np.int32),
                    ('dT'                           , floattype),
                    ('Pos'                          , (floattype, 3)),
                    ('Vel'                          , (floattype, 3)),
                    ('Spin'                         , (floattype, 3)),
                    ('Len'                          , np.int32),
                    ('LenMax'                       , np.int32),
                    ('Mvir'                         , floattype),
                    ('Rvir'                         , floattype),
                    ('Vvir'                         , floattype),
                    ('Vmax'                         , floattype),
                    ('VelDisp'                      , floattype),
                    ('DiscRadii'                    , (floattype, Nannuli+1)), 
                    ('ColdGas'                      , floattype),
                    ('StellarMass'                  , floattype),
                    ('MergerBulgeMass'              , floattype),
                    ('InstabilityBulgeMass'          , floattype),
                    ('HotGas'                       , floattype),
                    ('EjectedMass'                  , floattype),
                    ('BlackHoleMass'                , floattype),
                    ('IntraClusterStars'            , floattype),
                    ('DiscGas'                      , (floattype, Nannuli)),
                    ('DiscStars'                    , (floattype, Nannuli)),
                    ('SpinStars'                    , (floattype, 3)),
                    ('SpinGas'                      , (floattype, 3)),
                    ('SpinClassicalBulge'           , (floattype, 3)),
                    ('StarsInSitu'                  , floattype),
                    ('StarsInstability'             , floattype),
                    ('StarsMergeBurst'              , floattype),
                    ('DiscHI'                       , (floattype, Nannuli)),
                    ('DiscH2'                       , (floattype, Nannuli)),
                    ('DiscSFR'                      , (floattype, Nannuli)), 
                    ('MetalsColdGas'                , floattype),
                    ('MetalsStellarMass'            , floattype),
                    ('ClassicalMetalsBulgeMass'     , floattype),
                    ('SecularMetalsBulgeMass'       , floattype),
                    ('MetalsHotGas'                 , floattype),
                    ('MetalsEjectedMass'            , floattype),
                    ('MetalsIntraClusterStars'      , floattype),
                    ('DiscGasMetals'                , (floattype, Nannuli)),
                    ('DiscStarsMetals'              , (floattype, Nannuli)),
                    ('SfrFromH2'                    , floattype),
                    ('SfrInstab'                    , floattype),
                    ('SfrMergeBurst'                , floattype),
                    ('SfrDiskZ'                     , floattype),
                    ('SfrBulgeZ'                    , floattype),
                    ('DiskScaleRadius'              , floattype),
                    ('CoolScaleRadius'              , floattype), 
                    ('StellarDiscScaleRadius'       , floattype),
                    ('Cooling'                      , floattype),
                    ('Heating'                      , floattype),
                    ('LastMajorMerger'              , floattype),
                    ('LastMinorMerger'              , floattype),
                    ('OutflowRate'                  , floattype),
                    ('infallMvir'                   , floattype),
                    ('infallVvir'                   , floattype),
                    ('infallVmax'                   , floattype)
                    ]
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    return Galdesc



def darksage_out_single(fname, fields=[], Nannuli=30):
    # Read a single Dark Sage output file, returning all the galaxy data
    # fname is the full name for the file to read, including its path
    # fields is the list of fields you want to read in.  If empty, will read all fields.
    
    Galdesc = galdtype_darksage(Nannuli)
    if len(fields)==0: fields=list(Galdesc.names)
    
    fin = open(fname, 'rb')  # Open the file
    Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
    NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
    GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
    G = np.fromfile(fin, Galdesc, NtotGals) # Read all the galaxy data
    G = G[fields]
    return G 



def darksage_snap(fpre, filelist, fields=[], Nannuli=30):
    # Read full Dark Sage snapshot, going through each file and compiling into 1 array
    # fpre is the name of the file up until the _ before the file number
    # filelist contains all the file numbers you want to read in
    
    Galdesc = galdtype_darksage()
    Glist = []
    Ngal = np.array([],dtype=np.int32)
    G = darksage_out_single(fpre+'_'+str(filelist[0]), fields, Nannuli)
    
    for i in filelist[1:]:
        G1 = darksage_out_single(fpre+'_'+str(i), fields)
        G = np.append(G, G1)
    return G



def massfunction(mass, Lbox, range=[8,12.5], c='k', lw=2, ls='-', label='', ax=None):
    masslog = np.log10(mass[(mass>0)*np.isfinite(mass)])
    N, edges = np.histogram(masslog, bins=np.arange(range[0],range[1]+0.1,0.1))
    binwidth = edges[1]-edges[0]
    x = edges[:-1] + binwidth/2
    y = N/(binwidth*Lbox**3)
    
    if ax is None: ax = plt.gca()
    
    if len(label)>0:
        ax.plot(x, y, c+ls, linewidth=lw, label=label)
    else:
        ax.plot(x, y, c+ls, linewidth=lw)


def schechter(phistar, Mstar, alpha, Mlog=False, range=[7,12], Npoints=2000, logM=None):
    if Mlog: Mstar = 10**Mstar
    if logM is None: logM = np.linspace(range[0],range[1],Npoints)
    M = 10**logM
    Phi = np.log(10.) * (phistar) * (M/Mstar)**(alpha+1) * np.exp(-M/Mstar)
    return Phi, logM


def stellar_massfunction_obsdata(h=0.678, ax=None):
    B = np.array([
                  [7.05, 1.3531e-01, 6.0741e-02],
                  [7.15, 1.3474e-01, 6.0109e-02],
                  [7.25, 2.0971e-01, 7.7965e-02],
                  [7.35, 1.7161e-01, 3.1841e-02],
                  [7.45, 2.1648e-01, 5.7832e-02],
                  [7.55, 2.1645e-01, 3.9988e-02],
                  [7.65, 2.0837e-01, 4.8713e-02],
                  [7.75, 2.0402e-01, 7.0061e-02],
                  [7.85, 1.5536e-01, 3.9182e-02],
                  [7.95, 1.5232e-01, 2.6824e-02],
                  [8.05, 1.5067e-01, 4.8824e-02],
                  [8.15, 1.3032e-01, 2.1892e-02],
                  [8.25, 1.2545e-01, 3.5526e-02],
                  [8.35, 9.8472e-02, 2.7181e-02],
                  [8.45, 8.7194e-02, 2.8345e-02],
                  [8.55, 7.0758e-02, 2.0808e-02],
                  [8.65, 5.8190e-02, 1.3359e-02],
                  [8.75, 5.6057e-02, 1.3512e-02],
                  [8.85, 5.1380e-02, 1.2815e-02],
                  [8.95, 4.4206e-02, 9.6866e-03],
                  [9.05, 4.1149e-02, 1.0169e-02],
                  [9.15, 3.4959e-02, 6.7898e-03],
                  [9.25, 3.3111e-02, 8.3704e-03],
                  [9.35, 3.0138e-02, 4.7741e-03],
                  [9.45, 2.6692e-02, 5.5029e-03],
                  [9.55, 2.4656e-02, 4.4359e-03],
                  [9.65, 2.2885e-02, 3.7915e-03],
                  [9.75, 2.1849e-02, 3.9812e-03],
                  [9.85, 2.0383e-02, 3.2930e-03],
                  [9.95, 1.9929e-02, 2.9370e-03],
                  [10.05, 1.8865e-02, 2.4624e-03],
                  [10.15, 1.8136e-02, 2.5208e-03],
                  [10.25, 1.7657e-02, 2.4217e-03],
                  [10.35, 1.6616e-02, 2.2784e-03],
                  [10.45, 1.6114e-02, 2.1783e-03],
                  [10.55, 1.4366e-02, 1.8819e-03],
                  [10.65, 1.2588e-02, 1.8249e-03],
                  [10.75, 1.1372e-02, 1.4436e-03],
                  [10.85, 9.1213e-03, 1.5816e-03],
                  [10.95, 6.1125e-03, 9.6735e-04],
                  [11.05, 4.3923e-03, 9.6254e-04],
                  [11.15, 2.5463e-03, 5.0038e-04],
                  [11.25, 1.4298e-03, 4.2816e-04],
                  [11.35, 6.4867e-04, 1.6439e-04],
                  [11.45, 2.8294e-04, 9.9799e-05],
                  [11.55, 1.0617e-04, 4.9085e-05],
                  [11.65, 3.2702e-05, 2.4546e-05],
                  [11.75, 1.2571e-05, 1.2571e-05],
                  [11.85, 8.4589e-06, 8.4589e-06],
                  [11.95, 7.4764e-06, 7.4764e-06]
                  ], dtype=np.float32)
    if ax is None: ax = plt.gca()
    ax.fill_between(B[:,0]+np.log10(0.7**2)-np.log10(h**2), (B[:,1]+B[:,2])*h**3, (B[:,1]-B[:,2])*h**3, facecolor='purple', alpha=0.2)
    ax.plot([1,1], [1,2], color='purple', linewidth=8, alpha=0.3, label=r'Baldry et al.~(2008)') # Just for the legend



def HIH2_massfunction_obsdata(h=0.678, HI=True, H2=True, K=True, OR=False, ax=None, Z=True, M=False, B=False):
    Zwaan = np.array([[6.933,   -0.333],
                      [7.057,   -0.490],
                      [7.209,   -0.698],
                      [7.365,   -0.667],
                      [7.528,   -0.823],
                      [7.647,   -0.958],
                      [7.809,   -0.917],
                      [7.971,   -0.948],
                      [8.112,   -0.927],
                      [8.263,   -0.917],
                      [8.404,   -1.062],
                      [8.566,   -1.177],
                      [8.707,   -1.177],
                      [8.853,   -1.312],
                      [9.010,   -1.344],
                      [9.161,   -1.448],
                      [9.302,   -1.604],
                      [9.448,   -1.792],
                      [9.599,   -2.021],
                      [9.740,   -2.406],
                      [9.897,   -2.615],
                      [10.053,  -3.031],
                      [10.178,  -3.677],
                      [10.335,  -4.448],
                      [10.492,  -5.083]])
        
    Martin_data = np.array([[6.302,    -0.504],
                              [6.500,    -0.666],
                              [6.703,    -0.726],
                              [6.904,    -0.871],
                              [7.106,    -1.135],
                              [7.306,    -1.047],
                              [7.504,    -1.237],
                              [7.703,    -1.245],
                              [7.902,    -1.254],
                              [8.106,    -1.414],
                              [8.306,    -1.399],
                              [8.504,    -1.476],
                              [8.705,    -1.591],
                              [8.906,    -1.630],
                              [9.104,    -1.695],
                              [9.309,    -1.790],
                              [9.506,    -1.981],
                              [9.707,    -2.141],
                              [9.905,    -2.317],
                              [10.108,    -2.578],
                              [10.306,    -3.042],
                              [10.509,    -3.780],
                              [10.703,    -4.534],
                              [10.907,    -5.437]])     
                    
    Martin_mid = Martin_data[:,1] + 3*np.log10(h/0.7)
    Martin_x = Martin_data[:,0] + 2*np.log10(0.7/h)
    Martin_high = np.array([-0.206, -0.418, -0.571, -0.725, -1.003, -0.944, -1.144, -1.189, -1.189, -1.358, -1.344, -1.417, -1.528, -1.586, -1.651, -1.753, -1.925, -2.095, -2.281, -2.537, -3.003, -3.729, -4.451, -5.222]) + 3*np.log10(h/0.7)
    Martin_low = np.array([-0.806, -0.910, -0.885, -1.019, -1.268, -1.173, -1.313, -1.314, -1.320, -1.459, -1.443, -1.530, -1.647, -1.669, -1.736, -1.838, -2.021, -2.191, -2.359, -2.621, -3.098, -3.824, -4.618, -5.663]) + 3*np.log10(h/0.7)
      
    Keres_high = np.array([-1.051, -1.821, -1.028, -1.341, -1.343, -1.614, -1.854, -2.791,  -3.54 , -5.021]) + 3*np.log10(h)
    Keres_mid = np.array([-1.271, -1.999, -1.244, -1.477, -1.464, -1.713, -1.929, -2.878,   -3.721, -5.22 ]) + 3*np.log10(h)
    Keres_low = np.array([-1.706, -2.302, -1.71 , -1.676, -1.638, -1.82 , -2.033, -2.977,   -4.097, -5.584]) + 3*np.log10(h)
    Keres_M = np.array([  6.953,   7.353,   7.759,   8.154,   8.553,   8.96 ,   9.365,  9.753,  10.155,  10.558]) - 2*np.log10(h)
      
    ObrRaw_high = np.array([-0.905, -1.122, -1.033, -1.1  , -1.242, -1.418, -1.707, -2.175, -2.984, -4.868]) + 3*np.log10(h)
    ObrRaw_mid = np.array([-1.116, -1.308, -1.252, -1.253, -1.373, -1.509, -1.806, -2.261,  -3.198, -5.067]) + 3*np.log10(h)
    ObrRaw_low = np.array([-1.563, -1.602, -1.73 , -1.448, -1.537, -1.621, -1.918, -2.369,  -3.556, -5.413]) + 3*np.log10(h)
    ObrRaw_M = np.array([ 7.301,  7.586,  7.862,  8.133,  8.41 ,  8.686,  8.966,  9.242,    9.514,  9.788]) - 2*np.log10(h)
      
    HI_x = Zwaan[:,0] - 2*np.log10(h)
    HI_y = 10**Zwaan[:,1] * h**3

    Boselli_const_XCO = np.array([[7.39189, -3.06989, -3.32527, -2.86828],
                                [7.78378, -2.45161, -2.54570, -2.37097],
                                [8.18919, -1.91398, -1.96774, -1.84677],
                                [8.62162, -2.12903, -2.20968, -2.03495],
                                [9.01351, -2.41129, -2.51882, -2.31720],
                                [9.41892, -2.62634, -2.80108, -2.53226],
                                [9.81081, -2.73387, -2.85484, -2.54570],
                                [10.2297, -3.64785, -5.97312, -3.36559]])

    Boselli_var_XCO = np.array([[7.59030, -3.19086, -3.58065, -2.98925],
                              [7.98113, -2.55914, -2.72043, -2.45161],
                              [8.37197, -2.22312, -2.30376, -2.14247],
                              [8.78976, -1.94086, -1.99462, -1.90054],
                              [9.18059, -1.98118, -2.06183, -1.90054],
                              [9.59838, -2.72043, -2.92204, -2.62634], 
                              [9.98922, -3.67473, -5.98656, -3.31183]])
    
    if ax is None: ax = plt.gca()         
    if HI and Z: ax.plot(HI_x, HI_y, '-', color='g', lw=8, alpha=0.4, label=r'Zwaan et al.~(2005)')
    if HI and M: ax.fill_between(Martin_x, 10**Martin_high, 10**Martin_low, color='c', alpha=0.4)
    if HI and M: ax.plot([0,1], [1,1], 'c-', lw=8, alpha=0.4, label=r'Martin et al.~(2010)')

    if H2 and K: ax.fill_between(Keres_M, 10**Keres_high, 10**Keres_low, color='teal', alpha=0.4)
    if H2 and K: ax.plot([0,1], [1,1], '-', color='teal', lw=8, alpha=0.4, label=r'Keres et al.~(2003)')

    if H2 and OR: ax.fill_between(ObrRaw_M, 10**ObrRaw_high, 10**ObrRaw_low, color='darkcyan', alpha=0.4)
    if H2 and OR: ax.plot([0,1], [1,1], '-', color='darkcyan', lw=8, alpha=0.4, label=r'Obreschkow \& Rawlings (2009)')
                        
    if H2 and B: ax.fill_between(Boselli_const_XCO[:,0]+2*np.log10(0.7/h), 10**Boselli_const_XCO[:,2]*(h/0.7)**3/0.4, 10**Boselli_const_XCO[:,3]*(h/0.7)**3/0.4, color='orange', alpha=0.4)
    if H2 and B: ax.plot([0,1], [1,1], '-', color='orange', lw=8, alpha=0.4, label=r'Boselli et al.~(2014), const.~$X_{\rm CO}$')
    if H2 and B: ax.fill_between(Boselli_var_XCO[:,0]+2*np.log10(0.7/h), 10**Boselli_var_XCO[:,2]*(h/0.7)**3/0.4, 10**Boselli_var_XCO[:,3]*(h/0.7)**3/0.4, color='violet', alpha=0.4)
    if H2 and B: ax.plot([0,1], [1,1], '-', color='violet', lw=8, alpha=0.4, label=r'Boselli et al.~(2014), var.~$X_{\rm CO}$')

    ax.set_xlabel(r'$\log_{10}(M_{\mathrm{H}\,\huge\textsc{i}}\ \mathrm{or}\ M_{\mathrm{H}_2}\ [\mathrm{M}_{\bigodot}])$')
    ax.set_ylabel(r'$\Phi\ [\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1}]$')
    ax.axis([8,11.5,1e-6,1e-1])



def savepng(filename, xsize=1024, ysize=None, fig=None, transparent=False):
    # Save a figure as a PNG with a normalised size / aspect ratio
    xpix = 2560
    ypix = 1440
    ss = 27
    
    if ysize==None: ysize = int(xsize*9./16)
    
    mydpi = np.sqrt(xpix**2 + ypix**2)/ss 
    xinplot = xsize*(9./7.)/mydpi
    yinplot = ysize*(9./7.)/mydpi
    if fig is None: fig = plt.gcf()
    fig.set_size_inches(xinplot,yinplot)
    fig.set_dpi(mydpi)
    
    filename = str(filename)
    if filename[-4:] != '.png':
        filename = filename+'.png'
    fig.savefig(filename, dpi=mydpi, bbox_inches='tight', transparent=transparent)


def Brown_HI_fractions(h):
    
    logM = np.array([[9.2209911, 9.6852989, 10.180009, 10.665453, 11.098589],
         [9.2085762, 9.6402225, 10.141238, 10.599669, 11.026575],
         [9.2121296, 9.6528578, 10.139588, 10.615245, 11.054683]]) + 2*np.log10(0.7/h) + np.log10(0.61/0.66)
         
    logHIfrac = np.array([[ 0.37694988,  0.0076254,  -0.45345795, -0.90604609, -1.39503932],
           [ 0.15731917, -0.16941574, -0.6199488,  -0.99943721, -1.30476058],
           [ 0.19498822, -0.27559358, -0.74410361, -1.12869251, -1.49363434]]) - np.log10(0.61/0.66)

    Ngal = np.array([[120, 318, 675, 1132, 727],
                      [3500, 4359, 3843, 2158, 268],
                      [2203, 3325, 2899, 1784, 356]])

    logM_all = np.log10((10**logM[0,:] * Ngal[0,:] + 10**logM[1,:] * Ngal[1,:] + 10**logM[2,:] * Ngal[2,:]) / (Ngal[0,:]+Ngal[1,:]+Ngal[2,:]))
    
    logHIfrac_all = np.log10((10**logHIfrac[0,:] * Ngal[0,:] + 10**logHIfrac[1,:] * Ngal[1,:] + 10**logHIfrac[2,:] * Ngal[2,:]) / (Ngal[0,:]+Ngal[1,:]+Ngal[2,:]))
    
    logHIfrac_all_err = np.array([np.log10(1.511+0.011)-np.log10(1.511),
                                  np.log10(0.643+0.011)-np.log10(0.643),
                                  np.log10(0.232+0.005)-np.log10(0.232),
                                  np.log10(0.096+0.002)-np.log10(0.096),
                                  np.log10(0.039+0.001)-np.log10(0.039)]) # taken straight from Table 1 of Brown+15

    return logM_all, logHIfrac_all, logHIfrac_all_err

def hist_Nmin(x, bins, Nmin, hard_bins=np.array([])):
    Nhist, bins = np.histogram(x, bins)
    while len(Nhist[Nhist<Nmin])>0:
        ii = np.where(Nhist<Nmin)[0][0]
        if (ii==0 or (ii!=len(Nhist)-1 and Nhist[ii+1]<Nhist[ii-1])) and np.all(~((bins[ii+1] <= 1.01*hard_bins) * (bins[ii+1] >= 0.99*hard_bins))):
            bins = np.delete(bins,ii+1)
        elif np.all(~((bins[ii] <= 1.01*hard_bins) * (bins[ii] >= 0.99*hard_bins))):
            bins = np.delete(bins,ii)
        else:
            print 'hard_bins prevented routines.hist_Nmin() from enforcing Nmin.  Try using wider input bins.'
            Nhist, bins = np.histogram(x, bins)
            break
        Nhist, bins = np.histogram(x, bins)
    if bins[0]<np.min(x): bins[0] = np.min(x)
    if bins[-1]>np.max(x): bins[-1] = np.max(x)
    return Nhist, bins


def meanbins(x, y, xmeans, tol=0.02, itmax=100):
    # Find bins in some dataset (x,y) which will have mean values of x matching xmeans
    fnan = np.isfinite(x)*np.isfinite(y)
    x, y = x[fnan], y[fnan]
    N = len(xmeans)
    bins = np.zeros(N+1)
    mean_x, mean_y = np.zeros(N), np.zeros(N)
    bins[0] = np.min(x)
    bins[-1] = np.max(x)
    for i in range(1,N+1):
        xleft = np.sort(x[x>=bins[i-1]])
        cumav = np.cumsum(xleft) / np.arange(1,len(xleft)+1)
        diff = abs(xmeans[i-1]-cumav)
        arg = np.where(diff==np.min(diff))[0][0]
        if i!=N and arg<len(diff)-1: arg += 1
        bins[i] = xleft[arg]
        m = cumav[arg]
        mean_x[i-1] = m
        f = (xleft<bins[i]) if i!=N else (xleft<=bins[i])
        mean_y[i-1] = np.mean(y[np.in1d(x,xleft[f])])
    return bins, mean_x, mean_y


def percentiles(x, y, low=0.16, med=0.5, high=0.84, bins=20, addMean=False, xrange=None, yrange=None, Nmin=10, weights=None, hard_bins=np.array([]), outBins=False):
    # Given some values to go on x and y axes, bin them along x and return the percentile ranges
    f = np.isfinite(x)*np.isfinite(y)
    if xrange is not None: f = (x>=xrange[0])*(x<=xrange[1])*f
    if yrange is not None: f = (y>=yrange[0])*(y<=yrange[1])*f
    x, y = x[f], y[f]
    if type(bins)==int:
        if len(x)/bins < Nmin: bins = len(x)/Nmin
        indices = np.array(np.linspace(0,len(x)-1,bins+1), dtype=int)
        bins = np.sort(x)[indices]
    elif Nmin>0: # Ensure a minimum number of data in each bin
        Nhist, bins = hist_Nmin(x, bins, Nmin, hard_bins)
    Nbins = len(bins)-1
    y_low, y_med, y_high = np.zeros(Nbins), np.zeros(Nbins), np.zeros(Nbins)
    x_av, N = np.zeros(Nbins), np.zeros(Nbins)
    if addMean: y_mean = np.zeros(Nbins)
    for i in range(Nbins):
        f = (x>=bins[i])*(x<bins[i+1])
        if len(f[f])>2:
            if weights is None:
                [y_low[i], y_med[i], y_high[i]] = np.percentile(y[f], [100*low, 100*med, 100*high], interpolation='linear')
            else:
                [y_low[i], y_med[i], y_high[i]] = weighted_percentile(y[f], [low, med, high], weights[f])
            x_av[i] = np.mean(x[f])
            N[i] = len(x[f])
            if addMean: y_mean[i] = np.mean(y[f])
    fN = (N>0) if Nmin>0 else np.array([True]*Nbins)
    if not addMean and not outBins:
        return x_av[fN], y_high[fN], y_med[fN], y_low[fN]
    elif not addMean and outBins:
        return x_av[fN], y_high[fN], y_med[fN], y_low[fN], bins
    elif addMean and  not outBins:
        return x_av[fN], y_high[fN], y_med[fN], y_low[fN], y_mean[fN]
    else:
        return x_av[fN], y_high[fN], y_med[fN], y_low[fN], y_mean[fN], bins



def Tremonti04(h):
    x_obs = np.array([8.52, 8.57, 8.67, 8.76, 8.86, 8.96, 9.06, 9.16, 9.26, 9.36, 9.46, 9.57, 9.66, 9.76, 9.86, 9.96, 10.06, 10.16, 10.26, 10.36, 10.46, 10.56, 10.66, 10.76, 10.86, 10.95, 11.05, 11.15, 11.25, 11.30])
    y_low = np.array([8.25, 8.25, 8.28, 8.32, 8.37, 8.46, 8.56, 8.59, 8.60, 8.63, 8.66, 8.69, 8.72, 8.76, 8.80, 8.83, 8.85, 8.88, 8.92, 8.94, 8.96, 8.98, 9.00, 9.01, 9.02, 9.03, 9.03, 9.04, 9.03, 9.03])
    y_high= np.array([8.64, 8.64, 8.65, 8.70, 8.73, 8.75, 8.82, 8.82, 8.86, 8.88, 8.92, 8.94, 8.96, 8.99, 9.01, 9.05, 9.06, 9.09, 9.10, 9.11, 9.12, 9.14, 9.15, 9.15, 9.16, 9.17, 9.17, 9.18, 9.18, 9.18])
    x_obs += np.log10(1.5/1.8) + 2*np.log10(0.7/h) # Accounts for difference in Kroupa & Chabrier IMFs and the difference in h
    return x_obs, y_low, y_high




def BH_bulge_obs(h=0.678, ax=None):
    M_BH_obs = (0.7/h)**2*1e8*np.array([39, 11, 0.45, 25, 24, 0.044, 1.4, 0.73, 9.0, 58, 0.10, 8.3, 0.39, 0.42, 0.084, 0.66, 0.73, 15, 4.7, 0.083, 0.14, 0.15, 0.4, 0.12, 1.7, 0.024, 8.8, 0.14, 2.0, 0.073, 0.77, 4.0, 0.17, 0.34, 2.4, 0.058, 3.1, 1.3, 2.0, 97, 8.1, 1.8, 0.65, 0.39, 5.0, 3.3, 4.5, 0.075, 0.68, 1.2, 0.13, 4.7, 0.59, 6.4, 0.79, 3.9, 47, 1.8, 0.06, 0.016, 210, 0.014, 7.4, 1.6, 6.8, 2.6, 11, 37, 5.9, 0.31, 0.10, 3.7, 0.55, 13, 0.11])
    M_BH_hi = (0.7/h)**2*1e8*np.array([4, 2, 0.17, 7, 10, 0.044, 0.9, 0.0, 0.9, 3.5, 0.10, 2.7, 0.26, 0.04, 0.003, 0.03, 0.69, 2, 0.6, 0.004, 0.02, 0.09, 0.04, 0.005, 0.2, 0.024, 10, 0.1, 0.5, 0.015, 0.04, 1.0, 0.01, 0.02, 0.3, 0.008, 1.4, 0.5, 1.1, 30, 2.0, 0.6, 0.07, 0.01, 1.0, 0.9, 2.3, 0.002, 0.13, 0.4, 0.08, 0.5, 0.03, 0.4, 0.38, 0.4, 10, 0.2, 0.014, 0.004, 160, 0.014, 4.7, 0.3, 0.7, 0.4, 1, 18, 2.0, 0.004, 0.001, 2.6, 0.26, 5, 0.005])
    M_BH_lo = (0.7/h)**2*1e8*np.array([5, 2, 0.10, 7, 10, 0.022, 0.3, 0.0, 0.8, 3.5, 0.05, 1.3, 0.09, 0.04, 0.003, 0.03, 0.35, 2, 0.6, 0.004, 0.13, 0.1, 0.05, 0.005, 0.2, 0.012, 2.7, 0.06, 0.5, 0.015, 0.06, 1.0, 0.02, 0.02, 0.3, 0.008, 0.6, 0.5, 0.6, 26, 1.9, 0.3, 0.07, 0.01, 1.0, 2.5, 1.5, 0.002, 0.13, 0.9, 0.08, 0.5, 0.09, 0.4, 0.33, 0.4, 10, 0.1, 0.014, 0.004, 160, 0.007, 3.0, 0.4, 0.7, 1.5, 1, 11, 2.0, 0.004, 0.001, 1.5, 0.19, 4, 0.005])
    M_sph_obs = (0.7/h)**2*1e10*np.array([69, 37, 1.4, 55, 27, 2.4, 0.46, 1.0, 19, 23, 0.61, 4.6, 11, 1.9, 4.5, 1.4, 0.66, 4.7, 26, 2.0, 0.39, 0.35, 0.30, 3.5, 6.7, 0.88, 1.9, 0.93, 1.24, 0.86, 2.0, 5.4, 1.2, 4.9, 2.0, 0.66, 5.1, 2.6, 3.2, 100, 1.4, 0.88, 1.3, 0.56, 29, 6.1, 0.65, 3.3, 2.0, 6.9, 1.4, 7.7, 0.9, 3.9, 1.8, 8.4, 27, 6.0, 0.43, 1.0, 122, 0.30, 29, 11, 20, 2.8, 24, 78, 96, 3.6, 2.6, 55, 1.4, 64, 1.2])
    M_sph_hi = (0.7/h)**2*1e10*np.array([59, 32, 2.0, 80, 23, 3.5, 0.68, 1.5, 16, 19, 0.89, 6.6, 9, 2.7, 6.6, 2.1, 0.91, 6.9, 22, 2.9, 0.57, 0.52, 0.45, 5.1, 5.7, 1.28, 2.7, 1.37, 1.8, 1.26, 1.7, 4.7, 1.7, 7.1, 2.9, 0.97, 7.4, 3.8, 2.7, 86, 2.1, 1.30, 1.9, 0.82, 25, 5.2, 0.96, 4.9, 3.0, 5.9, 1.2, 6.6, 1.3, 5.7, 2.7, 7.2, 23, 5.2, 0.64, 1.5, 105, 0.45, 25, 10, 17, 2.4, 20, 67, 83, 5.2, 3.8, 48, 2.0, 55, 1.8])
    M_sph_lo = (0.7/h)**2*1e10*np.array([32, 17, 0.8, 33, 12, 1.4, 0.28, 0.6, 9, 10, 0.39, 2.7, 5, 1.1, 2.7, 0.8, 0.40, 2.8, 12, 1.2, 0.23, 0.21, 0.18, 2.1, 3.1, 0.52, 1.1, 0.56, 0.7, 0.51, 0.9, 2.5, 0.7, 2.9, 1.2, 0.40, 3.0, 1.5, 1.5, 46, 0.9, 0.53, 0.8, 0.34, 13, 2.8, 0.39, 2.0, 1.2, 3.2, 0.6, 3.6, 0.5, 2.3, 1.1, 3.9, 12, 2.8, 0.26, 0.6, 57, 0.18, 13, 5, 9, 1.3, 11, 36, 44, 2.1, 1.5, 26, 0.8, 30, 0.7])
    core = np.array([1,1,0,1,1,0,0,0,1,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,1,0,1,0])
    yerr2, yerr1 = np.log10(divide(M_BH_obs+M_BH_hi, M_BH_obs)), -np.log10(divide(M_BH_obs-M_BH_lo, M_BH_obs))
    xerr2, xerr1 = np.log10(divide(M_sph_obs+M_sph_hi, M_sph_obs)), -np.log10(divide(M_sph_obs-M_sph_lo, M_sph_obs))
    if ax is None: ax = plt.gca()
    ax.errorbar(np.log10(M_sph_obs[core==1]), np.log10(M_BH_obs[core==1]), yerr=[yerr1[core==1],yerr2[core==1]], xerr=[xerr1[core==1],xerr2[core==1]], color='c', alpha=0.5, label=r'Scott et al.~(2013) S\`{e}rsic', ls='none', lw=2, ms=0)
    ax.errorbar(np.log10(M_sph_obs[core==0]), np.log10(M_BH_obs[core==0]), yerr=[yerr1[core==0],yerr2[core==0]], xerr=[xerr1[core==0],xerr2[core==0]], color='purple', alpha=0.3, label=r'Core-S\`{e}rsic', ls='none', lw=2, ms=0)




def Leroygals(HI=False, H2=False, HighVvir=True, LowVvir=False, ax=None, SFR=False, h=0.678, c='k', alpha=0.5, lw=2):
    # Plot galaxy surface density profiles for select galaxies.  Stars done by default
    N628 = np.array([[0.2, 0.5, 0.9, 1.2, 1.6, 1.9, 2.3, 2.7, 3.0, 3.4, 3.7, 4.1, 4.4, 4.8, 5.1, 5.5, 5.8, 6.2, 6.5, 6.9, 7.3, 7.6, 8.0, 8.3, 8.7, 9.0, 9.4, 9.7, 10.1, 10.4, 10.8, 11.1, 11.5, 11.9, 12.2],
     [1.6, 2.1, 2.6, 3.1, 3.7, 4.6, 5.3, 5.8, 6.1, 6.5, 7.3, 7.9, 8.1, 7.9, 8.2, 8.5, 8.6, 8.6, 8.8, 8.8, 8.6, 8.2, 7.6, 7.1, 6.7, 6.5, 6.0, 5.2, 4.5, 4.1, 3.9, 3.9, 4.0, 4.3, 4.6],
     [0.3, 0.3, 0.4, 0.4, 0.3, 0.3, 0.4, 0.5, 0.5, 0.5, 0.7, 0.8, 0.8, 0.9, 1.0, 1.0, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 0.6, 0.6, 0.5, 0.4, 0.5, 0.4, 0.4, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5],
     [22.7, 20.2, 16.1, 12.7, 11.4, 11.1, 11.1, 10.6, 8.9, 7.2, 6.2, 5.9, 5.4, 4.3, 3.1, 2.1, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [1.2, 1.3, 1.2, 0.8, 1.1, 1.2, 1.7, 1.9, 1.5, 1.2, 1.5, 1.7, 1.5, 1.1, 0.8, 0.7, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [1209.4, 557.8, 313.6, 231.9, 194.3, 163.5, 143.9, 123.5, 107.5, 151.0, 81.6, 68.0, 61.6, 48.3, 41.8, 37.0, 33.2, 37.0, 52.9, 19.5, 18.9, 18.7, 12.9, 17.6, 17.0, 10.8, 8.0, 7.5, 5.0, 4.1, 3.6, 3.9, 4.4, 9.5, 5.8],
     [18.3, 4.8, 1.0, 0.5, 0.5, 0.7, 0.8, 0.5, 0.4, 10.5, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.4, 2.3, 6.1, 0.1, 0.1, 0.7, 0.1, 1.3, 1.6, 0.4, 0.1, 0.2, 0.1, 0.0, 0.0, 0.1, 0.2, 0.9, 0.2],
     [105.1, 92.3, 76.7, 65.5, 62.2, 72.4, 90.2, 90.7, 71.9, 57.9, 55.8, 59.6, 59.9, 48.8, 37.4, 33.5, 30.2, 23.5, 17.4, 13.6, 11.6, 9.8, 7.5, 5.4, 4.1, 3.2, 2.5, 1.8, 1.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
     [14, 9.9, 5.1, 4.2, 3.7, 12.2, 23.5, 21.3, 11.9, 8.5, 11.3, 14.1, 15.2, 11, 6.6, 8.7, 10, 6.5, 3.1, 1.9, 2.3, 2.2, 1.5, 0.9, 0.6, 0.4, 0.3, 0.3, 0.2, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0]])
                     
                     
    N3184 = np.array([[0.3, 0.8, 1.3, 1.9, 2.4, 3.0, 3.5, 4.0, 4.6, 5.1, 5.7, 6.2, 6.7, 7.3, 7.8, 8.3, 8.9, 9.4, 10.0, 10.5, 11.0, 11.6, 12.1, 12.6, 13.2, 13.7, 14.3],
    [3.7, 3.2, 3.3, 3.8, 4.7, 5.5, 5.7, 5.9, 6.5, 7.3, 7.5, 7.8, 8.1, 8.0, 7.3, 7.0, 7.0, 6.7, 6.1, 5.4, 5.0, 4.6, 4.0, 3.3, 2.9, 2.8, 2.7],
    [0.5, 0.3, 0.3, 0.3, 0.5, 0.4, 0.4, 0.5, 0.4, 0.3, 0.5, 0.6, 0.6, 0.5, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2],
    [44.2, 20.8, 14.5, 11.9, 12.6, 12.6, 11.0, 9.6, 7.4, 6.2, 5.5, 4.3, 2.7, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [9.6, 3.3, 2.0, 1.6, 2.0, 2.1, 2.1, 0.9, 1.0, 0.6, 0.9, 0.8, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [701.8, 270.5, 200.3, 146.5, 121.1, 113.0, 100.2, 94.2, 83.6, 74.6, 96.2, 61.2, 46.3, 34.6, 27.5, 22.4, 19.3, 14.9, 12.9, 9.5, 7.9, 6.6, 5.1, 4.7, 3.7, 3.0, 3.8],
    [22.8, 1.1, 0.6, 1.0, 0.5, 0.4, 0.4, 0.5, 0.3, 0.3, 6.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.1],
    [282.1, 85.4, 47.1, 47.1, 48.6, 50.8, 51.3, 50.0, 42.5, 38.4, 36.5, 37.6, 29.9, 20.4, 13.9, 9.9, 8.4, 6.9, 4.3, 3.0, 2.4, 1.7, 1.1, 0.0, 0.0, 0.0, 0.0],
    [78.8, 23.5, 6.5, 4.8, 6.9, 7.7, 8.8, 10.4, 6.5, 4.7, 4.5, 7.2, 5.0, 2.3, 1.4, 1.1, 1.1, 0.8, 0.4, 0.6, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0]])

    N3521 = np.array([[0.3, 0.8, 1.3, 1.8, 2.3, 2.9, 3.4, 3.9, 4.4, 4.9, 5.4, 6.0, 6.5, 7.0, 7.5, 8.0, 8.6, 9.1, 9.6, 10.1, 10.6, 11.2, 11.7, 12.2, 12.7, 13.2, 13.7, 14.3, 14.8, 15.3],
    [4.5, 4.8, 5.7, 7.1, 8.3, 8.7, 8.7, 8.9, 9.9, 10.6, 10.5, 10.2, 10.2, 9.4, 8.6, 8.4, 8.5, 8.6, 8.5, 8.2, 8.1, 8.1, 8.2, 8.3, 8.2, 8.0, 7.9, 7.7, 7.1, 6.5],
    [0.2, 0.3, 0.6, 0.7, 0.8, 0.8, 0.6, 0.8, 0.9, 0.6, 0.5, 0.4, 0.5, 0.3, 0.4, 0.5, 0.5, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9, 0.9, 0.9, 0.9, 1.0, 0.9, 0.8, 0.6],
    [25.7, 35.4, 43.4, 44.6, 41.5, 36.8, 30.6, 24.6, 22.2, 21.0, 17.5, 12.5, 8.4, 5.2, 3.2, 2.1, 1.7, 1.6, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [4.7, 5.5, 3.0, 1.5, 2.1, 2.7, 2.6, 2.3, 2.0, 2.1, 2.4, 2.5, 2.0, 1.4, 0.9, 0.6, 0.5, 0.5, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [4545.9, 1442.2, 929.8, 589.1, 462.9, 381.3, 322.9, 250.8, 212.3, 192.2, 169.7, 134.7, 106.7, 82.4, 66.1, 55.4, 47.7, 41.7, 35.8, 30.5, 27.1, 25.4, 22.9, 20.0, 17.4, 15.7, 14.4, 13.3, 12.0, 10.9],
    [287.3, 23.0, 8.9, 5.4, 2.9, 2.0, 1.9, 1.8, 1.4, 1.3, 1.1, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
    [106, 152.7, 183.4, 186.8, 176.9, 160.3, 138.2, 109, 89.2, 80.1, 70.3, 53.5, 36.1, 23.5, 16.2, 12.9, 11.6, 10.7, 9.4, 7.5, 6.1, 5.3, 4.6, 3.9, 3.1, 2.5, 2.1, 1.8, 1.4, 1.1],
    [27.9, 22.8, 11.5, 8.6, 8.7, 10.9, 16.3, 15.1, 11.5, 9.3, 9.2, 9.2, 7.6, 5.2, 3.6, 3.0, 3.0, 2.5, 2.0, 1.3, 1.3, 1.3, 1.2, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2]])

    N5194 = np.array([[0.2, 0.6, 1.0, 1.4, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1, 4.5, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.3, 8.7, 9.1, 9.5, 9.9, 10.3, 10.7],
    [4.5, 5.5, 6.1, 6.1, 6.7, 7.9, 8.5, 7.5, 6.2, 6.1, 7.2, 9.1, 11.2, 12.8, 12.7, 11.1, 9.4, 8.4, 7.8, 7.8, 7.8, 7.8, 7.8, 7.3, 6.4, 5.8, 5.1, 4.5],
    [0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.9, 0.9 ,0.7, 0.5, 0.7, 0.9, 0.9, 0.8, 0.8, 0.8, 0.9, 0.9, 0.9, 1.0, 1.1, 1.2, 1.2, 1.1, 1.0, 1.0, 0.9, 0.9],
    [197.4, 207.7, 181.6, 134.5, 106.8, 94.8, 72.6, 40.9, 19.1, 14.6, 22.6, 33.6, 35.9, 28.0, 14.9, 4.3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [35.0, 33.5, 30.5, 41.2, 33.8, 19.2, 15.7 ,12.2, 5.6, 3.3, 7.5, 10.3, 9.6, 7.0, 4.7, 8.2, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [4912.2, 2352.5, 1251.1, 703.1, 471.8, 417.5, 394.9, 334.1, 286.9, 253.2, 236.7, 224.2, 227.8, 206.0, 176.5, 148.9, 106.3, 77.5, 64.3, 50.2, 45.0, 46.4, 53.2, 69.2, 71.2, 86.5, 210.3, 102.7],
    [111.6, 15.3, 6.5, 3.8, 1.3, 1.6, 1.9, 1.5, 1.3, 2.0, 1.7, 1.2, 1.4, 1.1, 0.9, 0.8, 0.4, 0.3, 0.6, 0.3, 0.3, 0.3, 0.6, 1.8, 1.3, 1.7, 11.3, 3.1],
    [1164, 1026.9, 772.8, 463.7, 295.4, 296.7, 304.1, 234.9, 154.7, 113, 113.5, 148.7, 188, 199.4, 180.2, 141.1, 95.9, 58.7, 36.7, 24.9, 19.5, 19.5, 20, 19, 14.8, 10.8, 8.3, 6.6],
    [96.6, 138.1, 134.6, 84.7, 45.3, 49.5, 49.9, 34.2, 19.4, 8.7, 19.8, 39.9, 48.7, 41.2, 31.7, 27.6, 19.2, 9.9, 5.7, 4.3, 3.9, 4.3, 5.5, 6.8, 5.1, 2.8, 1.8, 1.3]])

    if ax is None: ax = plt.gca()

    if HighVvir and not SFR:
        if HI and H2:
            ax.errorbar(N628[0,:], N628[1,:]+N628[3,:], N628[2,:]+N628[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3184[0,:], N3184[1,:]+N3184[3,:], N3184[2,:]+N3184[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3521[0,:], N3521[1,:]+N3521[3,:], N3521[2,:]+N3521[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5194[0,:], N5194[1,:]+N5194[3,:], N5194[2,:]+N5194[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        elif HI and not H2:
            ax.errorbar(N628[0,:], N628[1,:], N628[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3184[0,:], N3184[1,:], N3184[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3521[0,:], N3521[1,:], N3521[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5194[0,:], N5194[1,:], N5194[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        elif H2 and not HI:
            ax.errorbar(N628[0,:], N628[3,:], N628[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3184[0,:], N3184[3,:], N3184[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3521[0,:], N3521[3,:], N3521[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5194[0,:], N5194[3,:], N5194[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        else:
            ax.errorbar(N628[0,:], N628[5,:]*0.61/0.66, N628[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3184[0,:], N3184[5,:]*0.61/0.66, N3184[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3521[0,:], N3521[5,:]*0.61/0.66, N3521[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5194[0,:], N5194[5,:]*0.61/0.66, N5194[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
    elif HighVvir and SFR:
        ax.errorbar(N628[0,:], N628[7,:]*1e-4*0.63/0.67, N628[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N3184[0,:], N3184[7,:]*1e-4*0.63/0.67, N3184[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N3521[0,:], N3521[7,:]*1e-4*0.63/0.67, N3521[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N5194[0,:], N5194[7,:]*1e-4*0.63/0.67, N5194[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)

    N3351 = np.array([[0.2, 0.7, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.7, 5.1, 5.6, 6.1, 6.6, 7.1, 7.6, 8.1, 8.6, 9.1, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5],
    [1.5, 1.0, 0.0, 0.0, 1.2, 2.0, 2.6, 2.7, 2.9, 3.1, 3.0, 2.9, 2.7, 2.7, 2.7, 2.8, 3.0, 3.3, 3.5, 3.6, 3.5, 3.1, 2.5, 2.1, 1.8, 1.4],
    [0.3, 0.2, 1.0, 1.0, 0.1, 0.2, 0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1],
    [161.4, 72.1, 15.7, 3.0, 2.1, 3.6, 4.8, 4.3, 3.1, 2.3, 1.6, 1.3, 1.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [29.9, 22.6, 6.3, 1.5, 0.6, 0.7, 0.8, 1.0, 0.7, 0.4, 0.4, 0.4, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [4525.1, 973.0, 418.7, 245.3, 189.8, 713.4, 162.3, 124.0, 93.3, 69.6, 55.5, 49.1, 44.5, 38.8, 32.9, 26.3, 20.4, 16.1, 12.6, 10.7, 9.2, 6.7, 5.9, 5.3, 3.7, 3.4],
    [98.2, 10.4, 2.7, 1.7, 1.2, 0.8, 0.5, 0.4, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.0, 0.1, 0.2, 0.1, 0.2],
    [2559.3, 1052, 194.1, 44.5, 30.6, 33.1, 37.2, 32.3, 23.1, 17.1, 13.6, 11.2, 9.7, 9.0, 8.6, 6.9, 5.2, 4.3, 3.5, 3, 2.6, 1.8, 1.1, 0.0, 0.0, 0.0],
    [473.1, 337, 68.7, 7.9, 1.8, 2, 2.7, 2.8, 1.7, 1.3, 1.4, 1.2, 1.1, 0.8, 1.1, 1.1, 0.5, 0.4, 0.3, 0.3, 0.3, 0.3, 0.1, 0.0, 0.0, 0.0]])

    N3627 = np.array([[0.2, 0.7, 1.1, 1.6, 2.0, 2.5, 2.9, 3.4, 3.8, 4.3, 4.7, 5.2, 5.6, 6.1, 6.5, 7.0, 7.4, 7.9, 8.3, 8.8, 9.2, 9.7, 10.1, 10.6, 11.0, 11.5, 11.9, 12.4, 12.8, 13.3, 13.8, 14.2, 14.7, 15.1, 15.6, 16.0, 16.5],
    [3.0, 3.5, 4.1, 4.9, 5.7, 6.4, 7.3, 7.7, 7.2, 6.5, 6.0, 5.3, 4.5, 3.7, 3.2, 2.9, 2.3, 1.7, 1.3, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.6, 0.4, 0.6, 0.7, 1.0, 1.5, 1.6, 1.2, 1.1, 1.3, 1.2, 1.0, 0.9, 0.9, 0.9, 0.6, 0.4, 0.3, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [173.3, 109.6, 59.1, 37.3, 38.7, 47.9, 49.6, 39.2, 22.2, 13.4, 10.5, 8.4, 7.2, 6.1, 5.4, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [30.1, 23.0, 11.8, 8.9, 16.3, 23.9, 24.3, 19.2, 12.2, 7.0, 4.4, 3.8, 3.7, 3.6, 3.3, 2.8, 1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [5428.1, 1553.9, 846.5, 592.0, 478.9, 289.1, 254.9, 275.1, 216.8, 172.9, 158.3, 132.8, 113.8, 94.8, 78.5, 69.2, 53.4, 41.3, 35.0, 28.7, 24.3, 21.9, 22.2, 57.6, 59.3, 16.5, 10.8, 8.7, 7.5, 8.6, 9.0, 5.4, 4.6, 4.3, 4.3, 3.7, 3.0],
    [273.0, 14.4, 4.6, 3.5, 3.3, 4.6, 4.7, 2.5, 1.6, 1.0, 1.2, 0.9, 0.8, 0.6, 0.5, 0.6, 0.4, 0.2, 0.3, 0.1, 0.1, 0.1, 0.3, 6.3, 6.7, 0.5, 0.1, 0.1, 0.1, 0.5, 0.6, 0.1, 0.0, 0.1, 0.1, 0.1, 0.0],
    [251.1, 205.9, 167.2, 197.1, 317.6, 430.4, 431.9, 327.5, 192.2, 127.9, 101, 62.5, 37.6, 27.7, 25.6, 25.3, 19.2, 11.5, 7, 4.6, 3.3, 2.4, 1.7, 1.3, 1.1, 0.0, 0.0, 0.0, 0.0, 1.1, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [27.4, 24, 15.7, 56.5, 152.2, 209.4, 181.4, 128.4, 72.2, 43.2, 40.8, 19.6, 8.9, 8.1, 9.9, 12.8, 9.7, 4.6, 1.8, 0.9, 0.6, 0.4, 0.3, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.8, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    N5055 = np.array([[0.2, 0.7, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.7, 5.1, 5.6, 6.1, 6.6, 7.1, 7.6, 8.1, 8.6, 9.1, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.4, 14.9, 15.4, 15.9, 16.4, 16.9, 17.4, 17.9, 18.4, 18.9, 19.3, 19.8, 20.3, 20.8],
    [5.6, 5.8, 5.9, 5.9, 6.2, 6.6, 6.4, 6.2, 6.5, 7.2, 8.2, 8.7, 8.7, 8.5, 8.5, 8.6, 8.5, 7.9, 7.5, 7.4, 7.2, 7.3, 7.5, 7.3, 6.7, 6.4, 6.2, 5.7, 5.1, 4.5, 4.2, 3.9, 3.7, 3.2, 2.8, 2.6, 2.5, 2.4, 2.1, 1.7, 1.4, 1.3, 1.2],
    [0.7, 0.5, 0.4, 0.2, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.7, 0.5, 0.5, 0.6, 0.5, 0.5, 0.5, 0.4, 0.4, 0.5, 0.7, 0.8, 0.9, 0.7, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.5, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2],
    [142.7, 98.8, 62.2, 43.7, 36.6, 32.1, 25.3, 20.3, 19.1, 18.6, 18.8, 17.3, 13.4, 10.9, 10.3, 10.0, 8.7, 6.2, 4.3, 3.1, 2.1, 1.5, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [19.4, 16.1, 8.5, 3.9, 2.3, 2.5, 2.7, 2.3, 2.1, 2.1, 2.1, 1.8, 1.1, 1.0, 1.0, 1.2, 1.3, 1.0, 0.8, 0.6, 0.5, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [4742.4, 1627.7, 987.4, 758.0, 569.7, 417.9, 325.4, 264.3, 230.7, 194.9, 169.5, 150.6, 133.4, 109.1, 94.0, 84.5, 75.3, 62.6, 52.0, 44.5, 40.8, 36.8, 33.4, 60.4, 59.2, 24.8, 21.2, 18.7, 17.3, 18.0, 13.8, 13.2, 11.8, 10.7, 9.8, 9.2, 8.4, 8.4, 7.7, 7.4, 7.1, 13.2, 12.4],
    [251.1, 11.4, 4.6, 2.2, 1.7, 1.1, 0.9, 0.8, 0.7, 0.5, 0.4, 0.4, 0.4, 0.3, 0.2, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 5.4, 5.7, 0.1, 0.1, 0.1, 0.1, 0.5, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.2, 0.2, 1.0, 1.4],
    [249.8, 267.2, 238.7, 192.7, 162.5, 135.2, 104.4, 81.3, 70.8, 67.2, 70.7, 73.4, 63, 51.4, 46.6, 43.3, 38.5, 29.6, 22.7, 20.1, 16.8, 12.6, 10.2, 8.6, 6.7, 5.3, 4.4, 3.7, 3.2, 2.8, 2.5, 2.2, 2, 1.6, 1.2, 1, 1.1, 1.4, 1.2, 0.0, 0.0, 0.0, 0.0],
    [28.4, 20.9, 13.2, 8.7, 5.9, 5.8, 5.9, 5.4, 5.2, 5.8, 7.3, 9.3, 6.4, 4.6, 5.3, 6.4, 6.4, 3.6, 2.1, 3.6, 3.7, 2.6, 2.0, 1.4, 0.8, 0.5, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.1, 0.2, 0.4, 0.9, 0.8, 0.0, 0.0, 0.0, 0.0]])

    N6946 = np.array([[0.1, 0.4, 0.7, 1.0, 1.3, 1.6, 1.9, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.1, 4.4, 4.7, 5.0, 5.3, 5.6, 5.9, 6.1, 6.4, 6.7, 7.0, 7.3, 7.6, 7.9, 8.2, 8.4, 8.7, 9.0, 9.3, 9.6, 9.9, 10.2, 10.4, 10.7, 11.0, 11.3, 11.6],
    [6.1, 6.4, 6.4, 5.9, 5.5, 5.5, 5.8, 6.4, 6.9, 7.4, 7.8, 8.2, 8.7, 9.3, 9.5, 9.5, 9.6, 9.6, 9.5, 9.3, 9.3, 9.3, 9.1, 8.8, 8.4, 8.1, 8.0, 8.0, 7.9, 7.4, 6.9, 6.3, 5.8, 5.3, 4.8, 4.5, 4.1, 3.9, 3.7, 3.6, 3.6],
    [1.1, 1.1, 1.0, 0.7 ,0.6, 0.5, 0.4, 0.4 ,0.4, 0.4, 0.4, 0.5, 0.6, 0.8, 1.0, 1.1, 1.1, 1.1, 1.1, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.9, 1.0, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4],
    [548.6, 390.7, 214.2, 110.4, 64.2, 46.4, 39.9, 37.9, 36.9, 35.2, 32.4, 29.7, 28.1, 27.5, 25.9, 22.8, 19.2, 15.6, 12.2, 9.4, 7.5, 6.3, 5.1, 4.1, 3.2, 2.3, 1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [68.4, 81.2, 54.9, 31.1, 17.4, 10.6, 7.0, 5.4, 4.9, 4.3, 3.5, 3.3, 3.9, 4.9, 5.7, 5.6, 5.1, 4.4, 3.7, 2.7, 1.9, 1.5, 1.4, 1.4, 1.4, 1.5, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [5937.7, 1125.9, 708.0, 496.6, 405.6, 390.4, 356.3, 313.6, 287.1, 258.4, 236.2, 212.1, 200.8, 276.6, 187.3, 159.7, 146.3, 127.2, 110.4, 97.4, 184.5, 105.3, 99.9, 85.8, 71.1, 59.1, 49.0, 40.7, 43.0, 35.2, 35.4, 30.4, 22.6, 26.8, 36.6, 86.9, 20.4, 17.3, 17.2, 17.8, 12.9],
    [348.7, 7.7, 4.5, 2.6, 1.8, 4.7, 1.9, 0.9, 1.0, 1.0, 1.5, 0.9, 1.5, 12.8, 1.8, 1.3, 1.0, 0.7, 0.6, 0.4, 10.5, 2.0, 4.4, 3.8, 1.1, 0.5, 0.4, 0.3, 1.0, 0.9, 0.9, 1.1, 0.5, 0.8, 1.9, 6.3, 0.5, 0.4, 0.5, 0.7, 0.3],
    [2286, 1549.6, 771.9, 373.6, 227.6, 186.1, 184.5, 197.2, 205.9, 202.3, 198.2, 201.7, 210.5, 221.8, 230.3, 230.6, 210.8, 175, 138.6, 110, 92, 79.9, 73.2, 74, 80.4, 77.7, 64.9, 45.4, 31.3, 23.1, 18.8, 15.7, 12.9, 10.4, 7.9, 5.6, 3.8, 2.7, 2.3, 2, 1.5],
    [553.4, 487.4, 265, 113.4, 48.6, 27, 21.2, 21.3, 25.9, 31.2, 24.7, 33.8, 31.9, 42.6, 60.7, 65.8, 57.1, 42.8, 31.8, 25.6, 23.9, 22.7, 22.9, 27.8, 35.7, 36.6, 28.7, 17.9, 9.9, 6.3, 5.9, 5.5, 4.8, 3.7, 2.4, 1.5, 1.3, 1.3, 1.3, 1.4, 1.4]])

    if LowVvir and not SFR:
        if HI and H2:
            ax.errorbar(N3351[0,:], N3351[1,:]+N3351[3,:], N3351[2,:]+N3351[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3627[0,:], N3627[1,:]+N3627[3,:], N3627[2,:]+N3627[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5055[0,:], N5055[1,:]+N5055[3,:], N5055[2,:]+N5055[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N6946[0,:], N6946[1,:]+N6946[3,:], N6946[2,:]+N6946[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        elif HI and not H2:
            ax.errorbar(N3351[0,:], N3351[1,:], N3351[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3627[0,:], N3627[1,:], N3627[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5055[0,:], N5055[1,:], N5055[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N6946[0,:], N6946[1,:], N6946[2,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        elif H2 and not HI:
            ax.errorbar(N3351[0,:], N3351[3,:], N3351[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3627[0,:], N3627[3,:], N3627[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5055[0,:], N5055[3,:], N5055[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N6946[0,:], N6946[3,:], N6946[4,:], elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        else:
            ax.errorbar(N3351[0,:], N3351[5,:]*0.61/0.66, N3351[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N3627[0,:], N3627[5,:]*0.61/0.66, N3627[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N5055[0,:], N5055[5,:]*0.61/0.66, N5055[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
            ax.errorbar(N6946[0,:], N6946[5,:]*0.61/0.66, N6946[6,:]*0.61/0.66, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
    elif LowVvir and SFR:
        ax.errorbar(N3351[0,:], N3351[7,:]*1e-4*0.63/0.67, N3351[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N3627[0,:], N3627[7,:]*1e-4*0.63/0.67, N3627[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N5055[0,:], N5055[7,:]*1e-4*0.63/0.67, N5055[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
        ax.errorbar(N6946[0,:], N6946[7,:]*1e-4*0.63/0.67, N6946[8,:]*1e-4*0.63/0.67, elinewidth=1, ecolor=c, alpha=alpha, color=c, lw=lw)
