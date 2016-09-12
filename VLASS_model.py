import argparse
import numpy as np
try:
    import scipy.optimize as opt
except:
    print('scipy not available. Will not run fit to estimate survey time.')

parser = argparse.ArgumentParser(description='Parameters for VLASS design calculator')
parser.add_argument('--fov', type=float, help='S-band primary beam fwhm (arcmin)', default=14.786)
#default fov fixed to give same answer as TIP for survey speed and time request
parser.add_argument('--decmin', type=int, help='Min Declination (deg)', default=-40)
parser.add_argument('--decmax', type=int, help='Max Declination (deg)', default=90)
parser.add_argument('--overhead', type=float, help='Multiplicative overhead factor', default=1.19)
parser.add_argument('--failurerate', type=float, help='Multiplicative failure rate factor', default=1.03)
parser.add_argument('--fullsens', type=float, help='Full survey sensitivity required (Jy/beam)', default=69e-6)
parser.add_argument('--nepoch', type=int, help='Number of epochs in survey', default=3)
parser.add_argument('--effbw', type=float, help='Effective (RFI-free) bandwidth (Hz)', default=1.5e9)
parser.add_argument('--nant', type=int, help='Number of antennas', default=26)
parser.add_argument('--nchan', type=int, help='Number of channels', default=1024)  # default = 64 ch * 16 spw
parser.add_argument('--rowsep', type=float, help='Row separation (arcmin)', default=7.2)
parser.add_argument('--tdump', type=float, help='Dump time (s)', default=0.45)
parser.add_argument('--drlimit', type=float, help='Data rate limit (MB/s)', default=25)

args = parser.parse_args()
fov = args.fov
decmin = args.decmin
decmax = args.decmax
overhead = args.overhead
failurerate = args.failurerate
fullsens = args.fullsens
nepoch = args.nepoch
effbw = args.effbw
nant = args.nant
rowsep = args.rowsep
nchan = args.nchan
tdump = args.tdump
drlimit = args.drlimit

# other parameters (implicit assumption of S band and B config)
freq = 3.0e9
bmax = 11.1e3   # longest baselines in meters for B config
sefd = 350      # average across the band from Emmanuel's recent measurements
eta = 0.92      # correlator efficiency
nprod = 4       # number of correlation products (RR, LL, RL, LR)

# print inputs
print('***********************************')
print('*** VLASS Calculator Parameters ***')
print('***********************************')
print('\tField of view: {0}'.format(fov))
print('\tDeclination range: {0} to {1}'.format(decmin, decmax))
print('\tFull sensitivity: {0}'.format(fullsens))
print('\tNumber of epochs: {0}'.format(nepoch))
print('\tNumber of antennas: {0}'.format(nant))
print('\tEffective bandwidth: {0}'. format(effbw))
print('\tNumber of channels: {0}'.format(nchan))
print('\tOTFM row separation: {0}'.format(rowsep))
print('\tAssumed overhead factor: {0}'.format(overhead))
print('\tAssumed failure factor: {0}'.format(failurerate))
print('\tCorrelator dump time: {0}'.format(tdump))
print('\tData rate limit: {0}'.format(drlimit))
print('Implicit assumptions:')
print('\t Freq: {0}'.format(freq))
print('\t Longest baseline: {0}'.format(bmax))
print('\t SEFD: {0}'.format(sefd))
print('\t Correlator efficiency: {0}'.format(eta))
print('\t Correlator products: {0}'.format(nprod))
print('\n')

def fit_extra(decmin, decmax):

    # fit time scaling factor for declinations south of -20 at middle of S band
    # scale is square of nominal Tsys over by true Tsys (largely due to spillover)
    # note: VLSS used a max scaling factor of 2 in time at -40 dec.
    
    # for uniform sensitivity:
#extra_time = np.array([ [-40, 2.6], [-35, 1.75], [-30, 1.3], [-25, 1.2], [-20, 1.10]])  # new est (w/pyephem) and explicit +-1.5-hr transit

    # for not quite uniform sensitivity at low dec (take sqrt of above numbers):
#extra_time = np.array([ [-40, 1.61], [-35, 1.32], [-30, 1.14], [-25, 1.1], [-20, 1.05]])  # new est (w/pyephem) and explicit +-1.5-hr transit
    # above nos^0.6
#extra_time = np.array([ [-40, 1.77], [-35, 1.40], [-30, 1.17], [-25, 1.12], [-20, 1.06]])  # new est (w/pyephem) and explicit +-1.5-hr transit
    # above nos^0.7
    extra_time = np.array([ [-40, 1.95], [-35, 1.48], [-30, 1.2], [-25, 1.14], [-20, 1.07]])  # new est (w/pyephem) and explicit +-1.5-hr transit

    # fit powerlaw to this
    powerlawneg = lambda dec, a, alpha: 1 + a * (dec/float(decmin))**alpha   # seems to behave about right
    #powerlawpos = lambda dec, a, alpha: 1 + a * (dec/float(decmax))**alpha   # to scale neg to Condon times at pos Dec, then amp => 1/3 amp
    p1, p1cov = opt.curve_fit(powerlawneg, extra_time[:,0], extra_time[:,1], p0 = (1, 1))
    print('\tExtra_time scaling model at (negative) declination: {0}'.format(p1))

    # sum up sky area to estimate effective time on sky. ignore Tsys effect at positive Dec, since it is small.
    nominalarea1 = []
    nominalarea2 = []
    effectivearea1 = []
    effectivearea2 = []

#for dec in range(-40,90,1): # actually -40.5! so do calc for -40.5 and -39.5 and average
    for dec in range(decmin,decmax,1):
        area = np.abs(np.cos(np.radians(dec)) * np.degrees(2*np.pi))
        nominalarea1.append(area)
        if dec <= 0:
            effectivearea1.append(area*powerlawneg(dec, p1[0], p1[1]))
        elif dec > 0:
            effectivearea1.append(area)

    for dec in range(decmin+1,decmax,1):
        area = np.abs(np.cos(np.radians(dec)) * np.degrees(2*np.pi))
        nominalarea2.append(area)
        if dec <= 0:
            effectivearea2.append(area*powerlawneg(dec, p1[0], p1[1]))
        elif dec > 0:
            effectivearea2.append(area)

    nomarea=(np.sum(nominalarea1)+np.sum(nominalarea2))/2.
    effarea=(np.sum(effectivearea1)+np.sum(effectivearea2))/2.0
    scaling=effarea/nomarea

    print('\tExample scaling as fcn of Dec:')
    print('\t-40: {0}'.format(powerlawneg(-40, p1[0], p1[1])))
    print('\t-35: {0}'.format(powerlawneg(-35, p1[0], p1[1])))
    print('\t-30: {0}'.format(powerlawneg(-30, p1[0], p1[1])))
    print('\t-25: {0}'.format(powerlawneg(-25, p1[0], p1[1])))
    print('\t-20: {0}'.format(powerlawneg(-20, p1[0], p1[1])))
    print('\tFor uniform sensitivity at all Dec:')
    print '\tTrue area %d. Effective area %d. Scaling in time %.3f' % (nomarea, effarea, scaling)
    print('\n')

    return nomarea, effarea

# calculate some things based on the above:
nbl = nant*(nant-1)/2  # number of baselines
resolution = 3600*np.degrees(2.997925e8/float(freq)/bmax)     # resolution in asec for freq in Hz, bmax in meters
sensitivity = np.sqrt(nepoch) * fullsens # sensitivity per epoch (Jy)
dt = nepoch * (sefd/(sensitivity*eta))**2/(2*nbl*effbw*2)    # required integration time in s (inverse of sensitivity eqn)
surveyspeed = 0.5665 * fov**2/dt        # ss in deg2/hr (or amin2/s), fov in amin, tint in s
scanrate = surveyspeed/rowsep
datarate = 45 * (nchan*nprod/16384.) / tdump    # does not include autocorrelations
min_tdump = 45 * (nchan*nprod/16384.) / drlimit    # maxdr is max data rate in MB/s; min tdump to stay within 25MB/s (not including autocorrelations)

print('***************************************')
print('*** Output Survey Design Parameters ***')
print('***************************************')
print('\tEffective total integration time per point (s): {0}'.format(dt))
print('\tSurvey speed (deg2/hr): {0}'.format(surveyspeed))
print('\tScan rate (arcmin/s): {0}'.format(scanrate))
print('\tFraction of beam slewed over per int: {0}'.format((scanrate*min_tdump)/fov))
print('\tData rate (MB/s): {0}'.format(datarate))
print('\tMin dump time at {0} MB/s data rate limit (s): {1}'.format(drlimit, min_tdump))
print('\n')

try:
    nomarea, effarea = fit_extra(decmin, decmax)
    nomtime = nepoch * nomarea/surveyspeed
    nomtime_overhead = nepoch * overhead * nomarea/surveyspeed
    total_time = nepoch * failurerate * overhead * effarea/surveyspeed

    print('***********************************')
    print('***           Time              ***')
    print('***********************************')
    print('\tNominal survey time (no overhead or failures; non-uniform sensitivity): {0}'.format(nomtime))
    print('\tNominal survey time with overhead (no failures; non-uniform sensitivity): {0}'.format(nomtime_overhead))
    print('\tTotal time: {0}'.format(total_time))
except:
    pass
