#matplotlib inline

import pylab as pl
import numpy as np
import scipy.optimize as opt

#inputs:

#fov=float(raw_input("S-band PB beam (arcmin): "))
#decmin=int(raw_input("Min declination (deg): "))
#decmax=int(raw_input("Max declination (deg): "))
#overhead=float(raw_input("(Multiplicative) overhead factor: "))
#fullsens=float(raw_input("Full survey sensitivity required (Jy/beam): "))
#effbw=float(raw_input("Effective RFI-free BW (Hz): "))
#nant=float(raw_input("N antennas (Hz): "))
#row_sep=float(raw_input("row separation (arcmin): ")

#fov=14.23	 # fwhm field of view (theta_pb) in arcminutes
fov=14.8	 # fwhm field of view (theta_pb) in arcminutes
decmin=-40
decmax=90
overhead=1.19
fullsens=69e-6
effbw=1.5e9
nspw=16		# number of spws, used in the data rate calculation
tdump=0.45
nant=26
row_sep=7.2
freq = 3.0e9
eta = 0.92    # correlator efficiency
bmax = 11.1e3    # longest baselines in meters for B config
#sefd = 370 # 3GHz from OSS center of the band
sefd = 350 # average across the band from Emmanuel's recent measurements

# calculate some things based on the above:

nbl=nant*(nant-1)/2  # number of baselines
resolution = 3600*np.degrees(2.997925e8/float(freq)/bmax)     # resolution in asec for freq in Hz, bmax in meters
sensitivity = np.sqrt(3) * fullsens # sensitivity per epoch (Jy)
dt = (sefd/(sensitivity*eta))**2/(2*nbl*effbw*2)    # required integration time in s (inverse of sensitivity eqn)
surveyspeed = 0.5665 * fov**2/dt        # ss in deg2/hr (or amin2/s), fov in amin, tint in s
scanrate = surveyspeed/row_sep
datarate = 45 * (nspw*64*4/16384.) / tdump    # does not include autocorrelations
min_tdump = 45*(nspw*64*4/16384.) / 25.    # maxdr is max data rate in MB/s; min tdump to stay within 25MB/s (not including autocorrelations)

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
print 'extra_time scaling at (negative) declination:',p1

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
print 'Example scaling as fcn of Dec:'
print '-40:', powerlawneg(-40, p1[0], p1[1])
print '-35:', powerlawneg(-35, p1[0], p1[1])
print '-30:', powerlawneg(-30, p1[0], p1[1])
print '-25:', powerlawneg(-25, p1[0], p1[1])
print '-20:', powerlawneg(-20, p1[0], p1[1])
#print '0:', powerlawneg(0, p1[0], p1[1])
#print '+70:', powerlawpos(70, (1/3.)*p1[0], p1[1])
#print '+80:', powerlawpos(80, (1/3.)*p1[0], p1[1])
#print '+90:', powerlawpos(90, (1/3.)*p1[0], p1[1])
print
print 'For uniform sensitivity at all Dec:'
print '\tTrue area %d. Effective area %d. Scaling in time %.2f' % (nomarea, effarea, scaling)

surveytime1 = 3 * overhead * effarea/surveyspeed
surveytime2 = 3 * overhead * nomarea/surveyspeed

print 'Total time including Tsys factor for low dec: %d' % surveytime1
print 'Total time for uniform survey speed: %d' % surveytime2
print 'Effective integration time (s): ', dt
print 'Survey speed: ',surveyspeed
print 'Scan rate in arcmin/s:', scanrate
print 'Min dump time for 25 MB/s rate limit:', min_tdump
print 'Fraction of beam slewed over per int:', (scanrate*min_tdump)/fov
