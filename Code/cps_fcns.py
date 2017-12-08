import pyfits as pf
from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import *
import scipy.stats as stat
from scipy.optimize import curve_fit,leastsq
from mpl_toolkits.mplot3d import Axes3D
from astropy.time import Time
from PyAstronomy.pyasl import *

def helcorr(obs_long, obs_lat, obs_alt, ra2000, dec2000, jd, debug=False):
  """
    Calculate barycentric velocity correction.
    
    This function calculates the motion of an observer in
    the direction of a star. In contract to :py:func:`baryvel`
    and :py:func:`baryCorr`, the rotation of the Earth is
    taken into account.
    
    .. note:: This function was ported from the REDUCE IDL package.
              See Piskunov & Valenti 2002, A&A 385, 1095 for a detailed
              description of the package and/or visit
              http://www.astro.uu.se/~piskunov/RESEARCH/REDUCE/
    
    .. warning:: Contrary to the original implementation the longitude
                 increases toward the East and the right ascension is
                 given in degrees instead of hours. The JD is given as is,
                 in particular, nothing needs to be subtracted.
    
    Parameters
    ----------
    obs_long : float
        Longitude of observatory (degrees, **eastern** direction is positive)
    obs_lat : float
        Latitude of observatory [deg]
    obs_alt : float
        Altitude of observatory [m]
    ra2000 : float
        Right ascension of object for epoch 2000.0 [deg]
    dec2000 : float
        Declination of object for epoch 2000.0 [deg]
    jd : float
        Julian date for the middle of exposure.
    
    Returns
    -------
    Barycentric correction : float
        The barycentric correction accounting for the rotation
        of the Earth, the rotation of the Earth's center around
        the Earth-Moon barycenter, and the motion of the Earth-Moon 
        barycenter around the center of the Sun [km/s].
    HJD : float
        Heliocentric Julian date for middle of exposure.

    Notes
    -----

    :IDL REDUCE - Documentation:


    Calculates heliocentric Julian date, barycentric and heliocentric radial
    velocity corrections from:
    
    INPUT:
    <OBSLON> Longitude of observatory (degrees, western direction is positive)
    <OBSLAT> Latitude of observatory (degrees)
    <OBSALT> Altitude of observatory (meters)
    <RA2000> Right ascension of object for epoch 2000.0 (hours)
    <DE2000> Declination of object for epoch 2000.0 (degrees)
    <JD> Julian date for the middle of exposure
    [DEBUG=] set keyword to get additional results for debugging
    
    OUTPUT:
    <CORRECTION> barycentric correction - correction for rotation of earth,
       rotation of earth center about the earth-moon barycenter, earth-moon 
       barycenter about the center of the Sun.
    <HJD> Heliocentric Julian date for middle of exposure
    
    Algorithms used are taken from the IRAF task noao.astutils.rvcorrect
    and some procedures of the IDL Astrolib are used as well.
    Accuracy is about 0.5 seconds in time and about 1 m/s in velocity.
    
    History:
    written by Peter Mittermayer, Nov 8,2003
    2005-January-13   Kudryavtsev   Made more accurate calculation of the sidereal time.
                                    Conformity with MIDAS compute/barycorr is checked.
    2005-June-20      Kochukhov Included precession of RA2000 and DEC2000 to current epoch

"""
  from PyAstronomy.pyaC import degtorad

  # This reverts the original longitude convention. After this,
  # East longitudes are positive
  obs_long = -obs_long

  if jd < 2.4e6:
    PE.warn(PE.PyAValError("The given Julian Date (" + str(jd) + ") is exceedingly small. Did you subtract 2.4e6?"))

  # Covert JD to Gregorian calendar date
  xjd = jd
  
  year, month, day, ut = tuple(daycnv(xjd))

  # Current epoch
  epoch = year + month/12. + day/365.

  # Precess ra2000 and dec2000 to current epoch, resulting ra is in degrees
  ra = ra2000
  dec = dec2000
  ra, dec = precess(ra, dec, 2000.0, epoch)  

  # Calculate heliocentric julian date
  rjd = jd-2.4e6
  hjd = helio_jd(rjd, ra, dec) + 2.4e6

  # DIURNAL VELOCITY (see IRAF task noao.astutil.rvcorrect)
  # convert geodetic latitude into geocentric latitude to correct
  # for rotation of earth
  dlat = -(11.*60.+32.743)*np.sin(2.0*degtorad(obs_lat)) \
         +1.1633*np.sin(4.0*degtorad(obs_lat)) - 0.0026*np.sin(6.0*degtorad(obs_lat))
  lat = obs_lat + dlat/3600.0

  # Calculate distance of observer from earth center
  r = 6378160.0 * (0.998327073+0.001676438*np.cos(2.0*degtorad(lat)) \
     -0.00000351 * np.cos(4.0*degtorad(lat)) + 0.000000008*np.cos(6.0*degtorad(lat))) \
     + obs_alt

  # Calculate rotational velocity (perpendicular to the radius vector) in km/s
  # 23.934469591229 is the sidereal day in hours for 1986
  v = 2.*np.pi * (r/1000.) / (23.934469591229*3600.)

  # Calculating local mean sidereal time (see astronomical almanach)
  tu = (rjd-51545.0)/36525.0
  gmst = 6.697374558 + ut + \
        (236.555367908*(rjd-51545.0) + 0.093104*tu**2 - 6.2e-6*tu**3)/3600.0
  lmst = idlMod(gmst-obs_long/15, 24)

  # Projection of rotational velocity along the line of sight
  vdiurnal = v*np.cos(degtorad(lat))*np.cos(degtorad(dec))*np.sin(degtorad(ra-lmst*15))

  # BARICENTRIC and HELIOCENTRIC VELOCITIES
  vh, vb = baryvel(xjd,0)

  # Project to line of sight
  vbar = vb[0]*np.cos(degtorad(dec))*np.cos(degtorad(ra)) + vb[1]*np.cos(degtorad(dec))*np.sin(degtorad(ra)) + \
         vb[2]*np.sin(degtorad(dec))
  vhel = vh[0]*np.cos(degtorad(dec))*np.cos(degtorad(ra)) + vh[1]*np.cos(degtorad(dec))*np.sin(degtorad(ra)) + \
         vh[2]*np.sin(degtorad(dec))
  
  # Use barycentric velocity for correction
  corr = (vdiurnal + vbar) 

  if debug:
    print ''
    print '----- HELCORR.PRO - DEBUG INFO - START ----'
    print '(obs_long (East positive),obs_lat,obs_alt) Observatory coordinates [deg,m]: ', -obs_long, obs_lat, obs_alt
    print '(ra,dec) Object coordinates (for epoch 2000.0) [deg]: ', ra,dec
    print '(ut) Universal time (middle of exposure) [hrs]: ', ut
    print '(jd) Julian date (middle of exposure) (JD): ', jd
    print '(hjd) Heliocentric Julian date (middle of exposure) (HJD): ', hjd
    print '(gmst) Greenwich mean sidereal time [hrs]: ', idlMod(gmst, 24)
    print '(lmst) Local mean sidereal time [hrs]: ', lmst
    print '(dlat) Latitude correction [deg]: ', dlat
    print '(lat) Geocentric latitude of observer [deg]: ', lat
    print '(r) Distance of observer from center of earth [m]: ', r
    print '(v) Rotational velocity of earth at the position of the observer [km/s]: ', v
    print '(vdiurnal) Projected earth rotation and earth-moon revolution [km/s]: ', vdiurnal
    print '(vbar) Barycentric velocity [km/s]: ', vbar
    print '(vhel) Heliocentric velocity [km/s]: ', vhel
    print '(corr) Vdiurnal+vbar [km/s]: ', corr
    print '----- HELCORR.PRO - DEBUG INFO - END -----'
    print ''
  
  return corr, hjd


#Read in CPS wavelength solutions for green and red chips
wv_green=genfromtxt('wvscale_green.csv',delimiter=',')
wv_red=genfromtxt('wvscale_red.csv',delimiter=',')
#Stack chips onto one array
wave=vstack(([wv_green,wv_red]))

cpstemp={'K5':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'K7':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M0':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M1':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M1.5':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M2':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M2.5':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M3':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.},
        'M4':{'f':zeros(4021),'w':wave,'rv':0.,'helio':0.}}

print 'Reading CPS templates...'
#Read in RV templates
img1=pf.open('HIRES-CPStemplates/rj59.2218.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij59.2218.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['K5']['f']=flux
cpstemp['K5']['h']=header
cpstemp['K5']['name']='HD 204587'
cpstemp['K5']['rv']=-84.533
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj38.235.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij38.235.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['K7']['f']=flux
cpstemp['K7']['h']=header
cpstemp['K7']['name']='HD 201092'
cpstemp['K7']['rv']=-64.07
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj72.438.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij72.438.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M0']['f']=flux
cpstemp['M0']['h']=header
cpstemp['M0']['name']='SAO 122446'
cpstemp['M0']['rv']=-12.51
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj72.168.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij72.168.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M1']['f']=flux
cpstemp['M1']['h']=header
cpstemp['M1']['name']='HD 165222'
cpstemp['M1']['rv']=32.67
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj61.725.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij61.725.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M1.5']['f']=flux
cpstemp['M1.5']['h']=header
cpstemp['M1.5']['name']='GL 87'
cpstemp['M1.5']['rv']=-2.99
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj82.722.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij82.722.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M2']['f']=flux
cpstemp['M2']['h']=header
cpstemp['M2']['name']='HIP 8051'
cpstemp['M2']['rv']=-25.88
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj59.993.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij59.993.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M2.5']['f']=flux
cpstemp['M2.5']['h']=header
cpstemp['M2.5']['name']='HD 180617'
cpstemp['M2.5']['rv']=35.737
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj53.291.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij53.291.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M3']['f']=flux
cpstemp['M3']['h']=header
cpstemp['M3']['name']='GL 687'
cpstemp['M3']['rv']=-28.58
img1.close()
img2.close()

img1=pf.open('HIRES-CPStemplates/rj45.605.fits') #green chip
img2=pf.open('HIRES-CPStemplates/ij45.605.fits') #red chip
green=img1[0].data
red=img2[0].data
header=img1[0].header
flux=vstack(([green,red]))
cpstemp['M4']['f']=flux
cpstemp['M4']['h']=header
cpstemp['M4']['name']='GL 876'
cpstemp['M4']['rv']=-1.59
img1.close()
img2.close()

def cps_shift_estimate(sci_path,spty1,spty2):
    #Read in science file from sci_path
    science={'f':zeros(4021),
            'w':wave,
            'helio':0.}
    str1=sci_path[:9]
    str2=sci_path[10:]
    img1=pf.open('CPS/'+'r'+str2) #green chip
    img2=pf.open('CPS/'+'i'+str2) #red chip
    green=img1[0].data
    red=img2[0].data
    header=img1[0].header
    flux=vstack(([green,red]))
    science['f']=flux
    science['h']=header
    img1.close()
    img2.close()
    
    #Create dictionary with science and both template spectra
    temp1={'w':zeros([26,4021]),'f':zeros([26,4021]),'h':cpstemp[spty1]['h'],\
        'rv':cpstemp[spty1]['rv']}
    temp2={'w':zeros([26,4021]),'f':zeros([26,4021]),'h':cpstemp[spty2]['h'],\
        'rv':cpstemp[spty2]['rv']}
    sci={'w':zeros([26,4021]),'f':zeros([26,4021]),'h':science['h']}
    for i in range(26):
        for j in range(4021):
            temp1['w'][i,j]=cpstemp[spty1]['w'][i,j]
            temp1['f'][i,j]=cpstemp[spty1]['f'][i,j]
            temp2['w'][i,j]=cpstemp[spty2]['w'][i,j]
            temp2['f'][i,j]=cpstemp[spty2]['f'][i,j]
            sci['w'][i,j]=science['w'][i,j]
            sci['f'][i,j]=science['f'][i,j]
    images={'science':sci,'temp1':temp1,'temp2':temp2}

    #Calculate HJD of middle of observation and heliocentric correction using
    #coordinates of Keck Observatory
    lat=19.828
    lon=-155.478
    alt=4160
    for image in images:
        times=[images[image]['h']['DATE_BEG'],images[image]['h']['DATE_END']]
        ra=images[image]['h']['RA'].split(':')
        dec=images[image]['h']['DEC'].split(':')
        ra=(float(ra[0])+float(ra[1])/60+float(ra[2])/3600)*15
        if dec[0][0]!='-':
            dec=float(dec[0])+float(dec[1])/60+float(dec[2])/3600
        else:
            dec=float(dec[0])-float(dec[1])/60-float(dec[2])/3600
        times[0]=times[0].replace('T',' ')
        times[1]=times[1].replace('T',' ')
        times=Time(times,format='iso',scale='utc')
        jd=(times.jd[0]+times.jd[1])/2
        bv,hjd=helcorr(lon,lat,alt,ra,dec,jd)
        hjd2=helio_jd(jd-2.4e6,ra,dec)
        images[image]['helio']=bv
        images[image]['hjd']=hjd
     
    r=1.28/2.998E5 #Use velocity spacing of 1.28 km/s/pixel
    s=zeros(1)
    t1=zeros(4021)
    t2=zeros(4021)
    for i in range(len(science['f'])):
        #Apply heliocentric correction
        for image in images:
            nflux,nwave=dopplerShift(images[image]['w'][i], \
                images[image]['f'][i],images[image]['helio'], \
                edgeHandling='firstlast')
            images[image]['f'][i]=nflux
        #Flatten spectra using 3rd order polynomial
        ps=polyfit(images['science']['w'][i],images['science']['f'][i],3)
        pt1=polyfit(images['temp1']['w'][i],images['temp1']['f'][i],3)
        pt2=polyfit(images['temp2']['w'][i],images['temp2']['f'][i],3)
    ##    ps2=polyfit(sci['w'][i],sci['f'][i],2)
    ##    pt12=polyfit(temp1['w'][i],temp1['f'][i],2)
    ##    pt22=polyfit(temp2['w'][i],temp2['f'][i],2)
    ##    ps3=polyfit(sci['w'][i],sci['f'][i],3)
    ##    pt13=polyfit(temp1['w'][i],temp1['f'][i],3)
    ##    pt23=polyfit(temp2['w'][i],temp2['f'][i],3)
    ##    ps7=polyfit(sci['w'][i],sci['f'][i],7)
    ##    pt17=polyfit(temp1['w'][i],temp1['f'][i],7)
    ##    pt27=polyfit(temp2['w'][i],temp2['f'][i],7)
    ##    fit5=[ps[0]*ele**5+ps[1]*ele**4+ps[2]*ele**3+ps[3]*ele**2+ps[4]*ele+ps[5] \
    ##        for ele in science['w'][i]]
    ##    fit2=[ps2[0]*ele**2+ps2[1]*ele+ps2[2] for ele in science['w'][i]]
    ##    fit3=[ps3[0]*ele**3+ps3[1]*ele**2+ps3[2]*ele+ps3[3] for ele in \
    ##        science['w'][i]]
    ##    fit7=[ps7[0]*ele**7+ps7[1]*ele**6+ps7[2]*ele**5+ps7[3]*ele**4+ps7[4]*ele**3\
    ##        +ps7[5]*ele**2+ps7[6]*ele+ps7[7] for ele in science['w'][i]]
        for j in range(4021):
    ##        #5th order polynomial
    ##        sci['f'][i,j]/=ps[0]*sci['w'][i,j]**5+ps[1]*sci['w'][i,j]**4+ps[2]* \
    ##            sci['w'][i,j]**3+ps[3]*sci['w'][i,j]**2+ps[4]*sci['w'][i,j]+ps[5]
    ##        temp1['f'][i,j]/=pt1[0]*temp1['w'][i,j]**5+pt1[1]*temp1['w'][i,j]**4+ \
    ##            pt1[2]*temp1['w'][i,j]**3+pt1[3]*temp1['w'][i,j]**2+pt1[4]*\
    ##            temp1['w'][i,j]+pt1[5]
    ##        temp2['f'][i,j]/=pt2[0]*temp2['w'][i,j]**5+pt2[1]*temp2['w'][i,j]**4+ \
    ##            pt2[2]*temp2['w'][i,j]**3+pt2[3]*temp2['w'][i,j]**2+pt2[4]*\
    ##            temp2['w'][i,j]+pt2[5]
            
            #3rd order polynomial 
            images['science']['f'][i,j]/=ps[0]*images['science']['w'][i,j]**3+\
                ps[1]*images['science']['w'][i,j]**2+ps[2]* \
                images['science']['w'][i,j]+ps[3]
            images['temp1']['f'][i,j]/=pt1[0]*images['temp1']['w'][i,j]**3+ \
                pt1[1]*images['temp1']['w'][i,j]**2+pt1[2]* \
                images['temp1']['w'][i,j]+pt1[3]
            images['temp2']['f'][i,j]/=pt2[0]*images['temp2']['w'][i,j]**3+ \
                pt2[1]*images['temp2']['w'][i,j]**2+pt2[2]* \
                images['temp2']['w'][i,j]+pt2[3]
            
    ##        #2nd order polynomial
    ##        sci['f'][i,j]/=ps[0]*sci['w'][i,j]**2+ps[1]*sci['w'][i,j]+ps[2]
    ##        temp1['f'][i,j]/=pt1[0]*temp1['w'][i,j]**2+pt1[1]*temp1['w'][i,j]+pt1[2]
    ##        temp2['f'][i,j]/=pt2[0]*temp2['w'][i,j]**2+pt2[1]*temp2['w'][i,j]+pt2[2]

            #Normalize to 0 continuum, rather than 1
            images['science']['f'][i,j]-=1
            images['temp1']['f'][i,j]-=1
            images['temp2']['f'][i,j]-=1
            
        #Interpolate to common logarithmic wavelength scale    
        fs=interp1d(images['science']['w'][i],images['science']['f'][i], \
            fill_value=0.,bounds_error=False)
        ft1=interp1d(images['temp1']['w'][i],images['temp1']['f'][i], \
            fill_value=0.,bounds_error=False)
        ft2=interp1d(images['temp2']['w'][i],images['temp2']['f'][i], \
            fill_value=0.,bounds_error=False)
        wv=images['science']['w'][i,0]*(1.+r)**linspace(0.,4020.,4021)
        if len(s)>1:
            s=vstack((s,fs(wv)))
            t1=vstack((t1,ft1(wv)))
            t2=vstack((t2,ft2(wv)))
        elif len(s)==1:
            s=fs(wv)
            t1=ft1(wv)
            t2=ft2(wv)

    #Calculate normalized 1-D cross-correlation functions for each order
    c1=zeros([len(s),8041])
    c2=zeros([len(s),8041])
    c12=zeros([len(s),8041])
    shift=range(-4020,4021)
    for i in range(len(s)):
        #CCF of science with primary template
        c1[i]=correlate(s[i],t1[i],mode='full')/ \
        (sqrt(mean(s[i]**2))*sqrt(mean(t1[i]**2))*len(s[i]))

        #CCF of science with secondary template 
        c2[i]=correlate(s[i],t2[i],mode='full')/ \
         (sqrt(mean(s[i]**2))*sqrt(mean(t2[i]**2))*len(s[i]))

        #CCF between templates 
        c12[i]=correlate(t1[i],t2[i],mode='full')/ \
          (sqrt(mean(t1[i]**2))*sqrt(mean(t2[i]**2))*len(t1[i]))
        
    return c1,c2,c12,images

def cps_todcor(c1,c2,c12,pshift,sshift,images):
    def gauss(x,*p):
        A,mu,sigma=p
        return A*exp(-(x-mu)**2/(2.*sigma**2))
    def gauss2D(height_x,center_x,width_x,height_y,center_y,width_y):
        return lambda x,y: height_x*exp(-((center_x-x)/width_x)**2/2)+ \
           height_y*exp(-((center_y-y)/width_y)**2/2)
    #Calculate velocity shift array using 1.28 km/s/pixel        
    v=range(-4020,4021)
    v=[1.28*vel for vel in v]

    #Calculate telluric offset using order near 7500 A
    delta1=list(c1[23]).index(max(c1[23]))
    delta2=list(c2[23]).index(max(c2[23]))
    delta1=v[delta1]
    delta2=v[delta2]

    #Account for differences in heliocentric velocities
    delta1-=images['science']['helio']-images['temp1']['helio']
    delta2-=images['science']['helio']-images['temp2']['helio']
    print 'cps deltas =',delta1,delta2

    #Remove orders that provide unreliabe RVs
    orders=list(images['science']['w'][:,0])
    for wave in orders:
        if abs(wave-5801)<15 or abs(wave-5898)<15 or abs(wave-6553)<15 or \
        abs(wave-6797)<15 or abs(wave-7213)<15 or abs(wave-7520)<15:
            i=orders.index(wave)
            c1=vstack(([c1[:i],c1[i+1:]]))
            c2=vstack(([c2[:i],c2[i+1:]]))
            c12=vstack(([c12[:i],c12[i+1:]]))
            orders.pop(i)

    print 'Calculating TODCOR...'
    #Calculate TODCOR
    R=zeros(len(orders),dtype=list)
    m=zeros([len(orders),3])
    sp=int(pshift)
    ss=int(sshift)
    #Set size of todcor box in velocity space (needs to be an odd number).
    #Processing time scales as size^2. E.g. vrange=79 will search in a 50 km/s
    #by 50 km/s box centered on the inital guess values for the primary and
    #secondary shifts.
    vrange=79#157 
    for k in range(len(R)):
        R[k]=zeros([vrange,vrange])
        for i in range(vrange):
            s2=i+4020+ss-(vrange-1)/2
            for j in range(vrange):
                s1=j+4020+sp-(vrange-1)/2
                R[k][i,j]=sqrt((c1[k,s1]**2-2*c1[k,s1]*c2[k,s2]* \
                    c12[k,s2-s1+4020]+c2[k,s2]**2)/(1-c12[k,s2-s1+4020]**2))
                if R[k][i,j]>m[k,2]:
                    m[k]=[i,j,R[k][i,j]]
        R[k][isnan(R[k])]=0.
        R[k][isinf(R[k])]=0.
        
    
    print 'Calculating velocity errors...'
    #Convert to velocity units
    v1=images['temp1']['rv']
    v2=images['temp2']['rv']
    vp=range(-(vrange-1)/2,(vrange-1)/2+1)
    vs=range(-(vrange-1)/2,(vrange-1)/2+1)
    for i in range(len(vp)):
        vp[i]=round((vp[i]+sp)*1.28,2)+v1
        vs[i]=round((vs[i]+ss)*1.28,2)+v2

    fwhm={'p':[],'s':[]}
    peak={'p':[],'s':[]}
    error={'p':[],'s':[]}
    center={'p':[],'s':[]}
    for i in range(len(m)):
        p1=[.05,v1+(sp-(vrange-1)/2+m[i,1])*1.28,10]
        p2=[.05,v2+(ss-(vrange-1)/2+m[i,0])*1.28,10]
        try:
            coeff1,var_matrix1=curve_fit(gauss,vp,R[i][m[i,0],:]- \
                min(R[i][m[i,0],:]),p0=p1)
            coeff2,var_matrix2=curve_fit(gauss,vs,R[i][:,m[i,1]]- \
                min(R[i][:,m[i,1]]),p0=p2)
        except RuntimeError:
            fwhm['p'].append('')
            fwhm['s'].append('')
            peak['p'].append('')
            peak['s'].append('')
            error['p'].append('')
            error['s'].append('')
            center['p'].append('')
            center['s'].append('')
##            X,Y=meshgrid(vp,vs)
##            hf=figure()
##            ha=hf.add_subplot(111,projection='3d')
##            ha.plot_surface(X,Y,R[i],cmap=cm.cubehelix)
##            xlabel('Primary Velocity')
##            ylabel('Secondary Velocity')
##            title(r'Order starting $\lambda$ = %f A' % orders[i])
##            show()
        else:
            fwhm['p'].append(2.355*coeff1[2])
            fwhm['s'].append(2.355*coeff2[2])
            peak['p'].append(coeff1[0])
            peak['s'].append(coeff2[0])
            center['p'].append(coeff1[1])
            center['s'].append(coeff2[1])
            fit1=gauss(vp,*coeff1)
            fit2=gauss(vs,*coeff2)
            #plot(vp,R[i][m[i,0],:]-min(R[i][m[i,0],:]),'b')
            #plot(vp,fit1,'--r',label='fit')
            #xlabel('Primary Velocity (km/s)')
            #title(r'Order starting $\lambda$ = %f A' % orders[i])
            #ylabel('Correlation')
            #legend(loc='best')
            #print m[8]
            #savefig('ccf_fit.pdf')
            #show()
            #plot(vs,R[i][:,m[i,1]]-min(R[i][:,m[i,1]]),'b')
            #plot(vs,fit2,'--r',label='fit')
            #xlabel('Secondary Velocity (km/s)')
            #ylabel('Correlation')
            #show()
            #close()
            out1=list(R[i][m[i,0],:m[i,1]-17])+list(R[i][m[i,0],m[i,1]+18:])
            out2=list(R[i][:m[i,0]-17,m[i,1]])+list(R[i][m[i,0]+18:,m[i,1]])
            std1=std(out1)
            std2=std(out2)
            r1=peak['p'][-1]/(sqrt(2)*std1)
            r2=peak['s'][-1]/(sqrt(2)*std2)
            sig1=(3./8.)*fwhm['p'][-1]/(1+r1)
            sig2=(3./8.)*fwhm['s'][-1]/(1+r2)
            error['p'].append(sig1)
            error['s'].append(sig2)
            fits={'fwhm':fwhm,'amp':peak,'center':center}
            p=[peak['p'][-1],center['p'][-1], \
               fwhm['p'][-1]/2.355,peak['s'][-1], \
               center['s'][-1],fwhm['s'][-1]/2.355]
##            #Plot TODCOR
##            X,Y=meshgrid(vp,vs)
##            Z=zeros([vrange,vrange])
##            f=gauss2D(*     p)
##            for k in range(vrange):
##                for j in range(vrange):
##                    Z[k,j]=f(X[k,j],Y[k,j])
##            hf=figure()
##            ha=hf.add_subplot(111,projection='3d')
##            ha.plot_surface(X,Y,R[i],cmap=cm.cubehelix)
##            ha.plot_wireframe(X,Y,Z)
##            xlabel('Primary Velocity')
##            ylabel('Secondary Velocity')
##            title(r'Order starting at $\lambda$ = %f A' % orders[i])
##            show()
##    params=[.4,m[1,1],fwhm['p'][1]/2.355,.3,m[1,0],fwhm['s'][1]/2.355]
##    errorfunction = lambda p: ravel(gauss2D(*p)(*indices(shape(R[1])))-R[1])
##    p,success = leastsq(errorfunction,params)

##    p0=[.2,m[1,1],10.,.1]
##    coeff,cov=curve_fit(gaussc,range(-78,79),R[1][20,:],p0=p0)
##    print coeff
##    print cov
       
        
        
##    ylab=[]
##    xlab=[]
##    for i in range(8):
##        ylab.append(v2+round((ss-(vrange-1)/2+20*i)*1.28,1))
##        xlab.append(v1+round((sp-(vrange-1)/2+20*i)*1.28,1))
##    
##    
##    fig=figure()
##    ax1=fig.add_subplot(111)
##    index=orders.index(5618.333593)
##    ax1.imshow(R[index],cmap=cm.cubehelix)
##    im=ax1.imshow(R[index],cmap=cm.cubehelix)
##    fig.colorbar(im)
##    ax1.set_yticks([0,20,40,60,80,100,120,140])
##    ax1.set_yticklabels(ylab)
##    ax1.set_xticks([0,20,40,60,80,100,120,140])
##    ax1.set_xticklabels(xlab)
##    xlabel(r'$v_1$ (km/s)')
##    ylabel(r'$v_2$ (km/s)')
##    scatter(m[index,1],m[index,0],c='k',marker='+')
##    title('Colormap Plot of R')
##    savefig('TODCOR_Plot2.pdf')
##    close()

    #Calculate centroid of R
    print 'Calculating centroid...'
    centroids=zeros([len(R),4])
    for i in range(len(R)):
        columns=zeros(len(R[i]))
        rows=zeros(len(R[i]))
        for j in range(len(R[i])):
            rows[j]=sum(R[i][j])
            columns[j]=sum(R[i][:,j])
            centroids[i,0]+=columns[j]*j
            centroids[i,1]+=rows[j]*j
        centroids[i,0]/=sum(columns)
        centroids[i,1]/=sum(rows)
        centroids[i,0]=v1+(sp-(vrange-1)/2+centroids[i,0])*1.28-delta1
        centroids[i,1]=v2+(ss-(vrange-1)/2+centroids[i,1])*1.28-delta2
        
##    f=open('rv_results','a')
    g=open('current_rvs','a')
##    print >>f, images['science']['target'],images['science']['hjd'],1,1
    print >>g, images['science']['h']['targname'],images['science']['hjd'],1,1           
    for i in range(len(m)):
        if error['p'][i]!='':
            centroids[i,2]=error['p'][i]
            centroids[i,3]=error['s'][i]
            m[i,1]=v1+(sp-(vrange-1)/2+m[i,1])*1.28-delta1
            m[i,0]=v2+(ss-(vrange-1)/2+m[i,0])*1.28-delta2
##            print >>f, orders[i],m[i,1],error['p'][i],m[i,0],error['s'][i]
##            print >>g, orders[i],center['p'][i],error['p'][i],m[i,0], \
##                  error['s'][i]           
##            print >>g, orders[i],center['p'][i]-delta1,error['p'][i], \
##                  center['s'][i]-delta2,error['s'][i]

##    print >>g, 'Centroids'
    for i in range(len(m)):
        if error['p'][i]!='':
            print >>g, orders[i],centroids[i,0],centroids[i,2],centroids[i,1], \
                  centroids[i,3]
##    f.close()
    g.close()

    return R,m,fits,vp,vs,centroids
    ##a1=stat.nanmedian(m[:,1])
    ##a2=stat.nanmedian(m[:,0])
    ##scatter(orders,m[:,1],c='b',label='Primary')
    ##scatter(orders,m[:,0],c='r',label='Seconday')
    ##hlines(a1,0,8000,colors='blue')
    ##hlines(a2,0,8000,colors='red')
    ##legend(loc='lower left')
    ##xlabel('Order Starting Wavelength')
    ##ylabel('Velocity (km/s')
    ##xlim(3000,8000)
    ##savefig('MG78457_velocities.pdf')
    ##close()

    ###Filter out bad orders (velocity > 2 sigma from median)
    ##for j in range(2):
    ##    i=0 
    ##    while i<len(m):
    ##        a1=stat.nanmedian(m[:,1])
    ##        a2=stat.nanmedian(m[:,0])
    ##        std1=std(m[:,1])
    ##        std2=std(m[:,0])
    ##        if abs(m[i,0]-a2)/std2>2:
    ##            m=vstack((m[:i],m[i+1:]))
    ##            orders.pop(i)
    ##        elif abs(m[i,1]-a1)/std1>2:
    ##            m=vstack((m[:i],m[i+1:]))
    ##            orders.pop(i)
    ##        else:
    ##            i+=1
    ##
    ##scatter(orders,m[:,1],c='b',label='Primary')
    ##scatter(orders,m[:,0],c='r',label='Seconday')
    ##hlines(a1,0,8000,colors='blue')
    ##hlines(a2,0,8000,colors='red')
    ##legend(loc='lower left')
    ##xlabel('Order Starting Wavelength')
    ##ylabel('Velocity (km/s')
    ##xlim(3000,8000)
    ##savefig('MG78457_Fixed_velocities.pdf')
    ##close()
    ##
    ##        
    ##print 'v_1=%.2f +/- %.2f km/s' % (a1,std1)
    ##print 'v_2=%.2f +/- %.2f km/s' % (a2,std2)


    
