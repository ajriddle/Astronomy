import pyfits as pf
from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import *
import scipy.stats as stat
from scipy.optimize import curve_fit,leastsq
from mpl_toolkits.mplot3d import Axes3D
from astropy.time import Time
from PyAstronomy.pyasl import *


def rmspace(var):
    while '' in var:
        var.pop(var.index(''))
    return var

def gauss(x,*p):
    A,mu,sigma=p
    return A*exp(-(x-mu)**2/(2.*sigma**2))

def gaussc(x,*p):
    A,mu,sigma,c=p
    return A*exp(-(x-mu)**2/(2.*sigma**2))+c

def gauss2D(height_x,center_x,width_x,height_y,center_y,width_y):
    return lambda x,y: height_x*exp(-((center_x-x)/width_x)**2/2)+ \
           height_y*exp(-((center_y-y)/width_y)**2/2)

print 'Reading archive templates...'
templates={'K7':{'f':genfromtxt('Archive_Templates/HI.20101030.44526_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20101030.44526_'\
                'temp1wave',delimiter=','),
                'rv':62.6,
                'helio':genfromtxt('Archive_Templates/HI.20101030.44526'\
                '_temp1_helio'),
                'name':'GJ 156'},
           'M0':{'f':genfromtxt('Archive_Templates/HI.20080120.54239_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20080120.54239_'\
                'temp1wave',delimiter=','),
                'rv':-13.9,
                'name':'HD 95650',
                'helio':genfromtxt('Archive_Templates/HI.20080120.54239'\
                '_temp1_helio')},
           'M0.5':{'f':genfromtxt('Archive_Templates/HI.20101030.40829_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20101030.40829_'\
                'temp1wave',delimiter=','),
                'rv':-37.9,
                'name':'GJ 96',
                'helio':genfromtxt('Archive_Templates/HI.20101030.40829'\
                '_temp1_helio')},
           'M1':{'f':genfromtxt('Archive_Templates/HI.20080811.39891_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20080811.39891_'\
                'temp1wave',delimiter=','),
                'rv':-71.1,
                'name':'GJ 908',
                'helio':genfromtxt('Archive_Templates/HI.20080811.39891'\
                '_temp1_helio')},
           'M1.5':{'f':genfromtxt('Archive_Templates/HI.20080120.42565_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20080120.42565_'\
                'temp1wave',delimiter=','),
                'rv':8.7,
                'name':'HD 36395',
                'helio':genfromtxt('Archive_Templates/HI.20080120.42565'\
                '_temp1_helio')},
           'M2.5':{'f':genfromtxt('Archive_Templates/HI.20090704.39575_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20090704.39575_'\
                'temp1wave',delimiter=','),
                'rv':35.737,
                'name':'HD 180617',
                'helio':genfromtxt('Archive_Templates/HI.20090704.39575'\
                '_temp1_helio')},
           'M3':{'f':genfromtxt('Archive_Templates/HI.20080810.34177_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20080810.34177_'\
                'temp1wave',delimiter=','),
                'rv':-.8,
                'name':'HD 173739',
                'helio':genfromtxt('Archive_Templates/HI.20080810.34177'\
                '_temp1_helio')},
           'M3.5':{'f':genfromtxt('Archive_Templates/HI.20101030.41219_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20101030.41219_'\
                'temp1wave',delimiter=','),
                'rv':30.57,
                'name':'GJ 109',
                'helio':genfromtxt('Archive_Templates/HI.20101030.41219'\
                '_temp1_helio')},
           'M4':{'f':genfromtxt('Archive_Templates/HI.20080120.54549_'\
                'temp1flux',delimiter=','),
                'w':genfromtxt('Archive_Templates/HI.20080120.54549_'\
                'temp1wave',delimiter=','),
                'rv':-31.1,
                'name':'GJ 447',
                'helio':genfromtxt('Archive_Templates/HI.20080120.54549'\
                '_temp1_helio')}}

def archive_shift_estimate(str1,str2,spty1,spty2):
#def calcrv(str1,str2,rv,spt):    
    #Define image dictionary
    #1, 2, and 3 are blue, green, and red chips, respectively
    #f and w are the flux and wavelength arrays for all orders and chips stacked

    images={'science':{'1':{}, 
                       '2':{}, 
                       '3':{}, 
                       'f':zeros(4056),
                       'w':zeros(4056)},
            'temp1':templates[spty1],
            'temp2':templates[spty2]}

    #Set paths for science and template images
    images['science']['folder']=str1
    images['science']['file']=str2

    
    key='science'
    try:
        images[key]['w']=genfromtxt('HIRES/'+images['science']['file']+key+ \
            'wave',delimiter=',')
        
        images[key]['f']=genfromtxt('HIRES/'+images['science']['file']+key+ \
            'flux',delimiter=',')
        
        header=genfromtxt('HIRES/'+images['science']['file']+'header', \
            delimiter=',',dtype=str)

        helio=genfromtxt('HIRES/'+images['science']['file']+key+'_helio')
    except IOError:
        print 'Reading science file...'
        for s in ['1','2','3']:
            path='HIRES/'+images[key]['file']+s+'_Flux.fits'
            try:
                img=pf.open(path)
                string=img[0].header['WV_0_01']
            except (IOError,KeyError):
                images[key][s]['skip']='yes'
            else:
                images[key][s]['skip']='no'
                images[key][s]['f']=img[0].data
                img.close()
                images[key][s]['h']=img[0].header
                images[key]['h']=img[0].header
                #Read in wavelength solution
                images[key][s]['c']=zeros([len(images[key][s]['f']),7])
                images[key][s]['w']=zeros([len(images[key][s]['f']),4056])
                i=0
                z=0
                keys=images[key][s]['h'].keys()
                while i<len(keys):
                    if 'WV_0' in keys[i]:
                        s1=images[key][s]['h'][keys[i]]
                        s2=images[key][s]['h'][keys[i+1]]
                        coeff1=rmspace(s1.split(' '))
                        coeff2=rmspace(s2.split(' '))
                        images[key][s]['c'][z,0]=float(coeff1[0])
                        images[key][s]['c'][z,1]=float(coeff1[1])
                        images[key][s]['c'][z,2]=float(coeff1[2])
                        images[key][s]['c'][z,3]=float(coeff1[3])
                        images[key][s]['c'][z,4]=float(coeff2[0])
                        images[key][s]['c'][z,5]=float(coeff2[1])
                        images[key][s]['c'][z,6]=float(coeff2[2])
                        i+=2
                        z+=1
                    else:
                        i+=1
                #Apply wavelength solution
                for i in range(4056):
                    for j in range(len(images[key][s]['c'])):
                        images[key][s]['w'][j,i]=images[key][s]['c'][j,0]+ \
                        images[key][s]['c'][j,1]*i+ \
                        images[key][s]['c'][j,2]*i**2+ \
                        images[key][s]['c'][j,3]*i**3+ \
                        images[key][s]['c'][j,4]*i**4+ \
                        images[key][s]['c'][j,5]*i**5+ \
                        images[key][s]['c'][j,6]*i**6
                #Stack orders from all chips onto one array
                images[key]['w']=vstack((images[key]['w'], \
                    images[key][s]['w']))
                images[key]['f']=vstack((images[key]['f'], \
                    images[key][s]['f']))
                try:
                    helio=images[key]['h']['heliovel']
                except KeyError:
                    pass
        images[key]['w']=images[key]['w'][1:]
        images[key]['f']=images[key]['f'][1:]

        savetxt('HIRES/'+images['science']['file'] \
            +key+'wave',images[key]['w'],delimiter=',')
        
        savetxt('HIRES/'+images['science']['file'] \
            +key+'flux',images[key]['f'],delimiter=',')

        f=open('HIRES/'+images['science']['file']+ \
               'header','w')
        #Calculate HJD of middle of observation
        times=[images['science']['h']['DATE_BEG'], \
            images['science']['h']['DATE_END']]
        ra=images['science']['h']['RA2000'].split(':')
        dec=images['science']['h']['DEC2000'].split(':')
        ra=(float(ra[0])+float(ra[1])/60+float(ra[2])/3600)*15
        if float(dec[0])>0:
            dec=float(dec[0])+float(dec[1])/60+float(dec[2])/3600
        else:
            dec=float(dec[0])-float(dec[1])/60-float(dec[2])/3600
        times[0]=times[0].replace('T',' ')
        times[1]=times[1].replace('T',' ')
        times=Time(times,format='iso',scale='utc')
        hjd=(times.jd[0]+times.jd[1])/2-2400000
        hjd=helio_jd(hjd,ra,dec)
        
        header=[images['science']['h']['targname'],hjd]
        print >>f, header[0],',',header[1]
        f.close()

        f=open('HIRES/'+images['science']['file']+key+'_helio','w') 
        print >>f, helio
        f.close()
    else:
        images['science']['target']=header[0]
        images['science']['hjd']=header[1]
        images['science']['helio']=helio
    

    #Continuum normalize each order and interpolate over common wavelength scale
    #for each order
    r=1.28/2.998E5
    temp1={'w':zeros([len(images['temp1']['w']),4056]), \
           'f':zeros([len(images['temp1']['f']),4056])}
    temp2={'w':zeros([len(images['temp2']['w']),4056]), \
           'f':zeros([len(images['temp2']['f']),4056])}
    science={'w':images['science']['w'],'f':images['science']['f']}
    for i in range(len(images['temp1']['f'])):
        for j in range(4056):
            temp1['w'][i,j]=images['temp1']['w'][i,j]
            temp1['f'][i,j]=images['temp1']['f'][i,j]
    for i in range(len(images['temp2']['f'])):
        for j in range(4056):
            temp2['w'][i,j]=images['temp2']['w'][i,j]
            temp2['f'][i,j]=images['temp2']['f'][i,j]  
       
    s=zeros(1)
    t1=zeros(4056)
    t2=zeros(4056)
    orders=[]
    orderdiff=15
    for i in range(len(science['w'])):
        for j in range(len(temp1['w'])):
            if abs(science['w'][i,0]-temp1['w'][j,0])<orderdiff:
                for k in range(len(temp2['w'])):
                    if abs(science['w'][i,0]-temp2['w'][k,0])<orderdiff:
                        ps=polyfit(science['w'][i,:],science['f'][i,:],1)
                        pt1=polyfit(temp1['w'][j,:],temp1['f'][j,:],1)
                        pt2=polyfit(temp2['w'][k,:],temp2['f'][k,:],1)
                        for l in range(4056):
                            science['f'][i,l]/=science['w'][i,l]*ps[0]+ps[1]
                            science['f'][i,l]-=1
                            temp1['f'][j,l]/=temp1['w'][j,l]*pt1[0]+pt1[1]
                            temp1['f'][j,l]-=1
                            temp2['f'][k,l]/=temp2['w'][k,l]*pt2[0]+pt2[1]
                            temp2['f'][k,l]-=1
                        fs=interp1d(science['w'][i,:],science['f'][i,:], \
                            fill_value=0.,bounds_error=False)
                        ft1=interp1d(temp1['w'][j,:],temp1['f'][j,:], \
                            fill_value=0.,bounds_error=False)
                        ft2=interp1d(temp2['w'][k,:],temp2['f'][k,:], \
                            fill_value=0.,bounds_error=False)
                        wv=science['w'][i,0]*(1.+r)**linspace(0.,4055.,4056)
                        orders.append(wv[0])
                        if len(s)>1:
                            s=vstack((s,fs(wv)))
                            t1=vstack((t1,ft1(wv)))
                            t2=vstack((t2,ft2(wv)))
                        elif len(s)==1:
                            s=fs(wv)
                            t1=ft1(wv)
                            t2=ft2(wv)
                        break
                    
    c1=zeros([len(s),8111])
    c2=zeros([len(s),8111])
    c12=zeros([len(s),8111])

    shift=range(-4055,4056)
    for i in range(len(s)):
        c1[i]=correlate(s[i],t1[i],mode='full')/ \
        (sqrt(mean(s[i]**2))*sqrt(mean(t1[i]**2))*len(s[i]))
         
        c2[i]=correlate(s[i],t2[i],mode='full')/ \
         (sqrt(mean(s[i]**2))*sqrt(mean(t2[i]**2))*len(s[i]))
         
        c12[i]=correlate(t1[i],t2[i],mode='full')/ \
          (sqrt(mean(t1[i]**2))*sqrt(mean(t2[i]**2))*len(t1[i]))

    return c1,c2,c12,images,orders

def archive_todcor(c1,c2,c12,pshift,sshift,images,orders):
    #Calculate telluric offset
    for i in range(len(c1)):
        if abs(orders[i]-7529)<15:
            delta=list(c1[i]).index(max(c1[i]))
            break
    v=range(-4055,4056)
    v=[1.28*vel for vel in v]
    
    try:
        delta=v[delta]  
    except (UnboundLocalError,TypeError):
        delta=0.
    else:
        delta-=float(images['science']['helio'])-float(images['temp1']['helio'])
    print delta
    for wave in orders:
        if abs(wave-5801)<15 or abs(wave-5898)<15 or abs(wave-6553)<15 or \
        abs(wave-6797)<15 or abs(wave-7213)<15 or abs(wave-7520)<15 or \
        abs(wave-8033)<15 or abs(wave-8229)<15 or abs(wave-5891)<15:
            i=orders.index(wave)
            c1=vstack(([c1[:i],c1[i+1:]]))
            c2=vstack(([c2[:i],c2[i+1:]]))
            c12=vstack(([c12[:i],c12[i+1:]]))
            orders.pop(i)
        
    print 'Calculating TODCOR...'
    ##shift=range(-4055,4056)
    ##plot(shift,c1[8],'b',label='Primary Shift')
    ##plot(shift,c2[8],'r',label='Seconday Shift')
    ##xlabel('Shift')
    ##xlim(-200,200)
    ##legend(loc='best')
    ##show()
    sp=int(pshift) 
    ss=int(sshift)
    R=zeros(len(orders),dtype=list)
    m=zeros([len(orders),3])
    vrange=79#157 #Set size of todcor box in velocity space; needs to be odd
    for k in range(len(R)):
        R[k]=zeros([vrange,vrange])
        for i in range(vrange):
#            show()
            s2=i+4055+ss-(vrange-1)/2
            for j in range(vrange):
                s1=j+4055+sp-(vrange-1)/2
                R[k][i,j]=sqrt((c1[k,s1]**2-2*c1[k,s1]*c2[k,s2]* \
                    c12[k,s2-s1+4055]+c2[k,s2]**2)/(1-c12[k,s2-s1+4055]**2))
                if R[k][i,j]>m[k,2]:
                    m[k]=[i,j,R[k][i,j]]
        R[k][isnan(R[k])]=0.
        R[k][isinf(R[k])]=0.
    
    print 'Calculating velocity errors...\n'
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
##            f=gauss2D(*p)
##            for k in range(vrange):
##                for j in range(vrange):
##                    Z[k,j]=f(X[k,j],Y[k,j])
##            hf=figure()
##            ha=hf.add_subplot(111,projection='3d')
##            ha.plot_surface(X,Y,R[i],cmap=cm.cubehelix)
##            ha.plot_wireframe(X,Y,Z)
##            xlabel('Primary Velocity')
##            ylabel('Secondary Velocity')
##            title(r'Order starting $\lambda$ = %f A' % orders[i])
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
        centroids[i,0]=v1+(sp-(vrange-1)/2+centroids[i,0])*1.28-delta
        centroids[i,1]=v2+(ss-(vrange-1)/2+centroids[i,1])*1.28-delta
        
    
##    f=open('rv_results','a')
    g=open('current_rvs','a')
##    print >>f, images['science']['target'],images['science']['hjd'],1,1
    print >>g, images['science']['target'],images['science']['hjd'],1,1           
    for i in range(len(m)):
        if error['p'][i]!='':
            centroids[i,2]=error['p'][i]
            centroids[i,3]=error['s'][i]
            m[i,1]=v1+(sp-(vrange-1)/2+m[i,1])*1.28-delta
            m[i,0]=v2+(ss-(vrange-1)/2+m[i,0])*1.28-delta
##            print >>f, orders[i],m[i,1],error['p'][i],m[i,0],error['s'][i]
##            print >>g, orders[i],center['p'][i],error['p'][i],m[i,0], \
##                  error['s'][i]           
##            print >>g, orders[i],center['p'][i]-delta,error['p'][i], \
##                  center['s'][i]-delta,error['s'][i]

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

##if __name__=='__main__':
##    list1=genfromtxt('workfile.csv',dtype=str,delimiter=',')
##    u=zeros(7,dtype=str)
##    for row in list1:
##        u=vstack((u,row[0].split('/')))
##    u=u[1:]
##
##    for q in range(len(u)):
##        g=u[q,1]
##        h=u[q,-1][:18]
##        print 'Starting epoch %i/%i' % (q+1,len(u)) 
##        images,R,m,fits,vp,vs=calcrv(g,h,int(list1[q,6]),int(list1[q,7])) 
    
    
