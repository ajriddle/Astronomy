from numpy import *
import scipy.stats as stat
from EB_Velocity import *
from matplotlib.pyplot import *
from corner import *
from PyAstronomy.pyasl import *

#Read in radial velocities
print 'Reading RV data'
d=genfromtxt('current_rvs',dtype=list)
orders=genfromtxt('order_residuals.txt',dtype=float)
rv={}
epoch=[]
v1=[]
v2=[]
e1=[]
e2=[]

for i in range(len(d)):
    if 'Cnc' in d[i,1]:
        key=d[i,2]
        rv[key]={'v':zeros(5),'m':[0.,0.],'e':[0.,0.],'bo':[]}
    else:
        for j in range(5):
            d[i,j]=float(d[i,j])
        rv[key]['v']=vstack(([rv[key]['v'],d[i]]))
        
#Do two iterations of sigma clipping  
for key in rv:
    rv[key]['v']=rv[key]['v'][1:]
    for i in range(len(rv[key]['v'])):
        for j in range(len(orders)):
            if abs(rv[key]['v'][i,0]-orders[j,0])<15:
                rv[key]['v'][i,2]=orders[j,2]
                rv[key]['v'][i,4]=orders[j,4]
                break
    for j in range(2):
        i=0 
        while i<len(rv[key]['v']):
            a1=stat.nanmedian(list(rv[key]['v'][:,1]))
            a2=stat.nanmedian(list(rv[key]['v'][:,3]))
            if key=='54467.371626':
                a1=81.57
                a2=30.73
            std1=std(rv[key]['v'][:,1])
            std2=std(rv[key]['v'][:,3])
            if abs(rv[key]['v'][i,3]-a2)/std2>2 or abs(rv[key]['v'][i,3]-a2)>30:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))              
            elif abs(rv[key]['v'][i,1]-a1)/std1>2 or  \
            abs(rv[key]['v'][i,1]-a1)>30:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))
            elif abs(rv[key]['v'][i,0]-6553)<10 or abs(rv[key]['v'][i,0]-8229) \
                <10:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))
            else:
                i+=1
    w1=[1/x**2 for x in rv[key]['v'][:,2]]
    w2=[1/x**2 for x in rv[key]['v'][:,4]]
    epoch.append(float(key))
    v1.append(average(list(rv[key]['v'][:,1]),weights=w1))
    v2.append(average(list(rv[key]['v'][:,3]),weights=w2))
    s1=0.
    s2=0.
    for i in range(len(rv[key]['v'])):
        s1+=w1[i]*(rv[key]['v'][i,1]-v1[-1])**2
        s2+=w2[i]*(rv[key]['v'][i,3]-v2[-1])**2
    
    e1.append((s1/sum(w1)/(len(w1)-1))**.5)
    e2.append((s2/sum(w2)/(len(w2)-1))**.5)
    if len(rv[key]['v'])==1:
        e1[-1]=rv[key]['v'][0][2]
        e2[-1]=rv[key]['v'][0][4]
    if key=='54467.371626' or key=='2455851.05523':
        v1[-1],v2[-1]=v2[-1],v1[-1]
        e1[-1],e2[-1]=e2[-1],e1[-1]
    if key=='54690.482060':
        v2[-1]=-56.26
    rv[key]['m']=[v1[-1],v2[-1]]
    rv[key]['e']=[e1[-1],e2[-1]]

#Remove epochs that are too close to eclipse (peak pulling can cause spurious
#velocity measurements when the correlation peaks are too close)
mvs=40 #Set minimum velocity separation allowed
i=0
while i<len(v1):
    if abs(v1[i]-v2[i])<mvs:
        v1.pop(i)
        v2.pop(i)
        e1.pop(i)
        e2.pop(i)
        epoch.pop(i)
    else:
        i+=1

##Secondary velocity is actually primary and vice versa
#vn1=zeros(len(v1))
#vn2=zeros(len(v2))
#en1=zeros(len(e1))
#en2=zeros(len(e2))
#for i in range(len(v1)):
#    vn1[i]=v2[i]
#    vn2[i]=v1[i]
#    en1[i]=e2[i]
#    en2[i]=e1[i]
#v1=vn1
#v2=vn2
#e1=en1
#e2=en2


def chi(y1,y2,D1,D2,err1,err2):
    """
    Calculate chi^2 given predicted values (y), observed values (D), and errors
    on observed values (err)
    """
    total=0.
    for i in range(len(D1)):
        total+=(y1[i]-D1[i])**2/err1[i]**2
        total+=(y2[i]-D2[i])**2/err2[i]**2
    return total

def conf(arr):
    s=sort(arr)
    n=len(s)
    i=int(.34*len(s))
    if n%2==1:
        minus=median(s)-s[(n-1)/2-i]
        plus=s[(n-1)/2+i]-median(s)
    else:
        minus=median(s)-s[n/2-1-i]
        plus=s[n/2+i]-median(s)
    minus=float('%s' % float('%.2g' % minus))
    plus=float('%s' % float('%.2g' % plus))
    return plus,minus

#Create chain arrays for each parameter with nw walkers and nj max jumps
#Fit parameters:
#a - semi-major axis of orbit,
#P - period
#e - eccentricity
#w - little omega - argument of periastron
#g - system velocity
#T - time of periastron
#q - mass ratio
nw=4
nj=5*10**5
params={'a':zeros([nw,nj+1],dtype=double),
        'P':zeros([nw,nj+1],dtype=double),
        'e':zeros([nw,nj+1],dtype=double),
        'w':zeros([nw,nj+1],dtype=double),
        'g':zeros([nw,nj+1],dtype=double),
        'T':zeros([nw,nj+1],dtype=double),
        'q':zeros([nw,nj+1],dtype=double)}
keys=params.keys()

#Use RVs from the literature?
use_lit_rvs=True
lit_file=genfromtxt('CUCnc_lit_RVs.txt',dtype=float)
if use_lit_rvs:
    print 'Including RVs from the literature'
    for i in range(len(lit_file)):
        #Convert literature RV epochs to HJD
        hjd=helio_jd(lit_file[i,0]-2400000,128.,19.39)+2400000
        epoch.append(hjd)
        v1.append(lit_file[i,1])
        v2.append(lit_file[i,3])
        e1.append(.254)
        e2.append(.482)

#Derived parameter chains
Msum=zeros([nw,nj+1],dtype=double)
M1=zeros([nw,nj+1],dtype=double)
M2=zeros([nw,nj+1],dtype=double)
chi2=zeros([nw,nj+1],dtype=double)

#Claculate initial guess for g and q using plot of v1 versus v2
coeffs=polyfit(v2,v1,1)

#Use previous timing/period measurements as priors?
use_eclipse_prior=True
num_eclipses=1
eclipses=zeros([num_eclipses,2])
eclipses[0]=[2450208.5068,0.00079]

use_period_prior=True
period_prior=[2.771468,0.000004]

#initial guesses for  fit parameters 
params['a'][:,0]=double(0.03628) #AU
params['P'][:,0]=double(2.771468) #days
params['e'][:,0]=double(.000)
params['w'][:,0]=double(270*pi/180.) #radians
params['q'][:,0]=double(-coeffs[0]) 
params['g'][:,0]=double(coeffs[1]/(1+params['q'][0,0])) #km/s
params['T'][:,0]=double(2450208.5068) #HJD

#Set jump sizes
jumps={'T':double(0.00079),
        'P':double(0.0000002),
        'a':double(0.000021),
        'e':double(0.004),
        'w':double(0.0048),
        'g':double(0.075),
        'q':double(0.001)}
        
#Set circular orbit flag?
circle=True
if circle:
    print 'Circular flag set'
    params['e'][:,0]=double(0.0)
    params['w'][:,0]=double(270*pi/180.)
    keys.pop(keys.index('e'))
    keys.pop(keys.index('w'))

#Calculate initial derived parameters
Msum[:,0]=params['a'][:,0]**3/(params['P'][:,0]/365.242)**2 #Solar masses
M1[:,0]=Msum[:,0]/(1+params['q'][:,0]) #Solar masses
M2[:,0]=Msum[:,0]*params['q'][:,0]/(1+params['q'][:,0]) #Solar masses
y1,y2=calc_velocity(params['a'][0,0],params['P'][0,0],params['e'][0,0], \
    params['g'][0,0],params['w'][0,0],params['q'][0,0],params['T'][0,0],epoch)
chi2[:,0]=chi(y1,y2,v1,v2,e1,e2)

if use_period_prior:
    print 'Using period prior'
    chi2[:,0]+=((params['P'][0,0]-period_prior[0])/period_prior[1])**2
if use_eclipse_prior:
    print 'Using eclipse prior'
    if tan(params['w'][0,0])>10000:
        t0=params['T'][0,0]
    else:
        #Calculate time of primary eclipse using e,w,P,T
        theta=2*arctan((tan(params['w'][0,0])+sqrt(tan(params['w'][0,0])**2-\
            params['e'][0,0]**2+1))/(params['e'][0,0]-1))
        E=2*arctan(sqrt((1-params['e'][0,0])/(1+params['e'][0,0]))*tan(theta/2))
        M=E-params['e'][0,0]*sin(E)
        t0=params['P'][0,0]*M/(2*pi)+params['T'][0,0]
    for i in range(len(eclipses)):
        chi2[:,0]+=(abs(t0-eclipses[i,0])%params['P'][0,0])**2\
            /eclipses[i,1]**2

#Store number of accepted and rejected jumps for each parameter
accepted={'T':0.,'P':0.,'a':0.,'e':0.,'w':0.,'g':0.,'q':0.}
rejected={'T':0.,'P':0.,'a':0.,'e':0.,'w':0.,'g':0.,'q':0.}
rates={}

from random import *
#Start jumping
print 'Starting chains'
step=1
for k in range(nj/100000):    
    for z in range(100000): 
        for i in range(4): #Iterate over 4 walkers
            #Randomly select parameter to change
            param=choice(keys)
            for key in params.keys():
                if key==param:
                    #Assign new Gaussian random value to param
                    params[key][i,step]=gauss(params[key][i,step-1],jumps[key])
#                    if param=='w' and params[key][i,step]<0:
#                        params['w'][i,step]+=2*pi
                    if param=='w' and params[key][i,step]>2*pi:
                        params['w'][i,step]-=2*pi
                    elif param=='e' and params[key][i,step]<0:
                        params['w'][i,step]+=pi
                        params['T'][i,step]+=params['P'][i,step]/2.
                        params['e'][i,step]=abs(params['e'][i,step])
                else:
                    #Keep all other parameters constant
                    params[key][i,step]=params[key][i,step-1]
            #Calculate new derived parameters
            Msum[i,step]=params['a'][i,step]**3/(params['P'][i,step]/365.242)**2 
            M1[i,step]=Msum[i,step]/(1+params['q'][i,step]) #Solar masses
            M2[i,step]=Msum[i,step]*params['q'][i,step]/(1+params['q'][i,step]) 
            y1,y2=calc_velocity(params['a'][i,step],params['P'][i,step],\
                params['e'][i,step],params['g'][i,step],params['w'][i,step],\
                params['q'][i,step],params['T'][i,step],epoch)
            chi2[i,step]=chi(y1,y2,v1,v2,e1,e2)
            if use_period_prior:
                chi2[i,step]+=((params['P'][i,step]-period_prior[0])\
                    /period_prior[1])**2
            if use_eclipse_prior:
                if tan(params['w'][i,step])>10000:
                    t0=params['T'][i,step]
                else:
                    #Calculate time of primary eclipse using e,w,P,T
                    theta=2*arctan((tan(params['w'][i,step])+\
                        sqrt(tan(params['w'][i,step])**2-params['e'][i,step]**2\
                        +1))/(params['e'][i,step]-1))
                    E=2*arctan(sqrt((1-params['e'][i,step])\
                        /(1+params['e'][i,step]))*tan(theta/2))
                    M=E-params['e'][i,step]*sin(E)
                    t0=params['P'][i,step]*M/(2*pi)+params['T'][i,step]
                for j in range(len(eclipses)):
                    chi2[i,step]+=(abs(t0-eclipses[j,0])%\
                        params['P'][i,step])**2/eclipses[j,1]**2
            #Calculate likelihood ratio
            ratio=exp(.5*(chi2[i,step-1]-chi2[i,step]))
            #Keep all jumps that lower chi2 (ratio>1) and keep some jumps that 
            #increase chi2 (if ratio>random deviate from [0,1))
            deviate=random()
##            if param=='T':
##                print params['T'][i,step],chi2[i,step],ratio,deviate
            if ratio>=deviate:
                accepted[param]+=1
            else:
                rejected[param]+=1
                #Revert back to previous values
                for key in params:
                    params[key][i,step]=params[key][i,step-1]
                Msum[i,step]=Msum[i,step-1]
                M1[i,step]=M1[i,step-1]
                M2[i,step]=M2[i,step-1]
                chi2[i,step]=chi2[i,step-1]
        step+=1
    print step-1
    
    #Adjust jump sizes to be close to 32%
    for key in keys:
        rates[key]=accepted[key]/(rejected[key]+accepted[key])
        if rates[key]<.28 or rates[key]>.38:
            jumps[key]*=rates[key]/.32
            print '%s jump size changed by a factor of %f' % (key,\
                rates[key]/.32)
        
    #Check for convergence using 1-sigma confidence interval
    converge=1
    for key in keys:
        m=vstack(([median(params[key][0,k*10**5:(k+1)*10**5]),\
                median(params[key][1,k*10**5:(k+1)*10**5]),
                median(params[key][2,k*10**5:(k+1)*10**5]),\
                median(params[key][3,k*10**5:(k+1)*10**5])]))
        plus,minus=conf(params[key][0,k*10**5:(k+1)*10**5])
        for b in range(1,4):
            if m[b]>m[0]+plus or m[b]<m[0]-minus:
                print '%s hasn\'t converged yet' % (key)
                converge=0
                break
#    if converge==1 and k>0:
#        break            
    
#Remove burn-in of 100,000 steps
for key in params:
    params[key]=params[key][:,100000:step]
Msum=Msum[:,100000:step]
M1=M1[:,100000:step]
M2=M2[:,100000:step]
chi2=chi2[:,100000:step]
##for i in range(4):
##    print 'M1 = %f +%f/-%f' % (median(M1[i]),conf(M1[i])[0],conf(M1[i])[1])
##    print 'M2 = %f +%f/-%f' % (media(M2[i]),conf(M2[i])[0],conf(M2[i])[1])
##    for key in keys:
##        print '%s = %f +%f/-%f' % (key,median(params[key][i]),\
##            conf(params[key][i])[0],conf(params[key][i])[1])
##    print '\n'

confidence={}
print 'Jump Acceptance Rates:'
#Combine chains from different walkers
M1=list(M1[0])+list(M1[1])+list(M1[2])+list(M1[3])
M2=list(M2[0])+list(M2[1])+list(M2[2])+list(M2[3])
chi2=list(chi2[0])+list(chi2[1])+list(chi2[2])+list(chi2[3])
for key in keys:
    rates[key]=accepted[key]/(rejected[key]+accepted[key])
    print key,rates[key]
    if key=='w':
        params[key]*=180./pi
    params[key]=list(params[key][0])+list(params[key][1])+list(params[key][2])\
    +list(params[key][3])
    confidence[key]=conf(params[key])

f=open('MCMC_results.txt','a')
print >>f, '\nSystem nj circle chi2 M1 M1_err M2 M2_err'
print >>f, 'CU Cnc %i %s %f %f +%f/-%f %f +%f/-%f' % (nj,circle,\
    median(chi2),median(M1),conf(M1)[0],conf(M1)[1],median(M2),conf(M2)[0],\
    conf(M2)[1])
print >>f, 'Period Prior:'
print >>f, period_prior
print >>f, 'Eclipse Prior:'
print >>f, eclipses
print >>f, 'Initial guesses a,P,e,w,q,g,T'
if circle:
    rates['e']='n/a'
    rates['w']='n/a'
    print >>f, params['a'][0],params['P'][0],0.0,1.5*pi,\
      params['q'][0],params['g'][0],params['T'][0]
else:
    print >>f, params['a'][0],params['P'][0],params['e'][0],params['w'][0],\
      params['q'][0],params['g'][0],params['T'][0]
print >>f, 'Jumps sizes a,P,e,w,q,g,T'
print >>f, jumps['a'],jumps['P'],jumps['e'],jumps['w'],jumps['q'],jumps['g'],\
      jumps['T']
print >>f, 'Jump acceptance rates a,P,e,w,q,g,T'
print >>f, rates['a'],rates['P'],rates['e'],rates['w'],rates['q'],rates['g'],\
      rates['T']
for key in keys:
    print >>f, '%s = %.8f +%.8f/-%.8f' % (key,median(params[key]),confidence[key][0],\
        confidence[key][1])
f.close()

print 'Final parameter values:'
for key in keys:
    print '%s = %.8f +%.8f/-%.8f' % (key,median(params[key]),\
        confidence[key][0],confidence[key][1])
print 'chi2 = %f +%f/-%f' % (median(chi2),conf(chi2)[0],conf(chi2)[1])
print 'M1 = %f +%f/-%f' % (median(M1),conf(M1)[0],conf(M1)[1])
print 'M2 = %f +%f/-%f' % (median(M2),conf(M2)[0],conf(M2)[1])

#Make RV plot
if circle:
    e=0.0
    w=1.5*pi
else:
    e=median(params['e'])
    w=median(params['w'])*pi/180.
a=median(params['a'])
g=median(params['g'])
q=median(params['q'])
P=median(params['P'])
T=median(params['T'])

keys=sort(rv.keys())
for key in keys:
    if abs(rv[key]['m'][1]-rv[key]['m'][0])>mvs:
        phase=((float(key)-T)%P)/P
        print 'CU Cnc & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\' % \
          (key[3:13],phase,rv[key]['m'][0],rv[key]['e'][0],rv[key]['m'][1],\
           rv[key]['e'][1])

###Calculate time of primary eclipse using e,w,P,T
##theta=2*arctan((tan(w)+sqrt(tan(w)**2-e**2+1))/(e-1))
##E=2*arctan(sqrt((1-e)/(1+e))*tan(theta/2))
##M=E-e*sin(E)
##t0=P*M/(2*pi)+T
        
phase=zeros(len(epoch))
for i in range(len(epoch)):
    phase[i]=((epoch[i]-T)%P)/P

y1,y2=calc_velocity(a,P,e,g,w,q,T,epoch)
resid1=zeros(len(epoch))
resid2=zeros(len(epoch))
for i in range(len(epoch)):
    resid1[i]=v1[i]-y1[i]
    resid2[i]=v2[i]-y2[i]
time=linspace(T,T+P,10000)
fit1,fit2=calc_velocity(a,P,e,g,w,q,T,time)
phasesec=linspace(0,1,10000)
ax1=subplot(211)
ax1.scatter(phase,v1,c='blue')
ax1.scatter(phase,v2,c='red')
ax1.plot(phasesec,fit1,'blue')
ax1.plot(phasesec,fit2,'red')
xlim(0,1)
ylabel('Velocity (km/s)')
title('CU Cnc')
ax2=subplot(212)
ax2.errorbar(phase,resid1,yerr=e1,fmt='none',ecolor='blue')
ax2.errorbar(phase,resid2,yerr=e2,fmt='none',ecolor='red')
if use_lit_rvs:
    #Mark my RV points with larger green points
    for i in range(len(epoch)):
        for j in range(len(rv.keys())):
            if abs(epoch[i]-float(rv.keys()[j]))<.0001:
                ax1.scatter(phase[i],v1[i],c='green',s=40)
                ax1.scatter(phase[i],v2[i],c='green',s=40)
                ax2.errorbar(phase[i],resid1[i],yerr=e1[i],fmt='none',ecolor=\
                    'green')
                ax2.errorbar(phase[i],resid2[i],yerr=e2[i],fmt='none',ecolor=\
                    'green')
ax2.hlines(0,0,1,linestyles='dashed')
xlim(0,1)
xlabel('Phase')
ylabel('O-E(km/s)')
pos=ax1.get_position().bounds
newpos=[pos[0],pos[1]-.2,pos[2],pos[3]+.2]
ax1.set_position(newpos)
pos=ax2.get_position().bounds
newpos=[pos[0],pos[1],pos[2],pos[3]-.2]
ax2.set_position(newpos)
savefig('CUCnc_circle.png')
close()

#Make triangle plot
am=array(matrix(params['a']).T)
Pm=array(matrix(params['P']).T)
gm=array(matrix(params['g']).T)
qm=array(matrix(params['q']).T)
Tm=array(matrix(params['T']).T)
Tm-=int(T)
#Pm-=int(P)
if circle:
    samples=hstack(([am,Pm,gm,Tm,qm]))
    figure=corner(samples,labels=[r'$a$ (AU)',r'$P$ (days)',\
        r'$\gamma$ (km/s)',r'$T_0$ (HJD-%i)'%int(T),r'$q$'])
else:
    em=array(matrix(params['e']).T)
    wm=array(matrix(params['w']).T)
    samples=hstack(([am,Pm,em,wm,gm,Tm,qm]))
    figure=corner(samples,labels=[r'$a$ (AU)',r'$P$ (days)',\
        r'$e$',r'$\omega$ (deg)',r'$\gamma$ (km/s)',r'$T_0$ (HJD-%i)'%int(T),\
        r'$q$'])
savefig('CUCnc_circle_triangle.png')
close()
