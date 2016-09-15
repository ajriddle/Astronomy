from numpy import *
from matplotlib.pyplot import *
import scipy.stats as stat
from EB_Velocity import *
from scipy.optimize import curve_fit
from EB_Velocity import *

def func(t,*a):
    A,f,p,const=a
    return A*sin(2*pi*f*t+p)+const
    
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
    
d=genfromtxt('current_rvs',dtype=list)
rv={}
epoch=[]
v1=[]
v2=[]
e1=[]
e2=[]

for i in range(len(d)):
    if 'M' in d[i,1] or 'E' in d[i,0]:
        key=d[i,2]
        rv[key]={'v':zeros(5),'m':[0.,0.],'e':[0.,0.],'bo':[]}
    else:
        for j in range(5):
            d[i,j]=float(d[i,j])
        rv[key]['v']=vstack(([rv[key]['v'],d[i]]))
        
#Do two iterations of sigma clipping  
for key in rv:
    rv[key]['v']=rv[key]['v'][1:]
    for j in range(2):
        #errorbar(rv[key]['v'][:,0],rv[key]['v'][:,1],rv[key]['v'][:,2],fmt='o'\
         #   ,c='b',label='Primary')
        #errorbar(rv[key]['v'][:,0],rv[key]['v'][:,3],rv[key]['v'][:,4],fmt='o'\
         #   ,c='r',label='Secondary')
        i=0 
        while i<len(rv[key]['v']):
            a1=stat.nanmedian(list(rv[key]['v'][:,1]))
            a2=stat.nanmedian(list(rv[key]['v'][:,3]))
            if key=='54467.371626':
                a1=81.57
                a2=30.73
            std1=std(rv[key]['v'][:,1])
            std2=std(rv[key]['v'][:,3])
            if (abs(rv[key]['v'][i,3]-a2)/std2>2 or abs(rv[key]['v'][i,3]-a2)>30)\
               and abs(rv[key]['v'][i,3]-a2)>3:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))              
            elif (abs(rv[key]['v'][i,1]-a1)/std1>2 or  \
            abs(rv[key]['v'][i,1]-a1)>30) and abs(rv[key]['v'][i,1]-a1)>3:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))
            elif abs(rv[key]['v'][i,0]-6553)<10 or abs(rv[key]['v'][i,0]-8229) \
                <10:
                rv[key]['bo'].append(rv[key]['v'][i,0])
                rv[key]['v']=vstack((rv[key]['v'][:i],rv[key]['v'][i+1:]))
            else:
                i+=1
        #xlabel('Order Wavelength (A)')
        #ylabel('Velocity (km/s)')
        #xmin,xmax=xlim()
        #hlines(a1,xmin,xmax,color='b',linestyles='solid')
        #hlines(a1+2*std1,xmin,xmax,color='b',linestyles='dashed')
        #hlines(a1-2*std1,xmin,xmax,color='b',linestyles='dashed')
        #hlines(a2,xmin,xmax,color='r',linestyles='solid')
        #hlines(a2+2*std2,xmin,xmax,color='r',linestyles='dashed')
        #hlines(a2-2*std2,xmin,xmax,color='r',linestyles='dashed')
        #legend(loc='best')
        #title(r'$v_1=%f \pm %f,v_2=%f \pm %f$' % (a1,std1,a2,std2))
        #savefig('%s_%i.pdf' % (key,j))
        #close()
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
    
    e1.append(max((s1/sum(w1)/(len(w1)-1))**.5,1))
    e2.append(max((s2/sum(w2)/(len(w2)-1))**.5,1))
    if len(rv[key]['v'])==1:
        e1[-1]=rv[key]['v'][0][2]
        e2[-1]=rv[key]['v'][0][4]
    rv[key]['m']=[v1[-1],v2[-1]]
    rv[key]['e']=[e1[-1],e2[-1]]
    if key=='54467.371626':
        v1[-1],v2[-1]=v2[-1],v1[-1]
        e1[-1],e2[-1]=e2[-1],e1[-1]
    if key=='54690.482060':
        v2[-1]=-56.26

##v1=array(matrix(v1).T)
##e1=array(matrix(e1).T)
##v2=array(matrix(v2).T)
##e2=array(matrix(e2).T)
##epoch=array(matrix(epoch).T)
##data=hstack(([epoch,v1,e1,v2,e2]))
##savetxt('RV_Data.txt',data)

#Guess fit parameters
#v1[0],v2[0]=v2[0],v1[0]
#e1[0],e2[0]=e2[0],e1[0]
coeffs=polyfit(v2,v1,1)
q=-coeffs[0]
g=coeffs[1]/(1+q)
P=4.46818
e=0.0
a=0.051
T=2454285.881
w=1.5*pi
time=linspace(min(epoch)-3,max(epoch)+3,100000)
y1,y2=calc_velocity(a,P,e,g,w,q,T,time)
errorbar(epoch,v1,yerr=e1,fmt='o',ecolor='blue',label='Primary')
errorbar(epoch,v2,yerr=e2,fmt='o',ecolor='green',label='Seconday')
plot(time,y1,'b')
plot(time,y2,'g')
ylabel('Velocity (km/s)')
xlabel('HJD')
legend(loc='best')
show()

#phase=zeros(len(epoch))
#for i in range(len(epoch)):
#    phase[i]=((epoch[i]-T)%P)/P    
#y1,y2=calc_velocity(a,P,e,g,w,q,T,epoch)
#chi2=chi(y1,y2,v1,v2,e1,e2)
#print chi2
#resid1=zeros(len(epoch))
#resid2=zeros(len(epoch))
#for i in range(len(epoch)):
#    resid1[i]=v1[i]-y1[i]
#    resid2[i]=v2[i]-y2[i]
#time=linspace(T,T+P,10000)
#fit1,fit2=calc_velocity(a,P,e,g,w,q,T,time)
#phasesec=linspace(0,1,10000)
#ax1=subplot(211)
#ax1.scatter(phase,v1,c='blue')
#ax1.scatter(phase,v2,c='red')
#ax1.plot(phasesec,fit1,'blue')
#ax1.plot(phasesec,fit2,'red')
#xlim(0,1)
#ylabel('Velocity (km/s)')
#ax2=subplot(212)
#ax2.errorbar(phase,resid1,yerr=e1,fmt='none',ecolor='blue')
#ax2.errorbar(phase,resid2,yerr=e2,fmt='none',ecolor='red')
#ax2.hlines(0,0,1,linestyles='dashed')
#xlim(0,1)
#xlabel('Phase')
#ylabel('O-E(km/s)')
#pos=ax1.get_position().bounds
#newpos=[pos[0],pos[1]-.2,pos[2],pos[3]+.2]
#ax1.set_position(newpos)
#pos=ax2.get_position().bounds
#newpos=[pos[0],pos[1],pos[2],pos[3]-.2]
#ax2.set_position(newpos)
#show()

#p0=[80.,0.67,0.,-50.]
#coeff,var_matrix=curve_fit(func,epoch,v1,p0=p0)
##for i in range(len(epoch)):
##    print epoch[i],v1[i],v2[i]

##a=zeros([3,1])
##a[0]=50.
##a[1]=50.
##a[2]=0.
##v=v2
#while (abs(delta[0])>0.05 or abs(delta[1])>0.5 or abs(delta[2])>1E-4 or  \
 #     abs(delta[3])>5) and iteration<100:

##P=2.91334166666667
##f=2*pi/P
##
##t0=4386.374
##t0+=50000
###theta=2*pi-mod(f*t0,2*pi)
##theta=-2*pi*t0/P
##
##phase=zeros(len(epoch))
##for i in range(len(epoch)):
##    #phase[i]=mod(f*epoch[i]+theta,2*pi)/(2*pi)
##    phase[i]=mod(epoch[i]-t0,P)/P

##vsys=-6.2
##a=[]
##for component in [[v1,e1],[v2,e2]]:
##    Y=[0.]
##    N=[0.]
##    S=0.
##    v=component[0]
##    err=component[1]
##    for i in range(len(v)):
##        v[i]-=vsys
##        w=1/err[i]**2
##        N+=w*sin(f*epoch[i]+theta)**2
##        Y+=w*sin(f*epoch[i]+theta)*v[i]
##    N=matrix(N)
##    Y=matrix(Y)
##    C=N.I
##    a.append(float(C*Y))

##Y=zeros([2,1])
##N=zeros([2,2])
##for i in range(len(v1)):
##    w1=1/e1[i]**2
##    #w2=1/e2[i]**2
##    N[0,0]+=w1*sin(f*epoch[i]+theta)**2
##    N[0,1]+=w1*sin(f*epoch[i]+theta)
##    N[1,1]+=w1
##    Y[0]+=w1*sin(f*epoch[i]+theta)*v1[i]
##    Y[1]+=w1*v1[i]
##N[1,0]=N[0,1]
##N=matrix(N)
##Y=matrix(Y)
##C=N.I
##a=array(C*Y)
##S=0.
##for i in range(len(v1)):
##    w1=1/e1[i]**2
##    #w2=1/e2[i]**2
##    S+=w1*(v1[i]-a[0]*sin(f*epoch[i]+theta)-a[1])**2
##C=S/(len(v1)-2-1)*array(C)

####Y=zeros([3,1])
####N=zeros([3,3])
####for i in range(len(v1)):
####    w1=1/e1[i]**2
####    w2=1/e2[i]**2
####    N[0,0]+=w1*sin(f*epoch[i]+theta)**2
####    N[0,1]+=w1*sin(f*epoch[i]+theta)
####    N[1,1]+=w1+w2
####    N[1,2]+=w2*sin(f*epoch[i]+theta)
####    N[2,2]+=w2*sin(f*epoch[i]+theta)**2
####    Y[0]+=w1*sin(f*epoch[i]+theta)*v1[i]
####    Y[1]+=w1*v1[i]+w2*v2[i]
####    Y[2]+=w2*sin(f*epoch[i]+theta)*v2[i]
####N[1,0]=N[0,1]
####N[2,1]=N[1,2]
####N=matrix(N)
####Y=matrix(Y)
####C=N.I
####a=array(C*Y)
####S=0.
####for i in range(len(v1)):
####    w1=1/e1[i]**2
####    w2=1/e2[i]**2
####    S+=w1*(v1[i]-a[0]*sin(f*epoch[i]+theta)-a[1])**2+ \
####        w2*(v2[i]-a[2]*sin(f*epoch[i]+theta)-a[1])**2
####C=S/(2*len(v1)-3-1)*array(C)

##Y=zeros([5,1])
##N=zeros([5,5])
##for i in range(len(v1)):
##    w1=1/e1[i]**2
##    w2=1/e2[i]**2
##    v1[i]-=vsys
##    v2[i]-=vsys
##    N[0,0]+=w1*sin(f*epoch[i])**2
##    N[0,1]+=w1*sin(f*epoch[i])*cos(f*epoch[i])
##    N[0,2]+=w1*sin(f*epoch[i])
##    N[1,1]+=w1*cos(f*epoch[i])**2
##    N[1,2]+=w1*cos(f*epoch[i])
##    N[2,2]+=w1+w2
##    N[2,3]+=w2*sin(f*epoch[i])
##    N[2,4]+=w2*cos(f*epoch[i])
##    N[3,3]+=w2*sin(f*epoch[i])**2
##    N[3,4]+=w2*sin(f*epoch[i])*cos(f*epoch[i])
##    N[4,4]+=w2*cos(f*epoch[i])**2
##    Y[0]+=w1*sin(f*epoch[i])*v1[i]
##    Y[1]+=w1*cos(f*epoch[i])*v1[i]
##    Y[2]+=w1*v1[i]+w2*v2[i]
##    Y[3]+=w2*sin(f*epoch[i])*v2[i]
##    Y[4]+=w2*cos(f*epoch[i])*v2[i]
##N[1,0]=N[0,1]
##N[2,0]=N[0,2]
##N[2,1]=N[1,2]
##N[3,2]=N[2,3]
##N[4,2]=N[2,4]
##N[4,3]=N[3,4]
##N=matrix(N)
##Y=matrix(Y)
##C=N.I
##a=array(C*Y)
##S=0.
##for i in range(len(v1)):
##    w1=1/e1[i]**2
##    w2=1/e2[i]**2
##    S+=w1*(v1[i]-a[0]*sin(f*epoch[i])-a[1]*cos(f*epoch[i])-a[2])**2+ \
##        w2*(v2[i]-a[3]*sin(f*epoch[i])-a[4]*cos(f*epoch[i])-a[2])**2
##C=S/(2*len(v1)-5-1)*array(C)
    
##a=zeros([2,2])
##z=0
##for component in [[v1,e1],[v2,e2]]:
##    Y=zeros([2,1])
##    N=zeros([2,2])
##    S=0.
##    v=component[0]
##    err=component[1]
##    for i in range(len(v)):
##        v[i]-=vsys
##        w=1/err[i]**2
##        #S+=w*(v[i]-a[0]*sin(f*epoch[i])-a[1]*cos(f*epoch[i])-a[2])**2
##        N[0,0]+=w*sin(f*epoch[i])**2
##        N[0,1]+=w*sin(f*epoch[i])*cos(f*epoch[i])
##        N[1,1]+=w*cos(f*epoch[i])**2  
##        Y[0]+=w*sin(f*epoch[i])*v[i]
##        Y[1]+=w*cos(f*epoch[i])*v[i]
##    N[1,0]=N[0,1]
##    N=matrix(N)
##    Y=matrix(Y)
##    C=N.I
##    a[z]=array((C*Y).T)    
##    A=sqrt(a[z,0]**2+a[z,1]**2)
##    phi=arctan2(a[z,1],a[z,0])
##    for i in range(len(v)):
##        w=1/err[i]**2
##        S+=w*(v[i]-a[z,0]*sin(f*epoch[i])-a[z,1]*cos(f*epoch[i]))**2
##    C=S/(len(v)-2-1)*array(C)
##    #print 'Period = %f +/- %f days' % (2*pi/a[2],2*pi/a[2]**2*C[2,2]**.5)
##    print 'K = %f +/- %f km/s' % (sqrt(a[z,0]**2+a[z,1]**2),((2*a[z,0]*a[z,1] \
##        *C[0,1]+a[z,0]**2*C[0,0]+a[z,1]**2*C[1,1])/(a[z,0]**2+a[z,1]**2))**.5)
##    print 'angular shift = %f' % (phi)
##    #print 'Velocity offset = %f +/- %f km/s' % (a[z,2],sqrt(C[2,2]))
##    print 'Reduced Chi^2 = %f' % (S/(len(v)-2-1))
####    if skip=='n':
####        print >>g, sqrt(a[z,0]**2+a[z,1]**2),((2*a[z,0]*a[z,1]*C[0,1]+a[z,0]**2\
####            *C[0,0]+a[z,1]**2*C[1,1])/(a[z,0]**2+a[z,1]**2))**.5,a[z,2], \
####            sqrt(C[2,2]),arctan2(a[0,1],a[0,0]),sqrt((1+(a[z,1]/a[z,0])**2)**-2\
####            *a[z,0]**-2*((a[z,1]/a[z,0])**2*C[0,0]+C[1,1]-a[z,1]/a[z,0]* \
####            C[0,1])),S/(len(epoch)-3-1),f
##    z+=1
##phi2=arctan2(a[1,1],a[1,0])
##phi1=arctan2(a[0,1],a[0,0])
##phase1=zeros(len(v1))
##phase2=zeros(len(v2))
##for i in range(len(v1)):
##    phase1[i]=mod(f*epoch[i]+phi1,2*pi)/(2*pi)
##    phase2[i]=mod(f*epoch[i]+phi2,2*pi)/(2*pi)
##
##A1=sqrt(a[0,0]**2+a[0,1]**2)
##A2=sqrt(a[1,0]**2+a[1,1]**2)
##x=linspace(0,2*pi,10000)
##fit1=[-A1*sin(t) for t in x]
##fit2=[A2*sin(t) for t in x]
##x=linspace(0,1,10000)
##ax1=subplot(211)
##ax1.plot(x,fit1,'k--')
###ax1.plot(x,fit2,'k--')
###ax1.scatter(phase,v2,c='r')
##ax1.scatter(phase,v1,c='b')
##resid1=[float(-A1*sin(f*t+phi)) for t in epoch]
##resid2=[float(A2*sin(f*t+phi)) for t in epoch]
##for i in range(len(v1)):
##    resid2[i]=v2[i]-resid2[i]
##    resid1[i]=v1[i]-resid1[i]
##ylabel('Velocity (km/s)')
##xlim(0,1)
###ylim(ymin=-100)
##title(d[0,1])
##
##ax2=subplot(212)
###ax2.errorbar(phase+.005,resid2,yerr=e2,fmt='none',ecolor='red')
##ax2.errorbar(phase,resid1,yerr=e1,fmt='none',ecolor='blue')
##ax2.hlines(0,0,1,linestyles='dashed')
##xlim(0,1)
###ylim(-2,1)
##xlabel('Phase')
##ylabel('Residual (km/s)')
##pos=ax1.get_position().bounds
##newpos=[pos[0],pos[1]-.2,pos[2],pos[3]+.2]
##ax1.set_position(newpos)
##pos=ax2.get_position().bounds
##newpos=[pos[0],pos[1],pos[2],pos[3]-.2]
##ax2.set_position(newpos)
##savefig('%s_Velocity_Fit.pdf' % d[0,1])
##close()
###show()

    
#Fit using same system offset and same phase shift   
##A1=sqrt(a[0]**2+a[1]**2)
##phi1=arctan2(a[1],a[0])
##A2=sqrt(a[3]**2+a[4]**2)
##phi2=arctan2(a[4],a[3])
##print 'K1 = %f +/- %f km/s' % (A1,((2*a[0]*a[1] \
##        *C[0,1]+a[0]**2*C[0,0]+a[1]**2*C[1,1])/(a[0]**2+a[1]**2))**.5)
##print 'Velocity offset = %f +/- %f km/s' % (a[2],sqrt(C[2,2]))
##print 'K2 = %f +/- %f km/s' % (A2,((2*a[3]*a[4] \
##        *C[3,4]+a[0]**2*C[3,3]+a[1]**2*C[4,3])/(a[0]**2+a[1]**2))**.5)
##print 'Reduced Chi^2 = %f' % (S/(2*len(v1)-5-1))

####x=linspace(0,2*pi,1000)
####fit=[a[2]*sin(t)+a[1] for t in x]
#####fit=[a[1]*sin(t) for t in x]
####x=linspace(0,1,1000)
####ax1=subplot(211)
####ax1.plot(x,fit,'k--')
####ax1.scatter(phase,v2,c='b')
####resid=[float(a[2]*sin(f*t+theta)+a[1]) for t in epoch]
#####resid=[float(a[1]*sin(f*t+theta)) for t in epoch]
####for i in range(len(resid)):
####    resid[i]=v2[i]-resid[i]
####x=linspace(0,2*pi,1000)
####fit=[a[0]*sin(t)+a[1] for t in x]
#####fit=[a[0]*sin(t) for t in x]
####x=linspace(0,1,1000)
####ax1.plot(x,fit,'k--')
####ax1.scatter(phase,v1,c='r')
####ylabel('Velocity (km/s)')
####xlim(0,1)
#####ylim(ymin=-100)
####title(d[0,1])
####
####ax2=subplot(212)
####ax2.errorbar(phase+.005,resid,yerr=e2,fmt='none',ecolor='blue')
####resid=[float(a[0]*sin(f*t+theta)+a[1]) for t in epoch]
#####resid=[float(a[0]*sin(f*t+theta)) for t in epoch]
####for i in range(len(resid)):
####    resid[i]=v1[i]-resid[i]
####ax2.errorbar(phase,resid,yerr=e1,fmt='none',ecolor='red')
####ax2.hlines(0,0,1,linestyles='dashed')
####xlim(0,1)
#####ylim(-2,1)
####xlabel('Phase')
####ylabel('Residual (km/s)')
####pos=ax1.get_position().bounds
####newpos=[pos[0],pos[1]-.2,pos[2],pos[3]+.2]
####ax1.set_position(newpos)
####pos=ax2.get_position().bounds
####newpos=[pos[0],pos[1],pos[2],pos[3]-.2]
####ax2.set_position(newpos)
####savefig('%s_Velocity_Fit.pdf' % (d[0,0]+d[0,1]))
####show()
####close()


####skip=raw_input('Skip rejected order save [y/n]? ')
####if skip=='n':
####    f=open('rejected_orders.txt','a')
####    for key in rv:
####        print >>f, d[0,1],key,len(rv[key]['v'])
####        print >>f, rv[key]['bo']
####    f.close()

####G=6.67E-11
####ka=abs(a[0])*1000
####kb=abs(a[2])*1000
####ea=sqrt(C[0,0])*1000
####eb=sqrt(C[2,2])*1000
####p=P*24*3600
####i=90*pi/180
####ma=(1+ka/kb)**-1*p/(2*pi*G)*(ka+kb)**3*sin(i)**-3
####mb=ma*(ka/kb)
####ema=sqrt((9*(ka+kb)**-2+(1+ka/kb)**-2/kb**2-6*(1+ka/kb)**-1*(ka+kb)**-1/kb)* \
####    ea**2+(9*(ka+kb)**-2+(1+ka/kb)**-2*ka**2/kb**4+6*(1+ka/kb)**-1*(ka+kb)**-1*\
####    ka/kb**2)*eb**2)
####ema*=ma
####emb=sqrt((ka/kb)**2*ema**2+(ma/kb)**2*ea**2+ma**2*(ka**2/kb**4)*eb**2)
####ma/=1.99E30
####mb/=1.99E30
####ema/=1.99E30
####emb/=1.99E30
####    
####print 'K1 = %f +/- %f km/s' % (-a[0],sqrt(C[0,0]))
####print 'K2 = %f +/- %f km/s' % (a[2],sqrt(C[2,2]))
####print 'Velocity offset = %f +/- %f km/s' % (a[1],sqrt(C[1,1]))
####print 'M1 = %f +/- %f solar masses' % (ma,ema)
####print 'M2 = %f +/- %f solar masses' % (mb,emb)
####print 'Reduced Chi^2 = %f' % (S/(2*len(v1)-3-1))#-1))
##
##skip=raw_input('Skip coefficient save [y/n]? ')
##if skip=='n':
##    g=open('fits.txt','a')
##    print >>g, d[0,1],-a[0],sqrt(C[0,0]),a[2],sqrt(C[2,2]),a[1],sqrt(C[1,1])#,\
##          #ma,ema,mb,emb
##    g.close()
