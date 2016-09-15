from numpy import *
from PyAstronomy.pyasl import *


stars=genfromtxt('EB_Targets.csv',delimiter=',',dtype=list,skiprows=1)
t0=zeros([len(stars),2])
P=zeros([len(stars),2])
coords=zeros([len(stars),2])

for i in range(len(stars)):
    ra=stars[i,4].split(' ')
    dec=stars[i,5].split(' ')
    ra=(float(ra[0])+float(ra[1])/60+float(ra[2])/3600)*15 #Convert to deg
    if int(dec[0])>=0:
        dec=float(dec[0])+float(dec[1])/60+float(dec[2])/3600
    elif int(dec[0])<0:
        dec=float(dec[0])-float(dec[1])/60-float(dec[2])/3600
    coords[i]=[ra,dec]
    
    if stars[i,10]!='': #Check for time of minimum
        if 'MJD' in stars[i,10]:
            t0[i,0]=float(stars[i,10][3:])+0.5 #Convert to reduced JD
            t0[i,0]=helio_jd(t0[i,0],ra,dec)
        elif 'HJD' in stars[i,10]:
            t0[i,0]=float(stars[i,10][3:])-2400000 #HJD-240 0000
        elif 'RJD' in stars[i,10]:
            t0[i,0]=float(stars[i,10][3:])-2400000 #Convert to reduced JD
            t0[i,0]=helio_jd(t0[i,0],ra,dec)
        P[i,0]=float(stars[i,8])
        try:
            t0[i,1]=float(stars[i,11])
            P[i,1]=float(stars[i,9])
        except ValueError:
            t0[i,1]=.014
            P[i,1]=1.

time=linspace(57138.5,57150.125,559)
phases=zeros([len(time),len(stars)])
errors=zeros(len(stars))
for i in range(len(stars)):
    if stars[i,10]!='':
        for j in range(len(time)):
            t=helio_jd(time[j],coords[i,0],coords[i,1])
            phases[j,i]=(t-t0[i,0])/P[i,0]-int((t-t0[i,0])/P[i,0])
        errors[i]=sqrt(t0[i,1]**2/P[i,0]**2+(57074.5-t0[i,0])**2/P[i,0]**4* \
                       P[i,1]**2)
phases=vstack(([stars[:,0],stars[:,4],stars[:,5],stars[:,1],stars[:,8],errors \
                ,phases]))

savetxt('April_Phases_2.csv',phases,fmt='%s',delimiter=',')

    
