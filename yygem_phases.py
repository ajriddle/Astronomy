from numpy import *
from PyAstronomy.pyasl import *

t0=2449345.112327
P=0.814282212

ra='07 34 37.584'.split(' ')
dec='+31 52 11.05'.split(' ')
ra=(float(ra[0])+float(ra[1])/60+float(ra[2])/3600)*15 #Convert to deg
if int(dec[0])>=0:
    dec=float(dec[0])+float(dec[1])/60+float(dec[2])/3600
elif int(dec[0])<0:
    dec=float(dec[0])-float(dec[1])/60-float(dec[2])/3600
coords=[ra,dec]
    
#Check for time of minimum
t0-=2400000 #HJD-240 0000

time=linspace(57417.5,57423.5,289)
h=zeros([len(time),1])
h[0]=18.0
for i in range(1,len(h)):
    h[i]=h[i-1]+.5
    if h[i]>=24:
         h[i]-=24

phases=zeros([len(time),1])
for j in range(len(time)):
    t=helio_jd(time[j],coords[0],coords[1])
    phases[j]=(t-t0)/P-int((t-t0)/P)
a=hstack(([h,phases]))
savetxt('YYGem_Phases.csv',a,fmt='%s',delimiter=',')