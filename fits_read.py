from pyfits import *
from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import *

def rmspace(var):
    while '' in var:
        var.pop(var.index(''))
    return var

img1=open('HI.20101030.40490_1_Flux.fits')
data1=img1[0].data
img1.close()
img2=open('HI.20101030.40490_2_Flux.fits')
data2=img2[0].data
img2.close()
img3=open('HI.20101030.40490_3_Flux.fits')
data3=img3[0].data
img3.close()
img4=open('HI.20101030.40829_1_Flux.fits')
data4=img4[0].data
img4.close()
img5=open('HI.20101030.40829_2_Flux.fits')
data5=img5[0].data
img5.close()
img6=open('HI.20101030.40829_3_Flux.fits')
data6=img6[0].data
img6.close()

wv1=zeros([20,4056])
wv2=zeros([13,4056])
wv3=zeros([8,4056])
wv4=zeros([21,4056])
wv5=zeros([13,4056])
wv6=zeros([8,4056])
c1=zeros([20,7])
c2=zeros([13,7])
c3=zeros([8,7])
c4=zeros([21,7])
c5=zeros([13,7])
c6=zeros([8,7])

z=0
i=img1[0].header.index('WV_0_01')
while i<img1[0].header.index('WV_0_01')+40:
    s1=img1[0].header[i]
    s2=img1[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c1[z,0]=float(coeff1[0])
    c1[z,1]=float(coeff1[1])
    c1[z,2]=float(coeff1[2])
    c1[z,3]=float(coeff1[3])
    c1[z,4]=float(coeff2[0])
    c1[z,5]=float(coeff2[1])
    c1[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img2[0].header.index('WV_0_01')
while i<img2[0].header.index('WV_0_01')+26:
    s1=img2[0].header[i]
    s2=img2[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c2[z,0]=float(coeff1[0])
    c2[z,1]=float(coeff1[1])
    c2[z,2]=float(coeff1[2])
    c2[z,3]=float(coeff1[3])
    c2[z,4]=float(coeff2[0])
    c2[z,5]=float(coeff2[1])
    c2[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img3[0].header.index('WV_0_01')
while i<img3[0].header.index('WV_0_01')+16:
    s1=img3[0].header[i]
    s2=img3[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c3[z,0]=float(coeff1[0])
    c3[z,1]=float(coeff1[1])
    c3[z,2]=float(coeff1[2])
    c3[z,3]=float(coeff1[3])
    c3[z,4]=float(coeff2[0])
    c3[z,5]=float(coeff2[1])
    c3[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img4[0].header.index('WV_0_01')
while i<img4[0].header.index('WV_0_01')+42:
    s1=img4[0].header[i]
    s2=img4[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c4[z,0]=float(coeff1[0])
    c4[z,1]=float(coeff1[1])
    c4[z,2]=float(coeff1[2])
    c4[z,3]=float(coeff1[3])
    c4[z,4]=float(coeff2[0])
    c4[z,5]=float(coeff2[1])
    c4[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img5[0].header.index('WV_0_01')
while i<img5[0].header.index('WV_0_01')+26:
    s1=img5[0].header[i]
    s2=img5[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c5[z,0]=float(coeff1[0])
    c5[z,1]=float(coeff1[1])
    c5[z,2]=float(coeff1[2])
    c5[z,3]=float(coeff1[3])
    c5[z,4]=float(coeff2[0])
    c5[z,5]=float(coeff2[1])
    c5[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img6[0].header.index('WV_0_01')
while i<img6[0].header.index('WV_0_01')+16:
    s1=img6[0].header[i]
    s2=img6[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c6[z,0]=float(coeff1[0])
    c6[z,1]=float(coeff1[1])
    c6[z,2]=float(coeff1[2])
    c6[z,3]=float(coeff1[3])
    c6[z,4]=float(coeff2[0])
    c6[z,5]=float(coeff2[1])
    c6[z,6]=float(coeff2[2])
    i+=2
    z+=1

#Calculate wavelength value at each pixel
for i in range(4056):
    for j in range(20):
        wv1[j,i]=c1[j,0]+c1[j,1]*i+c1[j,2]*i**2+c1[j,3]*i**3+c1[j,4]*i**4 \
            +c1[j,5]*i**5+c1[j,6]*i**6
    for j in range(13):
        wv2[j,i]=c2[j,0]+c2[j,1]*i+c2[j,2]*i**2+c2[j,3]*i**3+c2[j,4]*i**4 \
            +c2[j,5]*i**5+c2[j,6]*i**6
        wv5[j,i]=c5[j,0]+c5[j,1]*i+c5[j,2]*i**2+c5[j,3]*i**3+c5[j,4]*i**4 \
            +c5[j,5]*i**5+c5[j,6]*i**6
    for j in range(8):
        wv3[j,i]=c3[j,0]+c3[j,1]*i+c3[j,2]*i**2+c3[j,3]*i**3+c3[j,4]*i**4 \
            +c3[j,5]*i**5+c3[j,6]*i**6
        wv6[j,i]=c6[j,0]+c6[j,1]*i+c6[j,2]*i**2+c6[j,3]*i**3+c6[j,4]*i**4 \
            +c6[j,5]*i**5+c6[j,6]*i**6
    for j in range(21):
        wv4[j,i]=c4[j,0]+c4[j,1]*i+c4[j,2]*i**2+c4[j,3]*i**3+c4[j,4]*i**4 \
            +c4[j,5]*i**5+c4[j,6]*i**6
    
#Continuum normalize each order and interpolate over common wavelength scale
#for each order
for i in range(len(data1)):
    p=polyfit(wv1[i,:],data1[i,:])
    for j in range(4056):
        data1[i,j]/=wv1[i,j]*p[0]+p[1]
        data1[i,j]-=1
    f=interp1d(wv1[i,:],data1[i,:])
    wv=wv1[i,0]*(1.+r)**linspace(0.,4055.,4056)
    science1[i]=f(wv)
for i in range(len(data2)):
    p=polyfit(wv1[i,:],data2[i,:])
    for j in range(4056):
        data2[i,j]/=wv2[i,j]*p[0]+p[1]
        data2[i,j]-=1
    f=interp1d(wv2[i,:],data2[i,:])
    wv=wv2[i,0]*(1.+r)**linspace(0.,4055.,4056)
    science2[i]=f(wv)
for i in range(len(data3)):
    p=polyfit(wv3[i,:],data3[i,:])
    for j in range(4056):
        data3[i,j]/=wv3[i,j]*p[0]+p[1]
        data3[i,j]-=1
    f=interp1d(wv3[i,:],data3[i,:])
    wv=wv3[i,0]*(1.+r)**linspace(0.,4055.,4056)
    science3[i]=f(wv)
for i in range(len(data4)):
    p=polyfit(wv4[i,:],data4[i,:])
    for j in range(4056):
        data4[i,j]/=wv4[i,j]*p[0]+p[1]
        data4[i,j]-=1
    f=interp1d(wv4[i,:],data4[i,:])
    wv=wv4[i,0]*(1.+r)**linspace(0.,4055.,4056)
    standard1[i]=f(wv)
for i in range(len(data5)):
    p=polyfit(wv5[i,:],data5[i,:])
    for j in range(4056):
        data5[i,j]/=wv5[i,j]*p[0]+p[1]
        data5[i,j]-=1
    f=interp1d(wv5[i,:],data5[i,:])
    wv=wv5[i,0]*(1.+r)**linspace(0.,4055.,4056)
    standard2[i]=f(wv)
for i in range(len(data6)):
    p=polyfit(wv6[i,:],data6[i,:])
    for j in range(4056):
        data6[i,j]/=wv6[i,j]*p[0]+p[1]
        data6[i,j]-=1
    f=interp1d(wv6[i,:],data6[i,:])
    wv=wv6[i,0]*(1.+r)**linspace(0.,4055.,4056)
    standard3[i]=f(wv)

##p1=polyfit(wv1[2,:],data1[2,:],1)
##p2=polyfit(wv4[3,:],data4[3,:],1)
##for i in range(4056):
##	data1[2,i]/=wv1[2,i]*p1[0]+p1[1]
##	data1[2,i]-=1
##	data4[3,i]/=wv4[3,i]*p2[0]+p2[1]
##	data4[3,i]-=1


##cc=correlate(science,standard,mode='full')
##v=range(-4055,4056)
##for i in range(-4055,4056):
##    v[i]*=1.28
##plot(v,cc)
##xlabel('Velocity Shift (km/s)')
##ylabel('Correlation Coefficient')
##title('CC of GJ 3236 A with GJ 96')
##savefig('cc.pdf')
##close()
##
##plot(v[3055:5057],cc[3055:5057])
##xlabel('Velocity Shift (km/s)')
##ylabel('Correlation Coefficient')
##title('CC of GJ 3236 A with GJ 96')
##savefig('cc_short.pdf')
##close()


##plot(wv1[2,:],data1[2,:],label='GJ 3236 A (Science)')
##plot(wv4[3,:],data4[3,:],label='GJ 96 (Standard)')
##xlabel('Wavelength [A]')
##ylabel('Normalized Flux')
##legend(loc='best')
##savefig('Order_Plot.pdf')
##close()


##for i in range(20):
##	plot(wv1[i,:],data1[i,:],'b')
##for i in range(13):
##	plot(wv2[i,:],data2[i,:],'g')
##for i in range(8):
##	plot(wv3[i,:],data3[i,:],'r')
##xlabel('Wavelength [A]')
##ylabel('Counts')
##title('GJ 3236')
##savefig('GJ3236_spec.pdf')
##close()
##
