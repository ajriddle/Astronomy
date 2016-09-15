from pyfits import *
from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import *


def rmspace(var):
    while '' in var:
        var.pop(var.index(''))
    return var

#Science target 
img1=open('HIRES_sci_27480_39/extracted/makee/ccd1/fits/' \
    'HI.20080813.53907_1_Flux.fits')
data1=img1[0].data
img1.close()
img2=open('HIRES_sci_27480_39/extracted/makee/ccd2/fits/' \
    'HI.20080813.53907_2_Flux.fits')
data2=img2[0].data
img2.close()
img3=open('HIRES_sci_27480_39/extracted/makee/ccd3/fits/' \
    'HI.20080813.53907_3_Flux.fits')
data3=img3[0].data
img3.close()

#Template 1 
img4=open('HIRES_sci_9812_1/extracted/makee/ccd1/fits/' \
    'HI.20101030.40717_1_Flux.fits')
data4=img4[0].data
img4.close()
img5=open('HIRES_sci_9812_1/extracted/makee/ccd2/fits/' \
    'HI.20101030.40717_2_Flux.fits')
data5=img5[0].data
img5.close()
img6=open('HIRES_sci_9812_1/extracted/makee/ccd3/fits/' \
    'HI.20101030.40717_3_Flux.fits')
data6=img6[0].data
img6.close()

#Template 2 
img7=open('HIRES_sci_9812_1/extracted/makee/ccd1/fits/' \
    'HI.20101030.40717_1_Flux.fits')
data7=img7[0].data
img7.close()
img8=open('HIRES_sci_9812_1/extracted/makee/ccd2/fits/' \
    'HI.20101030.40717_2_Flux.fits')
data8=img8[0].data
img8.close()
img9=open('HIRES_sci_9812_1/extracted/makee/ccd3/fits/' \
    'HI.20101030.40717_3_Flux.fits')
data9=img8[0].data
img9.close()

wv1=zeros([(img1[0].header.index('VACUUM')-img1[0].header.index('WV_0_01'))/2, \
           4056])
#wv2=zeros([(img2[0].header.index('VACUUM')-img2[0].header.index('WV_0_01'))/2, \
#           4056])
wv3=zeros([(img3[0].header.index('VACUUM')-img3[0].header.index('WV_0_01'))/2, \
           4056])
wv4=zeros([(img4[0].header.index('VACUUM')-img4[0].header.index('WV_0_01'))/2, \
           4056])
wv5=zeros([(img5[0].header.index('VACUUM')-img5[0].header.index('WV_0_01'))/2, \
           4056])
wv6=zeros([(img6[0].header.index('VACUUM')-img6[0].header.index('WV_0_01'))/2, \
           4056])
wv7=zeros([(img7[0].header.index('VACUUM')-img7[0].header.index('WV_0_01'))/2, \
           4056])
wv8=zeros([(img8[0].header.index('VACUUM')-img8[0].header.index('WV_0_01'))/2, \
           4056])
wv9=zeros([(img9[0].header.index('VACUUM')-img9[0].header.index('WV_0_01'))/2, \
           4056])
c1=zeros([len(wv1),7])
#c2=zeros([len(wv2),7])
c3=zeros([len(wv3),7])
c4=zeros([len(wv4),7])
c5=zeros([len(wv5),7])
c6=zeros([len(wv6),7])
c7=zeros([len(wv7),7])
c8=zeros([len(wv8),7])
c9=zeros([len(wv9),7])

z=0
i=img1[0].header.index('WV_0_01')
while img1[0].header.index('WV_0_01')<=i<img1[0].header.index('VACUUM'):
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

##z=0
##i=img2[0].header.index('WV_0_01')
##while img2[0].header.index('WV_0_01')<=i<img2[0].header.index('VACUUM'):
##    s1=img2[0].header[i]
##    s2=img2[0].header[i+1]
##    coeff1=rmspace(s1.split(' '))
##    coeff2=rmspace(s2.split(' '))
##    c2[z,0]=float(coeff1[0])
##    c2[z,1]=float(coeff1[1])
##    c2[z,2]=float(coeff1[2])
##    c2[z,3]=float(coeff1[3])
##    c2[z,4]=float(coeff2[0])
##    c2[z,5]=float(coeff2[1])
##    c2[z,6]=float(coeff2[2])
##    i+=2
##    z+=1

z=0
i=img3[0].header.index('WV_0_01')
while img3[0].header.index('WV_0_01')<=i<img3[0].header.index('VACUUM'):
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
while img4[0].header.index('WV_0_01')<=i<img4[0].header.index('VACUUM'):
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
while img5[0].header.index('WV_0_01')<=i<img5[0].header.index('VACUUM'):
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
while img6[0].header.index('WV_0_01')<=i<img6[0].header.index('VACUUM'):
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

z=0
i=img7[0].header.index('WV_0_01')
while img7[0].header.index('WV_0_01')<=i<img7[0].header.index('VACUUM'):
    s1=img7[0].header[i]
    s2=img7[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c7[z,0]=float(coeff1[0])
    c7[z,1]=float(coeff1[1])
    c7[z,2]=float(coeff1[2])
    c7[z,3]=float(coeff1[3])
    c7[z,4]=float(coeff2[0])
    c7[z,5]=float(coeff2[1])
    c7[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img8[0].header.index('WV_0_01')
while img8[0].header.index('WV_0_01')<=i<img8[0].header.index('VACUUM'):
    s1=img8[0].header[i]
    s2=img8[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c8[z,0]=float(coeff1[0])
    c8[z,1]=float(coeff1[1])
    c8[z,2]=float(coeff1[2])
    c8[z,3]=float(coeff1[3])
    c8[z,4]=float(coeff2[0])
    c8[z,5]=float(coeff2[1])
    c8[z,6]=float(coeff2[2])
    i+=2
    z+=1

z=0
i=img9[0].header.index('WV_0_01')
while img9[0].header.index('WV_0_01')<=i<img9[0].header.index('VACUUM'):
    s1=img9[0].header[i]
    s2=img9[0].header[i+1]
    coeff1=rmspace(s1.split(' '))
    coeff2=rmspace(s2.split(' '))
    c9[z,0]=float(coeff1[0])
    c9[z,1]=float(coeff1[1])
    c9[z,2]=float(coeff1[2])
    c9[z,3]=float(coeff1[3])
    c9[z,4]=float(coeff2[0])
    c9[z,5]=float(coeff2[1])
    c9[z,6]=float(coeff2[2])
    i+=2
    z+=1

#Calculate wavelength value at each pixel
for i in range(4056):
    for j in range(len(wv1)):
        wv1[j,i]=c1[j,0]+c1[j,1]*i+c1[j,2]*i**2+c1[j,3]*i**3+c1[j,4]*i**4 \
            +c1[j,5]*i**5+c1[j,6]*i**6
##    for j in range(len(wv2)):
##        wv2[j,i]=c2[j,0]+c2[j,1]*i+c2[j,2]*i**2+c2[j,3]*i**3+c2[j,4]*i**4 \
##            +c2[j,5]*i**5+c2[j,6]*i**6
    for j in range(len(wv5)):
        wv5[j,i]=c5[j,0]+c5[j,1]*i+c5[j,2]*i**2+c5[j,3]*i**3+c5[j,4]*i**4 \
            +c5[j,5]*i**5+c5[j,6]*i**6
    for j in range(len(wv3)):
        wv3[j,i]=c3[j,0]+c3[j,1]*i+c3[j,2]*i**2+c3[j,3]*i**3+c3[j,4]*i**4 \
            +c3[j,5]*i**5+c3[j,6]*i**6
    for j in range(len(wv6)):
        wv6[j,i]=c6[j,0]+c6[j,1]*i+c6[j,2]*i**2+c6[j,3]*i**3+c6[j,4]*i**4 \
            +c6[j,5]*i**5+c6[j,6]*i**6
    for j in range(len(wv4)):
        wv4[j,i]=c4[j,0]+c4[j,1]*i+c4[j,2]*i**2+c4[j,3]*i**3+c4[j,4]*i**4 \
            +c4[j,5]*i**5+c4[j,6]*i**6
    for j in range(len(wv7)):
        wv7[j,i]=c7[j,0]+c7[j,1]*i+c7[j,2]*i**2+c7[j,3]*i**3+c7[j,4]*i**4 \
            +c7[j,5]*i**5+c7[j,6]*i**6
    for j in range(len(wv8)):
        wv8[j,i]=c8[j,0]+c8[j,1]*i+c8[j,2]*i**2+c8[j,3]*i**3+c8[j,4]*i**4 \
            +c8[j,5]*i**5+c8[j,6]*i**6
    for j in range(len(wv9)):
        wv9[j,i]=c9[j,0]+c9[j,1]*i+c9[j,2]*i**2+c9[j,3]*i**3+c9[j,4]*i**4 \
            +c9[j,5]*i**5+c9[j,6]*i**6
    
#Stack orders from each chip into one array
wvs=vstack(([wv1,wv3]))
wvt1=vstack(([wv4,wv5,wv6]))
wvt2=vstack(([wv7,wv8,wv9]))
datas=vstack(([data1,data3]))
datat1=vstack(([data4,data5,data6]))
datat2=vstack(([data7,data8,data9]))

#Continuum normalize each order and interpolate over common wavelength scale
#for each order
r=1.28/2.998E5

#Convert to velocity units    
v=range(-4055,4056)
for i in range(-4055,4056):
    v[i]=round(v[i]*1.28,2)
    
science=zeros(1)
temp1=zeros(4056)
temp2=zeros(4056)


for i in range(len(wvs)):
    for j in range(len(wvt1)):
        if abs(wvs[i,0]-wvt1[j,0])<10:
            for k in range(len(wvt2)):
                if abs(wvs[i,0]-wvt2[k,0])<10:
                    diff1=10
                    diff2=10
                    index1=0
                    index2=0
                    for l in range(4056):
                        if abs(wvs[i,0]-wvt1[j,l])<diff1:
                            diff1=abs(wvs[i,0]-wvt1[j,l])
                            index1=l
                        if abs(wvs[i,0]-wvt2[k,l])<diff2:
                            diff2=abs(wvs[i,0]-wvt2[k,l])
                            index2=l

wvs=wvs[:,:-max([index1,index2])+1]
wvt1=wvt1[:,max([index1,index2]):]
wvt2=wvt2[:,max([index1,index2]):]

for i in range(len(wvs)):
    for j in range(len(wvt1)):
        if abs(wvs[i,0]-wvt1[j,0])<10:
            for k in range(len(wvt2)):
                if abs(wvs[i,0]-wvt2[k,0])<10:                    
                    ps=polyfit(wvs[i,:],datas[i,:],1)
                    pt1=polyfit(wvt1[j,:],datat1[j,:],1)
                    pt2=polyfit(wvt2[k,:],datat2[k,:],1)
                    for l in range(4056):
                        datas[i,l]/=wvs[i,l]*ps[0]+ps[1]
                        datas[i,l]-=1
                        datat1[j,l]/=wvt1[j,l]*pt1[0]+pt1[1]
                        datat1[j,l]-=1
                        datat2[k,l]/=wvt2[k,l]*pt2[0]+pt2[1]
                        datat2[k,l]-=1
                    fs=interp1d(wvs[i,:],datas[i,:])
                    ft1=interp1d(wvt1[j,:],datat1[j,:])
                    ft2=interp1d(wvt2[k,:],datat2[k,:])
                    wv=wvs[i,0]*(1.+r)**linspace(0.,4055.,4056)
                    if len(science)>1:
                        science=vstack((science,fs(wv)))
                        temp1=vstack((temp1,ft1(wv)))
                        temp2=vstack((temp2,ft2(wv)))
                    elif len(science)==1:
                        science=fs(wv)
                        temp1=ft1(wv)
                        temp2=ft2(wv)
                    break
            break
                

##for i in range(len(wv2)):
##    for j in range(len(wv5)):
##        if abs(wv2[i,0]-wv5[j,0])<2:
##            for k in range(len(wv8)):
##                if abs(wv2[i,0]-wv8[k,0])<2:
##                    ps=polyfit(wv2[i,:],data2[i,:],1)
##                    pt1=polyfit(wv5[j,:],data5[j,:],1)
##                    pt2=polyfit(wv8[k,:],data8[k,:],1)
##                    for l in range(4056):
##                        data2[i,l]/=wv2[i,l]*ps[0]+ps[1]
##                        data2[i,l]-=1
##                        data5[j,l]/=wv5[j,l]*pt1[0]+pt1[1]
##                        data5[j,l]-=1
##                        data8[k,l]/=wv8[k,l]*pt2[0]+pt2[1]
##                        data8[k,l]-=1
##                    fs=interp1d(wv2[i,:],data2[i,:])
##                    ft1=interp1d(wv5[j,:],data5[j,:])
##                    ft2=interp1d(wv8[k,:],data8[k,:])
##                    wv=wv2[i,0]*(1.+r)**linspace(0.,4055.,4056)
##                    science=vstack((science,fs(wv)))
##                    temp1=vstack((temp1,ft1(wv)))
##                    temp2=vstack((temp2,ft2(wv)))
##                    break
##            break
        
#cc=correlate(science[2],temp1[2],mode='full')
