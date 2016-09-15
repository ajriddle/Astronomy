from numpy import *
from matplotlib.pyplot import *

s=zeros(4)
c1=zeros(4)
c2=zeros(4)
c3=zeros(4)
t=[]
x=[]
y=[]
r=genfromtxt('20140511/results',dtype=list)
for i in range(len(r)/4):
    s=vstack(([s,r[4*i,2:6]]))
    c1=vstack(([c1,r[4*i+1,2:6]]))
    c2=vstack(([c2,r[4*i+2,2:6]]))
    c3=vstack(([c3,r[4*i+3,2:6]]))
    t.append(r[4*i,-1])
    x.append(float(r[4*i,0]))
    y.append(float(r[4*i,1]))
s=s[1:]
c1=c1[1:]
c2=c2[1:]
c3=c3[1:]
r1=[]
r2=[]
r3=[]
ar=[]
d1=[]
d2=[]
d3=[]
ad=[]
for i in range(len(s)):
    for j in range(4):
        s[i,j]=float(s[i,j])
        c1[i,j]=float(c1[i,j])
        c2[i,j]=float(c2[i,j])
        c3[i,j]=float(c3[i,j])
    k=t[i].split(':')
    t[i]=(float(k[0])+float(k[1])/60.+float(k[2])/3600.+17-24)/24.+788.5
    r1.append(s[i,0]/c1[i,0])
    r2.append(s[i,0]/c2[i,0])
    r3.append(s[i,0]/c3[i,0])
    ar.append(average([r1[-1],r2[-1],r3[-1]]))
    d1.append(-s[i,2]+c1[i,2])
    d2.append(-s[i,2]+c2[i,2])
    d3.append(-s[i,2]+c3[i,2])
    ad.append(average([d1[-1],d2[-1],d3[-1]]))
#scatter(t,d1,c='k',marker='+')
#xlabel('UT (hr)')
#ylabel(r'$m_{comparison}-m_{target}$')
#savefig('NSVS1082016_20140511_Photometry.pdf')
#close()
t511=t
d511=d1

s=zeros(4)
c1=zeros(4)
c2=zeros(4)
c3=zeros(4)
t=[]
x=[]
y=[]
r=genfromtxt('20140716/results',dtype=list)
for i in range(len(r)/4):
    s=vstack(([s,r[4*i,2:6]]))
    c1=vstack(([c1,r[4*i+1,2:6]]))
    c2=vstack(([c2,r[4*i+2,2:6]]))
    c3=vstack(([c3,r[4*i+3,2:6]]))
    t.append(r[4*i,-1])
    x.append(float(r[4*i,0]))
    y.append(float(r[4*i,1]))
s=s[1:]
c1=c1[1:]
c2=c2[1:]
c3=c3[1:]
r1=[]
r2=[]
r3=[]
ar=[]
d1=[]
d2=[]
d3=[]
ad=[]
for i in range(len(s)):
    for j in range(4):
        s[i,j]=float(s[i,j])
        c1[i,j]=float(c1[i,j])
        c2[i,j]=float(c2[i,j])
        c3[i,j]=float(c3[i,j])
    k=t[i].split(':')
    t[i]=(float(k[0])+float(k[1])/60.+float(k[2])/3600.+17-24)/24.+854.5
    r1.append(s[i,0]/c1[i,0])
    r2.append(s[i,0]/c2[i,0])
    r3.append(s[i,0]/c3[i,0])
    ar.append(average([r1[-1],r2[-1],r3[-1]]))
    d1.append(-s[i,2]+c1[i,2])
    d2.append(-s[i,2]+c2[i,2])
    d3.append(-s[i,2]+c3[i,2])
    ad.append(average([d1[-1],d2[-1],d3[-1]]))
#scatter(t,d1,c='k',marker='+')
#xlabel('UT (hr)')
#ylabel(r'$m_{comparison}-m_{target}$')
#savefig('NSVS1082016_20140511_Photometry.pdf')
#close()
t716=t
d716=d1

s=zeros(4)
c1=zeros(4)
c2=zeros(4)
c3=zeros(4)
t=[]
x=[]
y=[]
r=genfromtxt('20140806/results',dtype=list)
for i in range(len(r)/4):
    s=vstack(([s,r[4*i,2:6]]))
    c1=vstack(([c1,r[4*i+1,2:6]]))
    c2=vstack(([c2,r[4*i+2,2:6]]))
    c3=vstack(([c3,r[4*i+3,2:6]]))
    t.append(r[4*i,-1])
    x.append(float(r[4*i,0]))
    y.append(float(r[4*i,1]))
s=s[1:]
c1=c1[1:]
c2=c2[1:]
c3=c3[1:]
r1=[]
r2=[]
r3=[]
ar=[]
d1=[]
d2=[]
d3=[]
ad=[]
for i in range(len(s)):
    for j in range(4):
        s[i,j]=float(s[i,j])
        c1[i,j]=float(c1[i,j])
        c2[i,j]=float(c2[i,j])
        c3[i,j]=float(c3[i,j])
    k=t[i].split(':')
    t[i]=(float(k[0])+float(k[1])/60.+float(k[2])/3600.+17-24)/24.+875.5
    r1.append(s[i,0]/c1[i,0])
    r2.append(s[i,0]/c2[i,0])
    r3.append(s[i,0]/c3[i,0])
    ar.append(average([r1[-1],r2[-1],r3[-1]]))
    d1.append(-s[i,2]+c1[i,2])
    d2.append(-s[i,2]+c2[i,2])
    d3.append(-s[i,2]+c3[i,2])
    ad.append(average([d1[-1],d2[-1],d3[-1]]))
#scatter(t,d1,c='k',marker='+')
#xlabel('UT (hr)')
#ylabel(r'$m_{comparison}-m_{target}$')
#savefig('NSVS1082016_20140511_Photometry.pdf')
#close()
t806=t
d806=d1

s=zeros(4)
c1=zeros(4)
c2=zeros(4)
c3=zeros(4)
t=[]
x=[]
y=[]
r=genfromtxt('20140807/results',dtype=list)
for i in range(len(r)/4):
    s=vstack(([s,r[4*i,2:6]]))
    c1=vstack(([c1,r[4*i+1,2:6]]))
    c2=vstack(([c2,r[4*i+2,2:6]]))
    c3=vstack(([c3,r[4*i+3,2:6]]))
    t.append(r[4*i,-1])
    x.append(float(r[4*i,0]))
    y.append(float(r[4*i,1]))
s=s[1:]
c1=c1[1:]
c2=c2[1:]
c3=c3[1:]
r1=[]
r2=[]
r3=[]
ar=[]
d1=[]
d2=[]
d3=[]
ad=[]
c=[]
for i in range(len(s)):
    for j in range(4):
        s[i,j]=float(s[i,j])
        c1[i,j]=float(c1[i,j])
        c2[i,j]=float(c2[i,j])
        c3[i,j]=float(c3[i,j])
    k=t[i].split(':')
    t[i]=(float(k[0])+float(k[1])/60.+float(k[2])/3600.+17-24)/24.+876.5
    r1.append(s[i,0]/c1[i,0])
    r2.append(s[i,0]/c2[i,0])
    r3.append(s[i,0]/c3[i,0])
    ar.append(average([r1[-1],r2[-1],r3[-1]]))
    d1.append(-s[i,2]+c1[i,2])
    d2.append(-s[i,2]+c2[i,2])
    d3.append(-s[i,2]+c3[i,2])
    ad.append(average([d1[-1],d2[-1],d3[-1]]))
    c.append(c1[i,2]-c2[i,2])
scatter(t,d1,c='k',marker='+')
scatter(t,c,c='red',marker='+')
xlabel('UT (hr)')
ylabel(r'$m_{comparison}-m_{target}$')
savefig('NSVS1082016_20140807_Photometry.pdf')
close()
t807=t
d807=d1

t=t511+t716+t806+t807
d=d511+d716+d806+d807
scatter(t,d)
xlabel('JD-2456000')
ylabel('mag difference')
show()

