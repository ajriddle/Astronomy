from numpy import *

v=zeros(24) #Valenti & Fischer (2005) CA data
s=zeros(22) #Sousa et al. (2011) HARPS volume-limited data
n=zeros(41) #Neves et al. (2009) HARPS GTO data
m=zeros(16) #Sousa et al. (2011) HARPS metal-poor sample
h=loadtxt('hip.tsv',dtype=str,delimiter=';') #Hipparcos catalog

f=open('V&F2005.tsv')
for line in f:
    if '#' not in line:
        r=line.split(';')
        l=r[1].split(' ')
        while '' in l:
            l.pop(l.index(''))
        t=''
        for i in range(len(l)):
            t+=l[i]
        r[1]=t
        v=vstack([v,r])
v=v[1:,:]
f.close()

f=open('Sousa2011.tsv')
for line in f:
    if '#' not in line:
        r=line.split(';')
        if ' ' in r[0]:
            l=r[0].split(' ')
            while '' in l:
                l.pop(l.index(''))
            l=l[0]+l[1]
            r[0]=l
        s=vstack([s,r])
s=s[1:,:]
f.close()


f=open('low-Z_HARPS.tsv')
for line in f:
    if '#' not in line:
        r=line.split(';')
        l=r[0].split(' ')
        while '' in l:
            l.pop(l.index(''))
        r[0]=l[0]
        m=vstack([m,r])
m=m[1:,:]
f.close()

f=open('Neves2009.tsv')
for line in f:
    if '#' not in line:
        r=line.split(';')
        if ' ' in r[2]:
            l=r[2].split(' ')
            while '' in l:
                l.pop(l.index(''))
            l=l[0]+l[1]
            r[2]=l
        n=vstack([n,r])
n=n[1:,:]
f.close()


v=column_stack([v,zeros([len(v),shape(h)[1]],dtype=list)])
for i in range(len(v)):
    if v[i,1][2:] in h[:,28]:
        v[i,-shape(h)[1]:]=h[h[:,28].index(v[i,1][2:]),:]

#total=0
#l=list(v[:,1])+list(s[:,0])+list(m[:,0])+list(n[:,2])
#for i in range(len(l)):
 #       total+=1.0/l.count(l[i])





