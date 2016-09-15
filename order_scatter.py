from numpy import *
from matplotlib.pyplot import *
import scipy.stats as stat

d=genfromtxt('rv_results',dtype=list)
rv={}
epoch=[]
v1=[]
v2=[]
e1=[]
e2=[]

for i in range(len(d)):
    if d[i,-1]=='1':
        key=d[i,2]
        rv[key]={'v':zeros(5),'m':[0.,0.],'e':[0.,0.],'bo':[]}
    else:
        for j in range(5):
            d[i,j]=float(d[i,j])
        rv[key]['v']=vstack(([rv[key]['v'],d[i]]))
    
for key in rv:
    rv[key]['v']=rv[key]['v'][1:]
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
            elif abs(rv[key]['v'][i,0]-6553)<10 or abs(rv[key]['v'][i,0]-8229)\
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
    rv[key]['m']=[v1[-1],v2[-1]]
    rv[key]['e']=[e1[-1],e2[-1]]
    if key=='54467.371626':
        v1[-1],v2[-1]=v2[-1],v1[-1]
        e1[-1],e2[-1]=e2[-1],e1[-1]
    if key=='54690.482060':
        v2[-1]=-56.26

separation=15
o=[0]
for key in rv.keys():
    for i in range(len(rv[key]['v'])):
        for j in range(len(o)):
            if abs(rv[key]['v'][i,0]-o[j])<separation:
                break
            elif j==len(o)-1:
                o.append(rv[key]['v'][i,0])
            
o.sort()
o=o[1:]
orders={}
for key in o:
    orders[key]=[0,0]

for key in rv.keys():
    for i in range(len(rv[key]['v'])):
        for wave in orders:
            if abs(float(wave)-rv[key]['v'][i,0])<separation:
                delta1=rv[key]['v'][i,1]-rv[key]['m'][0]
                delta2=rv[key]['v'][i,3]-rv[key]['m'][1]
                orders[wave]=vstack((orders[wave],[delta1,delta2]))
                break
x=[]
y=[0,0]
e=[0,0]
print len(o)
#o=vstack((o,zeros([4,len(o)])))
for key in orders:
    m1=sum(orders[key][:,0])/len(orders[key])
    m2=sum(orders[key][:,1])/len(orders[key])
    s1=0.
    s2=0.
    for i in range(len(orders[key])):
        s1+=(orders[key][i,0]-m1)**2
        s2+=(orders[key][i,1]-m2)**2
    s1/=(len(orders[key])+1)
    s2/=(len(orders[key])+1)
    s1=s1**.5
    s2=s2**.5
    x.append(float(key))
    y=vstack((y,[m1,m2]))
    e=vstack((e,[s1,s2]))
y=y[1:]
e=e[1:]

errorbar(x,y[:,0],e[:,0],fmt='none',ecolor='blue',label='Primary')
errorbar(x,y[:,1],e[:,1],fmt='none',ecolor='red',label='Secondary')
hlines(0,4000,9000,linestyles='dashed')
xlabel(r'Wavelength ($\AA$)')
ylabel('Residual (km/s)')
legend(loc='best')
savefig('Order_Residuals.png')
close()

f=open('order_residuals.txt','w')
for i in range(len(x)):
    print >>f, x[i],'\t',y[i,0],'\t',e[i,0],'\t',y[i,1],'\t',e[i,1],'\t',\
        len(orders[x[i]])
f.close()

