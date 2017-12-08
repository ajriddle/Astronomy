from numpy import *
from matplotlib.pyplot import *
from matplotlib.cm import *

eb=genfromtxt('EB_Characterization.csv',dtype='string',delimiter=',',skiprows=1)
time=[2456854.5]
while time[-1]<2456857.0:
    time.append(time[-1]+1./24)
phases=zeros([len(eb),len(time)])
for i in range(22,len(eb)):
    if eb[i,23]!='':
        for j in range(len(time)):
            phases[i,j]=((time[j]-float(eb[i,23]))%float(eb[i,2]))/float(eb[i,2])
for i in range(22,len(phases)):
    if eb[i,23]!='':
        for j in range(len(time)):
            if (0.4<=phases[i,j]<=0.6 or phases[i,j]>=0.9 or phases[i,j]<0.1) \
            and .583<=time[j]-int(time[j])<=.875:
                scatter(time[j],phases[i,j])
        xlabel('HJD')
        ylabel('Phase')
        title(eb[i,0])
        savefig('%s_Phases.pdf' % eb[i,0])
        close()

for i in range(22,len(phases)):
    if eb[i,23]!='':
        for j in range(len(time)):
            if (0.4<=phases[i,j]<=0.6 or phases[i,j]>=0.9 or phases[i,j]<0.1) \
            and .583<=time[j]-int(time[j])<=.875:
                scatter(time[j],phases[i,j],c=gist_rainbow(float(i-22)/14))
xlabel('HJD')
ylabel('Phase')
title('Phases')
savefig('Phases.pdf')
close()
