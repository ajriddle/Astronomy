from numpy import *

##f22=genfromtxt('febobs22.csv',delimiter=',',skiprows=1,dtype='string')
##f23=genfromtxt('febobs23.csv',delimiter=',',skiprows=1,dtype='string')
##f24=genfromtxt('febobs24.csv',delimiter=',',skiprows=1,dtype='string')
##f26=genfromtxt('febobs26.csv',delimiter=',',skiprows=1,dtype='string')
##f27=genfromtxt('febobs27.csv',delimiter=',',skiprows=1,dtype='string')
##f28=genfromtxt('febobs28.csv',delimiter=',',skiprows=1,dtype='string')
##log21=genfromtxt('2013-02-21-log.txt',dtype='string')
##log22=genfromtxt('2013-02-22-log.txt',dtype='string')
##log23=genfromtxt('2013-02-23-log.txt',dtype='string')
##log25=genfromtxt('2013-02-25-log.txt',dtype='string')
##log26=genfromtxt('2013-02-26-log.txt',dtype='string')
##log27=genfromtxt('2013-02-27-log.txt',dtype='string')

##m20=genfromtxt('march20obs.csv',delimiter=',',skiprows=1,dtype='string')
##m21=genfromtxt('march21obs.csv',delimiter=',',skiprows=1,dtype='string')
##m22=genfromtxt('march22obs.csv',delimiter=',',skiprows=1,dtype='string')
##m23=genfromtxt('march23obs.csv',delimiter=',',skiprows=1,dtype='string')
##m25=genfromtxt('march25obs.csv',delimiter=',',skiprows=1,dtype='string')
##log19=genfromtxt('2013-03-19-log.txt',dtype='string')
##log20=genfromtxt('2013-03-20-log.txt',dtype='string')
##log21=genfromtxt('2013-03-21-log.txt',dtype='string')
##log22=genfromtxt('2013-03-22-log.txt',dtype='string')
##log24=genfromtxt('2013-03-24-log.txt',dtype='string')

a13=genfromtxt('aprobs13.txt',delimiter='\t',skiprows=1,dtype='string')
a14=genfromtxt('aprobs14.txt',delimiter='\t',skiprows=1,dtype='string')
a15=genfromtxt('aprobs15.txt',delimiter='\t',skiprows=1,dtype='string')
a16=genfromtxt('aprobs16.txt',delimiter='\t',skiprows=1,dtype='string')
a17=genfromtxt('aprobs17.txt',delimiter='\t',skiprows=1,dtype='string')
a19=genfromtxt('aprobs19.txt',delimiter='\t',skiprows=1,dtype='string')
a20=genfromtxt('aprobs20.txt',delimiter='\t',skiprows=1,dtype='string')
a21=genfromtxt('aprobs21.txt',delimiter='\t',skiprows=1,dtype='string')
a22=genfromtxt('aprobs22.txt',delimiter='\t',skiprows=1,dtype='string')
log12=genfromtxt('2013-04-12-log.txt',dtype='string')
log13=genfromtxt('2013-04-13-log.txt',dtype='string')
log14=genfromtxt('2013-04-14-log.txt',dtype='string')
log15=genfromtxt('2013-04-15-log.txt',dtype='string')
log16=genfromtxt('2013-04-16-log.txt',dtype='string')
log18=genfromtxt('2013-04-18-log.txt',dtype='string')
log19=genfromtxt('2013-04-19-log.txt',dtype='string')
log20=genfromtxt('2013-04-20-log.txt',dtype='string')
log21=genfromtxt('2013-04-21-log.txt',dtype='string')



f=open('imedit_script','w')
for i in range(len(log12)):
    if log12[i,3] in a13[:,0]:
        n=list(a13[:,0]).index(log12[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log12[i,0],a13[n,3],a13[n,4],a13[n,6])
    else:
        print '1',log12[i,3],log12[i,0]
    
for i in range(len(log13)):
    if log13[i,3] in a14[:,0]:
        n=list(a14[:,0]).index(log13[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log13[i,0],a14[n,3],a14[n,4],a14[n,6])
    else:
        print '2',log13[i,3],log13[i,0]

for i in range(len(log14)):
    if log14[i,3] in a15[:,0]:
        n=list(a15[:,0]).index(log14[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log14[i,0],a15[n,3],a15[n,4],a15[n,6])
    else:
        print '3',log14[i,3],log14[i,0]

for i in range(len(log15)):
    if log15[i,3] in a16[:,0]:
        n=list(a16[:,0]).index(log15[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log15[i,0],a16[n,3],a16[n,4],a16[n,6])
    else:
        print '4',log15[i,3],log15[i,0]

for i in range(len(log16)):
    if log16[i,3] in a17[:,0]:
        n=list(a17[:,0]).index(log16[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log16[i,0],a17[n,3],a17[n,4],a17[n,6])
    else:
        print '5',log16[i,3],log16[i,0]

for i in range(len(log18)):
    if log18[i,3] in a19[:,0]:
        n=list(a19[:,0]).index(log18[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log18[i,0],a19[n,3],a19[n,4],a19[n,6])
    else:
        print '6',log18[i,3],log18[i,0]

for i in range(len(log19)):
    if log19[i,3] in a20[:,0]:
        n=list(a20[:,0]).index(log19[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log19[i,0],a20[n,3],a20[n,4],a20[n,6])
    else:
        print '7',log19[i,3],log119[i,0]

for i in range(len(log20)):
    if log20[i,3] in a21[:,0]:
        n=list(a21[:,0]).index(log20[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log20[i,0],a21[n,3],a21[n,4],a21[n,6])
    else:
        print '8',log20[i,3],log20[i,0]

for i in range(len(log21)):
    if log21[i,3] in a22[:,0]:
        n=list(a22[:,0]).index(log21[i,3])
        print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
        (log21[i,0],a22[n,3],a22[n,4],a22[n,6])
    else:
        print '9',log21[i,3],log21[i,0]

##for i in range(len(log21)):
##    n=list(f22[:,0]).index(log21[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log21[i,0],f22[n,3],f22[n,4],f22[n,6])
##    
##for i in range(len(log22)):
##    n=list(f23[:,0]).index(log22[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log22[i,0],f23[n,3],f23[n,4],f23[n,6])
##
##for i in range(len(log23)):
##    n=list(f24[:,0]).index(log23[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log23[i,0],f24[n,3],f24[n,4],f24[n,6])
##
##for i in range(len(log25)):
##    n=list(f26[:,0]).index(log25[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log25[i,0],f26[n,3],f26[n,4],f26[n,6])
##
##for i in range(len(log26)):
##    n=list(f27[:,0]).index(log26[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log26[i,0],f27[n,3],f27[n,4],f27[n,6])
##
##for i in range(len(log21)):
##    n=list(f22[:,0]).index(log21[i,3])
##    print >>f, 'plug image=%s rrate=%s drate=%s ptime=%s' % \
##    (log21[i,0],f22[n,3],f22[n,4],f22[n,6])

f.close()
