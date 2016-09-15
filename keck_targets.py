from numpy import *

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

f=genfromtxt('keck_targets.csv',dtype='string',delimiter=',',autostrip=True, \
skiprows=1)

i=0
while i <len(f):
    if ('HD' in f[i,1]) or ('GJ' in f[i,1]) or ('HIP' in f[i,1]) or  \
    ('hd' in f[i,1]) or ('hip' in f[i,1]) or ('gj' in f[i,1]) or \
    ('gl' in f[i,1]) or is_number(f[i,1])==True:
        f=vstack(([f[:i],f[i+1:]]))       
    else:
        i+=1
#a=open('keck_trimmed.csv')
savetxt('keck_trimmed.csv',f,fmt='%s',delimiter=',')
