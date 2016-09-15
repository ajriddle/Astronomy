from numpy import *

def calc_velocity(a,P,e,g,w,q,T,t):
    """
    Funtion to calculate radial velocities of binary components at a given epoch
    (t), given the orbital elements {a,P,e,g,w,q,T}
    inputs: a[AU],P[days],e,g[km/s],w,q,T[HJD],t[HJD]
    outputs: v1,v2 (km/s)
    """

    v1=[]
    v2=[]
    for i in range(len(t)):
        #Mean anomaly M = 2*pi(t-T)/P
        M=2*pi*((t[i]-T)%P)/P

        #Eccentric anomaly M = E - e sin(E)
        #Solve transcendental equation using g(E)=E-e sin(E)-M=0,
        #E_i+1=E_i-g(E_i)/g'(E_i),E_0=M
        E=zeros(10)
        E[0]=M
        for j in range(9):
            E[j+1]=E[j]-(E[j]-e*sin(E[j])-M)/(1-e*cos(E[j]))
        E=E[-1]

        #True anomaly (theta)
        theta=2*arctan2(sqrt((1+e)/(1-e))*sin(E/2),cos(E/2))

        #Calculate a1, semi-major axis of primary
        Msum=a**3/(P/365.242)**2 *1.99E30 #kg
        M1=Msum/(1+q)
        M2=Msum*q/(1+q)
        a2=M1/Msum*a*1.496E8 #km
        v2.append(2*pi*a2/(P*24*3600*sqrt(1-e**2))*(cos(theta+w)+e*cos(w))+g) #km/s
        v1.append(-q*v2[-1]+(1+q)*g)

    return v1,v2

