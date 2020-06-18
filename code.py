import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as pl
import os
os.environ["PROJ_LIB"]=r"C:\Users\Dell\anaconda3\Library\share"
from mpl_toolkits.basemap import Basemap


#Integrating Function
def func(f,t,e,k,u):
    return u*u*(1+e*np.cos(f))**2/k**3

def conv(i,theta,phi):
    sintheta_=-np.cos(theta)*np.cos(phi)*np.sin(i)+np.sin(theta)*np.cos(i)
    sintheta_[sintheta_>1]=1
    sintheta_[sintheta_<-1]=-1
    theta_=np.arcsin(sintheta_)
    
    cosphi_=(np.cos(theta)*np.cos(phi)*np.cos(i)+np.sin(theta)*np.sin(i))/np.cos(theta_)
    cosphi_[cosphi_>1]=1
    cosphi_[cosphi_<-1]=-1
    phi_=np.arccos(cosphi_)
    
    sinphi_=(np.cos(theta)*np.sin(phi))/np.cos(theta_)
    
    phi_[sinphi_<0]=2*np.pi-phi_[sinphi_<0]
    
    return theta_,phi_

#Inputs
lat=np.radians(24.67)
L=np.radians(83.5)
e=0.0167
k=4.47086*10**15
T_rot=23*3600+56*60
m_earth=5.972*10**24
m_sun=1.989*10**30
G=6.67408*10**-11
i=np.radians(23.5)

#Mean Solar Day Calculation
u=G*(m_sun+m_earth)
T_rev=2*np.pi*(k**3)/(u*u*(1-e*e)**1.5)
T_solar=1/(1/T_rot-1/T_rev)

#time frame
t=np.arange(0,T_rev,T_solar)

#Integration for anomaly
f=odeint(func,0,t,args=(e,k,u))
f=f.reshape(len(t))

#Ra, dec, Hour angle calculation
dec,ra=conv(i,0,np.pi-f)
ra=(ra-np.pi/2)%(2*np.pi)

H=(2*np.pi*t/T_rot+ra+L-np.pi/2)%(2*np.pi)

#Altitude and Azimuth
a,A=conv(lat-np.pi/2,dec,H)

#pLot
A[A>np.pi]=A[A>np.pi]-2*np.pi
Amin,Amax,amin,amax=np.degrees(np.min(A)),np.degrees(np.max(A)),np.degrees(np.min(a)),np.degrees(np.max(a))
m=Basemap(projection="stere",lon_0=(Amax+Amin)/2,lat_0=(amax+amin)/2,
                              width=200*10**5,height=120*10**5)
m.drawmapboundary(color='grey',linewidth=1.0,fill_color="black")
m.drawparallels(np.arange(-80.,81.,10.),color='white', textcolor='blue', 
                linewidth=1.0, labelstyle="+/-", labels=[1,1,1,0])
m.drawmeridians(np.arange(-180.,181.,10.),color='white', textcolor='green', 
                linewidth=1.0, labelstyle="+/-", labels=[0,0,0,1])

eqH=np.radians(np.arange(0.0,360.0,5))
eqa,eqA=conv(lat-np.pi/2,0,eqH)
eqA[eqA>np.pi]=eqA[eqA>np.pi]-2*np.pi

m.scatter(np.degrees(eqA),np.degrees(eqa),latlon=True,
          s=1,marker="o",color="white")

m.scatter(np.degrees(A),np.degrees(a),latlon=True,
          s=1,marker="o",color="red")

pl.show()