
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt
""" Dictionaire de variables

Fonctions:
-degarad(angle): Prend un angle comme argument et le transforme en degrés.
-calculdistance(rsoleil,anglerad): Calcule le rayon galactique à partir de
la distance entre le soleil et le centre galactique et de la longitude galactique.
-calculvitesse(rsoleil,vsoleil,distance,vitessemax): Calcule la vitesse de rotation à
à partir de la vitesse de rotation du soleil, sa distance au centre galactique,
les distances calculées et leurs vitesses.
-deltar(dr0,dl,r0,l): Calcule les incertitudes sur la distance.
-deltav(dvmax,v0,r0,dr,r,dv0,dr0): Calcule les incetitudes sur la vitesse.
-deltaadim(dx0,dx,x0,x): Calcule les incertitudes des variables adimensionnées à partir des
variables dimensionnées.
-calculmasse(rayon,vitesse): Calcule la masse à partir des distances et leurs Vitesses.
-deltam(dr,dv,r,v): Calcule les incertitudes sur la masse.

Constantes:
a: Distance du semi-grand axe galactique.
b: Rayon de la sphère de Plumer dans le modèle de Toomre.
r0: Distance du soleil au centre galactique.
v0: Vitesse de rotation du soleil par rapport au centre galactique.
dl: Incertitude sur les angles en radiants
dr0: Incertitude sur r0
dv0: Incertitude sur v0
dvmax: Incertitude sur les vitesses
kpc: kpc en mètres
km: km en mètres
G: Contante universelle de gravitation

Variables:

data: Stocke les données.
lrad: Longitude galactique en radians.
r: Rayon galactique.
dr: Incertitude sur r.
v: Vitesse de rotation.
dv: Incertitude sur v
rtheo: Rayon théorique.
vtheo: Vitesse théorique (du modèle statistique).
vadim: Vitesse adimensionée.
radim: Rayon adimensioné.
dvadim: Incertitudes sur la vitesse adimensionée.
dradim: Incertitudes sur le rayon adimensioné.
Masse: Masse calculée à partir de données.
dm: Incertitudes sur la masse.
Mtheo: Masse théorique (à partir du modèle statistique).
header: Tête du fichier de sortie.
dataout: Fichier de sortie.
"""

"""Fonctions"""
def degarad(angle):
    rad=angle*2*np.pi/360
    return np.array(rad)
def calculdistance(rsoleil,anglerad):
    r=rsoleil*np.sin(anglerad)
    return r
def calculvitesse(rsoleil,vsoleil,distance,vitessemax):
    vitesse=vitessemax+vsoleil*distance/rsoleil
    return vitesse
def deltar(dr0,dl,r0,l):
    dr=np.sqrt((np.sin(l)*dr0)**2+(r0*np.cos(l)*dl)**2)
    return dr
def deltav(dvmax,v0,r0,dr,r,dv0,dr0):
    dv=np.sqrt(dvmax**2+(v0/r0*dr)**2+((r/r0)*dv0)**2+(-v0*r/r0**2*dr0)**2)
    return dv
def deltaadim(dx0,dx,x0,x):
    delta=np.sqrt((1/x0*dx)**2+(-x/x0**2*dx0)**2)
    return delta
def calculmasse(rayon,vitesse):
    r=rayon*kpc
    v=vitesse*km
    M=(v)**2*(r**2+(a+b)**2)**(3/2)/(G*r**2)
    Ms=M/1.98847e30
    return Ms
def deltam(dr,dv,r,v):
    r=r*kpc
    v=v*km
    pv=2*v/(G*r**2)*(r**2+(a+b)**2)**(3/2)*dv
    pr=v**2/(G*r**3)*(3*(r**2+(a+b)**2)**(1/2)*r**2-2*(r**2+(a+b)**2)**(3/2))*dr
    dm=np.sqrt(pv**2+pr**2)/1.98847e30
    return dm
""" Constantes """
a=15#kpc
b=1#kpc
r0=8.5#Kpc
v0=220 #km/s
dl=degarad(0.5)#rads
dr0=0.05#kpc
dv0=0.5#km/s
dvmax=0.5#km/s
kpc=3.086e19#m
km=1e3#m
G=6.67430e-11#s.i

"""Traitement de données pour l'analyse des vitesses"""
#donnees.csv archive avec première colonne l, deuxième r et troisième vmax
data = genfromtxt('donnees.csv', delimiter = ';')

#On transforme l en radiants:
lrad = degarad(data[:,0])

#On calcule r à partir de l en radiants:
r = calculdistance(r0,lrad)

#On calcule les incertitudes sur r:
dr = deltar(dr0,dl,r0,lrad)

#On calcule v à partir de vmax
v=calculvitesse(r0,v0,r,data[:,1])

#On calcule les incertitudes sur v
dv=deltav(dvmax,v0,r0,dr,r,dv0,dr0)

#On introduit le modèle théorique calculé dans le logiciel R.
rtheo=np.linspace(2,9,100)
vtheo=-180.61+453.25*np.log(rtheo)-68.36*rtheo

""" Graphe des vitesses"""
plt.figure()
plt.scatter(r,v,marker="+",c="k",label="Vitesses de rotation")
plt.xlabel("r[kpc]")
plt.xticks(np.arange(1,9,step=0.5))
plt.ylabel("V[km/s]")
plt.legend()
plt.savefig("courbe_rotation.png",dpi=500)
plt.show()
plt.close()

""" Graphe Adimensioné """
#On adimensione les variable et calcule leurs incertitudes
radim=r/r0
vadim=v/v0
dradim=deltaadim(dr0,dr,r0,r)
dvadim=deltaadim(dv0,dv,v0,v)
#On introduit leur modèle théorique
rtheo=np.linspace(0.05,1.5,1000)
vtheo=3.5881+2.0602*np.log(rtheo)-2.6414*rtheo
#Graphe:
plt.figure()
plt.scatter(radim,vadim,marker="+",label="Données expérimentales")
plt.plot(rtheo,vtheo,"--",c="k",label="Modèle statistique")
plt.xlabel("$r/r_0$")
plt.ylabel("$v/v_0$")
plt.legend()
plt.savefig("Analyse_vitesse.png",dpi=500)
plt.show()

#On calcule masse pour les nuages de magellan et pour r0, v0
print(calculmasse(8.5,220))
print(deltam(0.05,0.5,8.5,220))
print(calculmasse(100,175))
print(deltam(0.05,0.5,100,175))

#On calcule la masse sur nos données
Masse=calculmasse(r,v)
#On calcule les incertitudes sur la masse
dM=deltam(dr,dv,r,v)
#On introduit le modèle théorique
rtheo=np.linspace(2,8.5,1000)
vtheo=-180.61+453.25*np.log(rtheo)-68.36*rtheo
Mtheo=calculmasse(rtheo,vtheo)
#Graphe:
plt.scatter(r,Masse,marker="+",label="Données expérimentales")
plt.plot(rtheo,Mtheo,"--",c="k",label="Modèle statistique")
plt.xlabel("$r[kpc]$")
plt.ylabel("$M/M_{\odot}$")
plt.legend(loc=4)
plt.savefig("Analyse_masse.png",dpi=500)
plt.show()


""" Création de fichier de sortie"""
dataout=[]
dataout=np.append(r,(dr,v,dv,radim,dradim,vadim,dvadim,Masse,dM))
dataout=np.reshape(dataout,(len(r),10),order='F')
dataout=np.append(data,dataout,axis=1)

header=np.array(['l',"vmax",'r',"dr",'v',"dv","radim","dradim","vadim","dvadim","M/Ms","dM/Ms"])
with open("dataout.csv","w") as f:
    writer=csv.writer(f)
    writer.writerow(header)
    writer.writerows(dataout)
f.close()
