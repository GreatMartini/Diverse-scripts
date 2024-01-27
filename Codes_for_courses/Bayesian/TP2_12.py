import numpy as np
import numpy.fft as fft
from scipy.stats import norm
from scipy.stats import chi2
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from getdist import plots, MCSamples

psd = np.load("./simu_cluster/psd.npy")
S = len(psd)
temps=np.linspace(0,10,S)

def bruit_colore():
    bruit = np.random.normal(0,1E-3,S)
    bruit_fourier = fft.fft(bruit)
    bruit_fourier2 = bruit_fourier*np.sqrt(psd)
    return np.real(fft.ifft(bruit_fourier2))



def covmat():
    N=10000
    cov = np.zeros((S,S))
    i = 0
    for i in range(N):
        vecteur = bruit_colore()[:]
        produit = vecteur[None,:] * vecteur[:,None]
        cov += produit
        i += 1
    return cov/N

inverse_covmat = np.linalg.inv(covmat())

r = np.load("./simu_cluster/r.npy")
y = np.load("./simu_cluster/y.npy")

def density_profil(r,theta):
    rf =  r  / theta[1]
    rho = (theta[0] / (((rf)**0.31)*((1 + (rf)**1.1))**((5.5 - 0.31)/1.1))) + theta[2] * norm.pdf(r, loc=theta[3], scale=theta[4])
    return rho

def khi2(rho):
    diff = (y - rho)
    return np.dot(diff[None,:],np.dot(inverse_covmat,diff[:,None]))

#########################################################################


def theta_prime(theta):
    dev=[1e-3,50,0.2,50,10]
    t=np.random.normal(theta,dev)
    print(t)
    return t
#Donne les priors sur theta ou theta prime
def log_prior(theta):
    #Conditions
    conditions=[False,False,False,False,False]
    #Conditions sur rho 0
    if (10>theta[0]>0):
        conditions[0]=True
    #Conditions sur rp
    if (theta[1]>0):
        conditions[1]=True
    #Conditions sur A:
    if (theta[2]>0):
        conditions[2]=True
    #Conditions sur mu
    if (500<=theta[3]<=5000):
        conditions[3]=True
    #Conditions sur sigma:
    if (10<=theta[4]<=1000):
        conditions[4]=True
    if (all(conditions)==True):
        return 0
    else:
        return -np.inf

def log_likelihood(r, theta):
    rho = density_profil(r,theta)
    diff = (y - rho)
    return -1/2*np.dot(diff[None,:],np.dot(inverse_covmat,diff[:,None]))

def log_posterior(r,theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(r, theta)

def log_acceptance(theta, thetaprime):
    A = log_posterior(r,thetaprime)
    #B = priorprime
    C = log_posterior(r,theta)
    #D = priori
    return A-C

def test(theta, thetaprime, logalpha):
    logu = np.log(np.random.uniform(0,1))
    if (logu <= logalpha):
        return thetaprime
    if (logu > logalpha):
        return theta


#Faire fonction log_posterieur, après log prior, revoyer -np.inf si le prior est -np.inf
#et sinon calculer le log acceptance
def metropolis(theta_init,N):
    theta=theta_init
    for i in range(N):
        thetaprime = theta_prime(theta)
        priorprime = log_prior(thetaprime)
        prior = log_prior(theta)
        logacc = log_acceptance(theta,thetaprime)
        theta = test(theta, thetaprime, logacc)
        if(i%1000==0):
            plt.plot(r,density_profil(r,theta), label="courbe_"+str(i))
    plt.legend()
    plt.show()
    return theta

thetai = [0.08, 200, 30, 1500, 300]
metropolis(thetai,10000)
