#integer function Rand_Poission(lambda)
#double precision, intent(in) :: lambda
# double precision             :: lambda_left, p, u
# int k

import numpy as np
import matplotlib.pyplot as plt
import random

def poission_dist(lamb): #testput = 0.2
    lambda_left = lamb
    k = 0.0
    p = 1.0
    Poisson_Step = 500

    while (p >= 1):
        k = k + 1
        u = np.random.rand()
        p = p*u
        while ((p < 1) and (lambda_left > 0)):
            if (lambda_left > Poisson_Step):
                p = p * np.exp(Poisson_Step)
                lambda_left = lambda_left - Poisson_Step
            else:
                p = p * np.exp(lambda_left)
                lambda_left = 0

    Rand_Poission = k - 1
    return Rand_Poission
    #print(Rand_Poission)

#  ! where step is the current time step
# ! returns the number of electrons allowed to be emitted in that time step
#  double precision function Gauss_Emission(step)
#step = range(0,5000)
#sigma = 300.0
#mu = 2000.0
#A = 5.0
#b = 1.0/(2.0*np.pi*sigma**2)

#    ! Gauss_Emission = IDNINT(  A * exp( -1.0d0*b*(step - mu)**2 )  )
#for i in step:
#    Gauss_Emission = A * np.exp( -1.0 * b * (i - mu)**2 )
#    Gauss_Emission_stak = np.floor(Gauss_Emission)
#    x = poission_dist(Gauss_Emission)
#    y = poission_dist(Gauss_Emission_stak)
#    plt.plot(i,x,'k.', ms=6)
#    plt.plot(i,y,'b.', ms=6)
#    plt.plot(i,Gauss_Emission,'r.', ms=6)
    

  #double precision function Photon_Emission(photon_energy)
runs = range(1,1000)
sigma = 1.0 #! Width / standard deviation
mu = 4700 #! Center
A = 1.0 #! Height
b = 1.0/(2.0*np.pi*sigma**2)
photon_energy = 4.7
freq_var = 0.01
hist_list = []
for i in runs:
    Photon_power = np.floor(photon_energy*10000)
    Frequency = np.floor(freq_var*10000)
    rand_photon = np.random.rand()
    photon_rand = (Photon_power-Frequency) + np.floor(((Photon_power+Frequency)-(Photon_power-Frequency))*rand_photon)
    Photon_Emission = poission_dist(photon_rand)/10000
    #photon_rand = random.randrange(1670,7730)
    #photon = photon_rand/1000
    #u = np.random.rand()
    #photon_rand = 4670 + np.floor((4730+1-4670)*u)
    #Photon_Emission = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (photon_rand - mu)**2 / (2 * sigma**2))
    #Photon_Emission_pois = poission_dist(photon_rand)
    #photon = Photon_Emission_pois/1000
    hist_list.append(Photon_Emission)
    #print(Photon_Emission)
    #plt.plot(photon,Photon_Emission_pois,'b.', ms=6)

n, bins, patches = plt.hist(hist_list, 50, density=True, facecolor='g', alpha=0.75)
plt.title('Histogram of Photon energies')
#plt.text(photon_energy, freq_var, r'$\mu=100,\ \sigma=15$')
plt.xlabel('Photon energy [eV]')
plt.ylabel('Nr of Photons')
#plt.xlim(4, 5)
#plt.ylim(0, 200)
plt.grid(True)
plt.show()      