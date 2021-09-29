import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

aminoAcids = np.array(["Tyr", "Phe", "Trp", "His",
                      "Met", "Arg", "Lys", "Val",
                      "Thr", "Leu", "Ser", "Ala",
                      "Pro", "Hyp", "Gly", "Asp"])

sigma = np.array([16.44157406*10**(-21), 10.74437747*10**(-21), 7.64724375*10**(-21), 1.950047156*10**(-21),
                 0.8794330312*10**(-21), 4.2824565*10**(-21), 0.4205984062*10**(-21), 0.09941416874*10**(-21),
                 0.1338267656*10**(-21), 0.09559054687*10**(-21), 0.1108850344*10**(-21), 0.1070614125*10**(-21),
                 0.1108850344*10**(-21), 0.1070614125*10**(-21), 0.09176692499*10**(-21), 0.09176692499*10**(-21)])

Phi = np.array([9.3*10**(-2), 5.7*10**(-2), 6.7*10**(-2), 26*10**(-2),
               26.2*10**(-2), 1.8*10**(-2), 14.5*10**(-2), 14.5*10**(-2),
               10.5*10**(-2), 13.9*10**(-2), 10.3*10**(-2), 1.8*10**(-2),
               6.7*10**(-2), 5.1*10**(-2), 5.9*10**(-2), 2*10**(-2)])


PowDen = 10000#lamp power density in W/m^2

t = np.linspace(0, 6, 1000) #time in seconds

PeV = PowDen / (1.6022*10**(-19)) # lamp power density in eV/sm^2
J = PeV / (6.42) #lamp photon flux in photons/sm^2

for i in range(0, len(aminoAcids)):
    print(aminoAcids[i])
    print(1/(J*sigma[i]*Phi[i])) #for obtaining the characteristic times of decomposition of the amino acids

#define sigma and phi for peptide bond

sigma_pb_193 = 2.485354219*10**(-21)
sigma_pb_222 = 1.4*10**(-22)
Phi_pb = 5.9*10**(-2)

#print(1/(J*sigma_pb_193*Phi_pb)) #for obtaining the characteristic time of decomposition of the peptide bond

Tyr_array = np.zeros(len(t)) #initialize y-values for Trp (curve with lowest peak probability)

for i in range(0, len(aminoAcids)):

    #total probability
    plt.plot(t, ((1-np.e**(-t*J*sigma_pb_193*Phi_pb))**2)*(np.e**(-t*J*sigma[i]*Phi[i])), label = "{}".format(aminoAcids[i]), linewidth = 1)

    #AA survival
    #plt.plot(t, np.e**(-t*J*sigma[i]*Phi[i]), label = "{}".format(aminoAcids[i]), linewidth = 1)

for i in range(0, len(t)):
    Tyr_array[i] = ((1-np.e**(-t[i]*J*sigma_pb_193*Phi_pb))**2)*(np.e**(-t[i]*J*sigma[3]*Phi[3]))

print(aminoAcids[3])
print("Tyr max prob is", np.amax(Tyr_array))

plt.xlabel("Time (hrs)", fontsize = 12)
plt.ylabel("Probability", fontsize = 12)

plt.yticks(np.arange(0, 1.2, step=.2))
plt.grid()
#plt.legend(frameon = False, bbox_to_anchor=(1,1), loc="upper left")
plt.show()
