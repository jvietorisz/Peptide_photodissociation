import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#define constants

L = .5*10**(-3) #width of the laser beam
P =  .01
          # laser power in W
rho = .032/(np.pi*(.5*L)**2)#P/(np.pi*(L/2)**2) #laser power density in W/m^2

mu_IR = 10**(5)
mu_UV = 10**(0) #linear attenuation coefficient of water at 10.6 um in meters
r_0 = 20*10**(-9) #tip radius in m

theta = 3*(np.pi/180) #tip angle in radians

#discretize distance from tip
l = np.linspace(0, L, 1001) #makes a 1 by 501 array for the length of capillary tip

dl = l[1]-l[0]           #define step along capillary axis (each step is 1 um)

r = np.zeros(1001)

q_IR = np.zeros(1001)       #dimesions of heat generation array match those of l
q_UV = np.zeros(1001)

absorbed_power_IR = np.zeros(1001)
absorbed_power_UV = np.zeros(1001)

sum_IR = np.zeros(1001)     #initialize the sum of the power absorbed in each cross section
sum_UV = np.zeros(1001)

for i in range(0, len(l)):
    r[i] = r_0 +l[i]*np.tan(theta)

    beta = 0        #initialize the angle inside each section
    n = 0           #initialize the number of steps h taken to adjust beta

    while beta < np.pi:

        h = r[i]*(np.cos(beta+np.pi/180)-np.cos(beta))
        absorbed_power_IR[i] = (h*dl)*rho*(1-np.e**(-2*mu_IR*np.sin(beta)*r[i]))
        absorbed_power_UV[i] = (h*dl)*rho*(1-np.e**(-2*mu_UV*np.sin(beta)*r[i]))

        sum_IR[i] += absorbed_power_IR[i]          #add each cell to cross section sum
        sum_UV[i] += absorbed_power_UV[i]

        beta = beta+np.pi/180

    q_IR[i] = -sum_IR[i]/(dl*np.pi*r[i]**2)
    q_UV[i] = -sum_UV[i]/(dl*np.pi*r[i]**2)


print("q_IR is", q_IR)
#finding best fit curve
from scipy.optimize import curve_fit

#scale values for fitting and plotting (undone in later step)

xdata = (10**6)*l
ydata_IR = q_IR*(10**(-10))
ydata_UV = q_UV*(10**(-10))

def fit_func_exp(x, a, b, c):
    y = a*np.exp(-b*(x))+c
    return y

#fit function for IR
popt_IR, pccov_IR = curve_fit(fit_func_exp, xdata, ydata_IR, maxfev = 5000)
a,b,c = popt_IR
print("a = ", a, "b = ", b, "c = ", c)

popt_UV, pccov_UV = curve_fit(fit_func_exp, xdata, ydata_UV, maxfev = 10000)
d,e,f = popt_UV
print("d = ", d, "e = ", e, "f = ", f)

fit_IR = (10**(10))*fit_func_exp((10**6)*l, a, b, c)
fit_UV = (10**(10))*fit_func_exp((10**6)*l, d, e, f) #scaled back to original

#plot best fit and array

plt.plot(l, q_IR, color = "red", label= "IR", linewidth = 2)
plt.plot(l, fit_IR, color = "pink", linestyle ="dashed", label= "IR fit", linewidth = 2)
plt.plot(l, q_UV, color = "blue", label= "UV", linewidth = 2)
plt.plot(l, fit_UV, color = "cyan", linestyle ="dashed", label= "UV fit", linewidth = 2)

plt.legend()
plt.xlabel("distance from tip ($m$)", fontsize = 12)
plt.ylabel("absorbed power density ($W m^{-3}$)", fontsize = 12)
plt.yscale("log")
plt.grid()
#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.legend()
#matplotlib.pyplot.axis([0, .1, 40, 100])
plt.show()
