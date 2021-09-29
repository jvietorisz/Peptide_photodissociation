import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.io

mat = scipy.io.loadmat("/Users/Jacob/Desktop/Python/temp_distribution.mat")
xdata_mat = mat["X"]
ydata_mat = mat["Tintrp"]

#define constants

pow_den = np.array([10000, 30000, 40000, 50000,
                    90000, 200000, 300000, 500000,
                    600000, 900000, 2000000, 4000000,
                    6000000, 8000000, 10000000])
max_temp_IR = np.array([52.44, 107, 135, 162.2,
                        271.9, 573.7, 848, 1397,
                        1671.2, 2494, 5.5122*10**(3), 1.0999*10**(4),
                        1.6487*10**(4), 2.197*10**(4), 2.7461*10**(4)])

max_temp_193 = np.array([25.007, 25.0209, 25.0278, 25.0348,
                        25.0626, 25.1391, 25.2087, 25.3478,
                        25.4173, 25.6260, 26.3910, 27.7820,
                        29.1730, 30.5640, 31.9550])

max_temp_222 = np.array([25.0007, 25.0021, 25.0028, 25.0035,
                         25.0063, 25.0139, 25.0209, 25.0348,
                         25.0417, 25.0626, 25.1391, 25.2782,
                         25.4173, 25.5565, 25.6956])

#convert maximum temperature to temp difference
temp_diff_IR = max_temp_IR - 25
temp_diff_193 = max_temp_193 - 25
temp_diff_222 = max_temp_222 - 25

#finding best fit curve
from scipy.optimize import curve_fit

def fit_func(x, a, b):
    y = a*x+b
    return y

popt_IR, pccov_IR = curve_fit(fit_func, pow_den, temp_diff_IR, p0 = [1, 0], maxfev = 10000)
a,b = popt_IR
print("a = ", a, "b = ", b)

popt_193, pccov_193 = curve_fit(fit_func, pow_den, temp_diff_193, p0 = [1, 0], maxfev = 10000)
c,d = popt_193
print("c = ", c, "d = ", d)

popt_222, pccov_222 = curve_fit(fit_func, pow_den, temp_diff_222, p0 = [1, 0], maxfev = 10000)
e,f = popt_222
print("e = ", e, "f = ", f)

#create x array for power density

x_values = np.linspace(5000, 10000000, 1000)

#create y array for fitted max temp
fit_IR = fit_func(x_values, a, b)
fit_193 = fit_func(x_values, c, d)
fit_222 = fit_func(x_values, e, f)



plt.plot(x_values, fit_IR, color = "red", label= "10.6 $\mu m$", linewidth = 2)
plt.plot(x_values, fit_193, color = "blue", linestyle ="dashed", label= "193 nm", linewidth = 2)
plt.plot(x_values, fit_222, color = "darkviolet", linestyle ="dotted", label= "222 nm", linewidth = 2)
plt.scatter(pow_den, temp_diff_IR, color = "hotpink")
plt.scatter(pow_den, temp_diff_193, color = "deepskyblue")
plt.scatter(pow_den, temp_diff_222, color = "mediumorchid")
#plt.legend()
plt.yscale("log")
plt.xscale("log")

plt.grid()
#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.legend()
#matplotlib.pyplot.axis([0, .1, 40, 100])
plt.show()
