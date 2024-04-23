import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

d = np.array([0.148,0.165,0.205,0.235,0.260,0.285,0.310,0.330,0.370,0.405,0.405,0.520,0.570,0.620])
r = np.array([2.2,3.5,4.4,5.8,6.4,8.8,9.6,14.2,15.9,17.4,25.6,32.7,33.1,44.4])/2.2046
logd = np.log(d)
logr = np.log(r)
slope, intercept, r_value, p_value, std_err = linregress(logd, logr)

x = np.linspace(0.130, 0.700, 200)
logx = np.log(x)
logy = slope*logx + intercept
"""
plt.figure(figsize=(12,8))
plt.loglog(d, r, 'k.', label='résistance selon d')
plt.minorticks_on()
plt.grid(True, which='both', linestyle=':', linewidth='0.5')
plt.legend()
plt.show()
"""
print(f'k = {slope}')
print(f'a = {np.exp(intercept)}')


d2 = np.array([0.172,0.194,0.222,0.248,0.272,0.294,0.314,0.416,0.444])
r2 = np.array([3.10,3.82,4.78,6.04,7.06,7.90,8.94,15.2,17.16])
logd2 = np.log(d2)
logr2 = np.log(r2)
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(logd2, logr2)

x2 = np.linspace(0.160, 0.500, 200)
logx2 = np.log(x2)
logy2 = slope2*logx2 + intercept2

plt.figure(figsize=(12,8))
plt.plot(logd2, logr2, 'k.', label='log(résistance) selon log(d)')
plt.plot(logx2, logy2, 'r-', label='linear fit in loglog')
plt.plot(logd, logr, 'k.', label='log(résistance) selon log(d)')
plt.plot(logx, logy, 'r-', label='linear fit in loglog')
plt.grid(True, linestyle=':', linewidth='0.5')
plt.legend()
plt.show()

print(f'k = {slope2}')
print(f'a = {np.exp(intercept2)}')