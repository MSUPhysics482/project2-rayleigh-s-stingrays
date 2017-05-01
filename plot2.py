import matplotlib.pyplot as plt 
import numpy as np
import os

file=np.loadtxt('data.dat').transpose()

wavelength=file[0]
after=file[2]

h = 6.626E-34      
c = 2.998E8
kb = 1.381E-23
T = 5778.0
x = np.arange(400.0,702.0,1.0)
y = 9.64**20.0/(x**5.0*(np.exp(h*c/(x*10.0**-9*kb*T))-1.0))

plt.figure(1)
ax = plt.axes()
plt.setp( ax.get_yticklabels(), visible=False)
plt.xlim(400,702)
plt.ylim(0,15000)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.title("Solar Radiation Spectrum After Rayleigh Scattering")

plt.bar(wavelength, after, 2.0, color="#6379A0", linewidth=0.0, label="After")
plt.plot(x, y, color='#000000', linewidth=1.0)

plt.legend(loc=4)

plt.savefig('plot2.pdf', format='pdf')
plt.savefig('plot2.png', format='png')
os.system('okular plot2.pdf')
