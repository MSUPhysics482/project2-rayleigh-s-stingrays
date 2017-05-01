import matplotlib.pyplot as plt 
import numpy as np
import os

file=np.loadtxt('data.dat').transpose()

wavelength=file[0]
unscattered=file[3]
scattered=file[4]

h = 6.626E-34      
c = 2.998E8
kb = 1.381E-23
T = 5778.0
x = np.arange(400.0,702.0,1.0)
y = 9.64**20.0/(x**5.0*(np.exp(h*c/(x*10.0**-9*kb*T))-1.0))

plt.figure(1)
ax = plt.axes()
plt.setp( ax.get_yticklabels(), visible=False)
plt.xticks(fontname="Arial")
plt.xlim(400,702)
plt.ylim(2000,15000)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.title("Solar Radiation Spectrum Before Rayleigh Scattering")

plt.bar(wavelength, scattered, 2.0, color='#FFD85B', bottom=unscattered, linewidth=0.0, label='Scattered Light')
plt.bar(wavelength, unscattered, 2.0, color='#FFA026', linewidth=0.0, label='Unscattered Light')
plt.plot(x, y, color='#000000', linewidth=1.0)

ax.annotate(r'$\propto \frac{1}{\lambda^5\left( e^{\frac{hc}{\lambda k_B T}}-1 \right)}$', xy=(620,10000), xytext=(620,12700), size=20, arrowprops=dict(arrowstyle='->'))

plt.legend(loc=4)

plt.savefig('plot1.pdf', format='pdf')
plt.savefig('plot1.png', format='png')
os.system('okular plot1.pdf')
