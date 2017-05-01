import pylab
import os

# TIME-AVG POYNTING VECTOR #

poynting = r'$\ \ \ \langle \mathbf{S} \rangle = \frac{\mu_0 \pi^2 q_0^2 d^2 c^3 \sin^2 \theta}{2r^2 \lambda^4} \hat{\mathbf{r}} $'
fig = pylab.figure()
text = fig.text(0, 0, poynting)
dpi = 600
fig.savefig('poynting.png', dpi=dpi)
bbox = text.get_window_extent()
width, height = bbox.size / float(dpi) + 0.2
fig.set_size_inches((width, height))
dy = 2*(bbox.ymin/float(dpi))/height
text.set_position((0, -dy))
fig.savefig('poynting.png', dpi=dpi)

# SCATTERING CROSS-SECTION #

cross_section = r'$ \ \ \sigma (\lambda) = \frac{2\pi^5 d^6}{3 \lambda^4}  \left( \frac{n^2-1}{n^2+2} \right)^2 = \frac{\alpha}{\lambda^4}$'
fig = pylab.figure()
text = fig.text(0, 0, cross_section)
fig.savefig('cross_section.png', dpi=dpi)
bbox = text.get_window_extent()
width, height = bbox.size / float(dpi) + 0.2
fig.set_size_inches((width, height))
dy = 2*(bbox.ymin/float(dpi))/height
text.set_position((0, -dy))
fig.savefig('cross_section.png', dpi=dpi)

# NUMBER DENSITY #

density = r'$ \ \ \ n(h) = n_0 exp \left( - \frac{Mgh}{RT_0} \right)$'
fig = pylab.figure()
text = fig.text(0, 0, density)
fig.savefig('density.png', dpi=dpi)
bbox = text.get_window_extent()
width, height = bbox.size / float(dpi) + 0.2
fig.set_size_inches((width, height))
dy = 2*(bbox.ymin/float(dpi))/height
text.set_position((0, -dy))
fig.savefig('density.png', dpi=dpi)

# PROBABILITY #

probability = r'$ \ \ \ P(\lambda) = \int_0^\infty n(h)\sigma(\lambda)dh = \frac{\alpha}{\lambda^4} \frac{n_0RT_0}{Mg}$'
fig = pylab.figure()
text = fig.text(0, 0, probability)
fig.savefig('probability.png', dpi=dpi)
bbox = text.get_window_extent()
width, height = bbox.size / float(dpi) + 0.2
fig.set_size_inches((width, height))
dy = 1.5*(bbox.ymin/float(dpi))/height
text.set_position((0, -dy))
fig.savefig('probability.png', dpi=dpi)


os.system('okular poynting.png cross_section.png density.png probability.png')