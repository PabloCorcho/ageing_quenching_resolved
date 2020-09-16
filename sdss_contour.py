# -*- encoding: utf-8 -*-
import numpy as np # esta es la librería de python que trabaja con arrays
import matplotlib.pyplot as plt #plotting library
from astropy.io import fits as pyfits

#-------------------------------------------------------------
#-------------------------------------------------------------

# DATA for fp contour
#data to plot the fp plot from paper SDSS
files='../../sdss_data/'
contour_file=files+'fp_g_r_0.9_contour.fits'
print('Opening file: '+contour_file)
fp_data=pyfits.open(contour_file)
fp=fp_data[0].data
bingr=fp_data[1].data
binew=fp_data[2].data
level90=fp_data[0].header['LEVELS']

#-------------------------------------------------------------
#-------------------------------------------------------------
#THE NEEDED CORRECTION TO PLOT IT PROPERLY
corr=3.0

# EW ticks and ticklabels
ew_ticks_names=['-2.9','-2.5', '-2','-1', '0','10','50','100']
ew_ticks_positions=[np.log10(-2.9+corr),np.log10(-2.5+corr),np.log10(-2.0+corr),np.log10(-1.0+corr),np.log10(corr),np.log10(10.0+corr),np.log10(50.0+corr),np.log10(100.0+corr)]
# g-r ticks and ticklabels
gr_ticks_names=['.0', '.2', '.4', '.6', '.8', '1', '1.2']
gr_ticks_positions=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
#plt.yticks(ew_ticks_positions,ew_ticks_names)

# fp axis limits
mingr=0.0
maxgr=1.2
minew=-1.3
maxew=2.0

plt.figure()

plt.yticks(ew_ticks_positions,ew_ticks_names)
plt.xticks(gr_ticks_positions,gr_ticks_names)
plt.ylabel('EW (H$_\\alpha$) - $\AA$', fontsize=6)
plt.xlabel('(g-r) color', fontsize=6)

levels=[level90]
plt.contour(bingr, binew, fp, levels, colors='gray')

plt.show()
