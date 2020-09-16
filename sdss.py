from astropy.io import fits as pyfits
from astropy.table import Table

sdss = Table.read('SDSS_EML_MAG_51.fits')

#VALUES FOR THE OBJECT - SDSS - OBTAINED THROUGH RUBEN'S PIPELINE
Haf_em_sdss_r=sdss['Ha_Flux']
Hac_em_sdss_r=sdss['Ha_baseline']
Haf_abs_sdss_r=sdss['HaSynAbsMod_Flux']# Integrated - Voigt
Hac_abs_sdss_r=sdss['HaSynAbsMod_baseline']
names_r=sdss['Name']

g_sdss_r2=sdss['Obs_AB_g']
r_sdss_r2=sdss['Obs_AB_r']

g_r_sdss_r2=g_sdss_r2-r_sdss_r2

ew_em_r=Haf_em_sdss_r/Hac_em_sdss_r
ew_abs_r=(Haf_abs_sdss_r/Hac_abs_sdss_r)
ew_r=ew_em_r+ew_abs_r

