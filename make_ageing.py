#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot individual panels for the figures in the resolved ageing paper

Created on Thu Jun  6 12:39:41 2019

@author: yago
"""

from __future__ import print_function, division

import glob
import numpy as np

from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt


califa = fits.open('data/Ageing_Database51_EML_MAG.fits')
sample_table = Table.read('data/sample_table.tex')
output_dir = '../paper/figures/'
fig_base_size = 2  # inches

# =============================================================================
# SDSS CONTOUR
# =============================================================================

corr=3.0
contour_file='data/fp_g_r_0.9_contour.fits'
fp_data=fits.open(contour_file)
fp=fp_data[0].data
bingr=fp_data[1].data
binew=fp_data[2].data
# binew = 10**binew - corr
level90=fp_data[0].header['LEVELS']
levels=[level90]

# =============================================================================
# 
# =============================================================================

califaid, g_sdss, r_sdss, ew_sdss = np.loadtxt('data/the_good_table_CALIFA_SDSS', usecols=(0, 8,9, 12), unpack=True)

califaid = np.array(califaid, dtype=int)
# ew_sdss = 10**ew_sdss

# =============================================================================
# COSMETIC 
# =============================================================================

ew_ticks_names=['-2.9', '-2', '0', '10', '100']
ew_ticks_positions=[np.log10(-2.9+corr), np.log10(-2.0+corr), np.log10(corr),
                    np.log10(10.0+corr), np.log10(100.0+corr)]

for filename in glob.glob('*/*csv'):
    panel = Table.read(filename)
    n_rows = np.max(panel['ROW'])
    n_cols = np.max(panel['COL'])
    print('\n {} ({}x{}):\n'.format(filename, n_rows, n_cols),
          *panel['GALAXY'])

    fig, axes = plt.subplots(n_rows, n_cols,
                             sharex=True, sharey=True, squeeze=False,
                             gridspec_kw={'wspace': 0, 'hspace': 0,
                                          'left': 0, 'right': 1,
                                          'top': 1, 'bottom': 0})
    fig.set_figheight(fig_base_size*n_rows)
    fig.set_figwidth(fig_base_size*n_cols)

    for galaxy in panel:
        ax = axes[n_rows-galaxy['ROW'], n_cols-galaxy['COL']]
        galaxy_name = galaxy['GALAXY']

        table_entry = sample_table[sample_table['NAME'] == galaxy_name]
        
        califa_object = table_entry['CALIFAID']
        sdss_galaxy = califaid == califa_object
                
        ax.annotate(table_entry['CALIFAID'].pformat(show_name=False)[0],
                    xy=(.1, .1),  xycoords='axes fraction')
#        ax.annotate(table_entry[table_entry.colnames[:3]].pformat(show_name=False)[0],
#                    xy=(.5, .1),  xycoords='axes fraction',
#                    horizontalalignment='center', verticalalignment='bottom')
#        ax.annotate(table_entry[table_entry.colnames[3:]].pformat(show_name=False)[0],
#                    xy=(.5, .09),  xycoords='axes fraction',
#                    horizontalalignment='center', verticalalignment='top')
        ax.tick_params(axis='both', direction='in',
                       right=True, top=True,
                       labelbottom=False, labelleft=False)
        ax.axhline(np.log10(0+corr), c='k', ls=':')

        califa_data = califa[galaxy_name].data
        r = califa_data['Obs_AB_r']
        g = califa_data['Obs_AB_g']
        EW_Ha_emission = califa_data['Ha_flux'] / califa_data['Ha_baseline']
        EW_Ha_absorption = califa_data['HaSynAbsMod_flux'] / califa_data['HaSynAbsMod_baseline']
        EW_Ha = EW_Ha_emission + EW_Ha_absorption
        good_SN = (califa_data['Ha_baseline'] > 30*califa_data['e_Ha_baseline'])

        sc = ax.scatter((g-r)[good_SN], np.log10(EW_Ha[good_SN]+corr), marker='s', s=0.1,
                        c=r[good_SN], cmap='gist_rainbow', vmin=20.5, vmax=22.5)
#                        c=r, cmap='gist_stern', vmin=18, vmax=24)
#                        c=r, cmap='gist_heat', vmin=18, vmax=24)

        # SDSS 90% contour
        
        ax.contour(bingr, binew, fp, levels, colors='gray')
        
        # SDSS values 
        ax.scatter(g_sdss[sdss_galaxy]-r_sdss[sdss_galaxy], ew_sdss[sdss_galaxy], s=35,
                   marker='^', c='gray', edgecolors='k')
        
        ax.set_xlim(.05, 1.15)
        ax.set_ylim(-1.3, 2.3)
        ax.set_yticks(ticks=ew_ticks_positions)
        ax.set_yticklabels(labels=ew_ticks_names)
        
        # ax.set_yscale('symlog', linthresh=3)
#        fig.colorbar()

    for ax in axes.flatten():
        if len(ax.get_children()) == 10:  # empty plot
            ax.axis('off')

#    fig.tight_layout()
    fig.savefig(output_dir+'ageing/'+filename[:-3]+'pdf')
    # plt.close()

plt.figure(figsize=(fig_base_size*1.2, fig_base_size*1.1))
plt.xlim(.05, 1.15)
plt.ylim(-1.3, 2.3)
plt.yscale('symlog')
plt.xlabel('g - r')
plt.ylabel(r'EW H$\alpha$')
plt.axhline(np.log10(0 + corr), c='k', ls=':')
plt.yticks(ticks=ew_ticks_positions, labels=ew_ticks_names)
cbar = plt.colorbar(sc)
cbar.ax.set_title(r'$\mu_r$')
plt.tight_layout()
plt.savefig(output_dir+'ageing/legend.pdf')
plt.close()

# <codecell> Bye
# -----------------------------------------------------------------------------
#                                                           ... Paranoy@ Rulz!
# -----------------------------------------------------------------------------
