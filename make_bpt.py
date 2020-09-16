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


califa = fits.open('Ageing_Database51_EML_MAG.fits')
sample_table = Table.read('sample_table.tex')
output_dir = '../paper/figures/'
fig_base_size = 2  # inches

# Kewwley et al. (2006)
N2_K06 = np.arange(-2, 0.05, .1)
O3_K06 = 1.3 + 0.61/(N2_K06-0.05)

for filename in glob.glob('mass_envir/*csv'):
    panel = Table.read(filename)
    n_rows = np.max(panel['ROW'])
    n_cols = np.max(panel['COL'])
    print('{} ({}x{}):'.format(filename, n_rows, n_cols),
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

        ax.plot(N2_K06, O3_K06, 'k:')

        califa_data = califa[galaxy_name].data
        r = califa_data['Obs_AB_r']
        Ha = califa_data['Ha_flux']
        Hb = califa_data['Hb_flux']
        OIII = califa_data['OIII_5007A_flux']
        NII = califa_data['NII_6584A_flux']
        good_SN = (califa_data['Hb_flux'] > 3*califa_data['e_Hb_flux'])

        sc = ax.scatter(np.log10(NII/Ha)[good_SN], np.log10(OIII/Hb)[good_SN],
                        marker='s', s=.1,
                        c=r[good_SN], cmap='gist_rainbow', vmin=20.5, vmax=22.5)
#                        c=r, cmap='gist_stern', vmin=18, vmax=24)
#                        c=r, cmap='gist_heat', vmin=18, vmax=24)
        ax.set_xlim(-1.5, 1.25)
        ax.set_ylim(-1.5, 1.25)
        ax.set_yscale('symlog')
#        fig.colorbar()

    for ax in axes.flatten():
        if len(ax.get_children()) == 10:  # empty plot
            ax.axis('off')

#    fig.tight_layout()
    fig.savefig(output_dir+'BPT/'+filename[:-3]+'pdf')
    plt.close()

plt.figure(figsize=(fig_base_size*1.2, fig_base_size))
plt.xlim(-1.5, 1.25)
plt.ylim(-1.5, 1.25)
plt.xlabel(r'log$_{10}$(NII/H$\alpha$)')
plt.ylabel(r'log$_{10}$(OIII/H$\beta$)')
plt.plot(N2_K06, O3_K06, 'k:')
cbar = plt.colorbar(sc)
cbar.ax.locator_params(nbins=5)
cbar.ax.set_title(r'$\mu_r$')
plt.tight_layout()
plt.savefig(output_dir+'BPT/legend.pdf')
plt.close()

# <codecell> Bye
# -----------------------------------------------------------------------------
#                                                           ... Paranoy@ Rulz!
# -----------------------------------------------------------------------------
