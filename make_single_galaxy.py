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

# %%
def running_median(array, map_array, neightbours=20):
    sort = np.argsort(map_array)
    sorted_array = array[sort]    
    running_median = np.zeros_like(sorted_array)
    running_median[:] = np.nan
    running_std = np.zeros_like(sorted_array)
    running_std[:] = np.nan
    
    for i in range(len(sorted_array)):   
        if (i>=neightbours/2)&(i<=len(sorted_array)-neightbours/2):
            running_median[i] = np.nanmedian(
                            sorted_array[int(i-neightbours/2):int(i+neightbours/2)])    
            running_std[i] = np.nanstd(
                            sorted_array[int(i-neightbours/2):int(i+neightbours/2)])    
        
    return running_median, running_std
    
# %%
for filename in glob.glob('morph_mass/*csv'):
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
        ax.axhline(0, c='k', ls=':')

        califa_data = califa[galaxy_name].data
        r = califa_data['Obs_AB_r']
        g = califa_data['Obs_AB_g']
        EW_Ha_emission = califa_data['Ha_flux'] / califa_data['Ha_baseline']
        EW_Ha_absorption = califa_data['HaSynAbsMod_flux'] / califa_data['HaSynAbsMod_baseline']
        EW_Ha = EW_Ha_emission + EW_Ha_absorption
        good_SN = (califa_data['Ha_baseline'] > 30*califa_data['e_Ha_baseline'])

        median, std = running_median(EW_Ha[good_SN].flatten(), 
                                     r[good_SN].flatten())
        # sc = ax.scatter((g-r)[good_SN], EW_Ha[good_SN], marker='s', s=0.1,
        #                 c=r[good_SN], cmap='gist_rainbow', vmin=20.5, vmax=22.5)
#                        c=r, cmap='gist_stern', vmin=18, vmax=24)
#                        c=r, cmap='gist_heat', vmin=18, vmax=24)
        
        sc = ax.scatter(r[good_SN], EW_Ha[good_SN].flatten(), c=(g-r)[good_SN],
                   s=0.1, cmap='jet', vmin=0.5, vmax=1)
        ax.plot(np.sort(r[good_SN].flatten()), median, 'k')         
        ax.set_xlim(18, 23)
        ax.set_ylim(-20, 200)
        ax.set_yscale('symlog')
#        fig.colorbar()

    for ax in axes.flatten():
        if len(ax.get_children()) == 10:  # empty plot
            ax.axis('off')

#    fig.tight_layout()
    fig.savefig(output_dir+'ageing_p/'+filename[:-3]+'pdf')
    plt.close()

plt.figure(figsize=(fig_base_size*1.2, fig_base_size*1.1))
plt.xlim(18, 23)
plt.ylim(-20, 200)
plt.yscale('symlog')
plt.xlabel('r')
plt.ylabel(r'EW H$\alpha$')
plt.axhline(0, c='k', ls=':')
plt.yticks(ticks=[-10, -1, 0, 1, 10, 100], labels=['-10', '-1', '0', '1', '10', '100'])
cbar = plt.colorbar(sc)
cbar.ax.set_title(r'$g-r$')
plt.tight_layout()
plt.savefig(output_dir+'ageing_p/legend.pdf')
plt.close()

# <codecell> Bye
# -----------------------------------------------------------------------------
#                                                           ... Paranoy@ Rulz!
# -----------------------------------------------------------------------------
