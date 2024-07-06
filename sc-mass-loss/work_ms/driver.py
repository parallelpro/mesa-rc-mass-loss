from __future__ import print_function
import numpy as np
import os
import sys
import h5py
import pandas as pd


# read in input params
params = pd.read_csv('grid_input_params.txt')

# get from input which track to run for this simulation
index_this_run = int(sys.argv[1])
idx = (params['index'] >= index_this_run) & (params['index'] < (index_this_run+1))

# for track in tracks[idx]:
for ip, p in params.loc[idx,:].iterrows():

    # get input params for this simulation
    # index, massinit, Xinit, Yinit, Zinit, amlt = p['index'], p['massinit'], p['Xinit'], p['Yinit'], p['Zinit'], p['amlt']
    # fehinit, fov_core, fov0_core, eta = p['fehinit'], p['fov_core'], p['fov0_core'], p['eta']
    index = p['index']

    output_history_name = f'index{index:06.0f}.history'
    output_final_model_name = f'index{index:06.0f}_botrgb.mod'

    print('Now calculating ', output_final_model_name)

    # # # Step 1: modify inlist 
    params = {'save_model_filename': output_final_model_name,
              'star_history_name': output_history_name,
              'initial_mass': p['massinit'],
              'initial_z': p['Zinit'],
              'initial_y': p['Yinit'],
              'mixing_length_alpha': p['amlt'],
              'overshoot_f': p['fov_core'],
              'overshoot_f0': p['fov0_core'],
              'sc_scaling_factor': p['eta'],
              }
    
    with open('inlist_template', 'r') as file:
        inlist_template = file.read()

    inlist = inlist_template.format(**params)

    with open('inlist_run', 'w') as file:
        file.write(inlist)

    # # # Step 2: run MESA.
    # os.system('\\rm -r LOGS; \\rm -r png; \\rm -r photos')

    print('------ MESA start ------')
    os.system('sh rn1 > mesa_terminal_output_index{:06.0f}.txt'.format(index))
    print('------ MESA done ------')


    filepath = 'LOGS/'

    # # # Step 3: create a .h5 file to store history and frequencies. Only run if history file exists.
    if os.path.exists(filepath+output_history_name):

        # # read in models
        track = pd.read_fwf(filepath+output_history_name, skiprows=5, infer_nrows=10000)

        # # append grid initial parameters as new columns
        for col in p.index:
            track[col] = p[col]

        # # append log properties
        track['luminosity'] = 10.0**track['log_L']
        track['radius'] = 10.0**track['log_R']
        track['Teff'] = 10.0**track['log_Teff']

        # # append seismic scaling quantities
        Dnu_sun, numax_sun, Teff_sun = 135.1, 3090., 5772.
        track['Dnu_int'] = track['delta_nu']
        track.drop(columns=['delta_nu'])
        track['Dnu_scaling'] = track['star_mass']**0.5 * track['radius']**-1.5 * Dnu_sun
        track['numax'] = track['star_mass'] * track['radius']**-2.0 * (track['Teff']/Teff_sun)**-0.5 * numax_sun
        
        # # append surface quantities
        Zsun, Xsun = 0.0134, 0.7381 # 0.0134, 0.7381, a09 # 0.0169, 0.7345, gs98
        track['mh'] = np.log10((1-track['surface_h1']-track['surface_he4']-track['surface_he3'])/track['surface_h1']) - np.log10(Zsun/Xsun)

        sumPaths = [filepath+f for f in os.listdir(filepath) if f.endswith('.txt')]
        sumDirs = np.array([f.split('/index')[0]+'/' for f in sumPaths])
        sumNames = np.array([f.split('/')[-1] for f in sumPaths])


        # # read in radial mode frequencies
        seismicCols = ['l', 'n_p', 'n_g', 'n_pg', 'E_p', 'E_g', 'E_norm', 'Re(freq)', 'Re(freq)_corr']
        # seismicCols = ['l', 'n_p', 'n_g', 'n_pg', 'E_norm', 'freq']
        seismicData = [[] for i in range(len(seismicCols))]

        for it, t in track.loc[:,:].iterrows():
            if t['flag_gyre']:
                # sumFile = filepath+'index{:06.0f}profile{:0.0f}.data.FGONG.sum'.format(index, profileIndex)
                sumFile = f"index{p['index']:06.0f}.history.model{t['model_number']:06.0f}.txt"
                if len(sumDirs[sumNames==sumFile])==0:
                    for i in range(len(seismicData)):
                        seismicData[i].append(np.nan)
                    track.loc[it, 'flag_gyre'] = 0
                else:
                    sumPath = filepath + sumFile
                    sum = pd.read_fwf(sumPath)
                    for i in range(len(seismicData)):
                        seismicData[i].append( sum[seismicCols[i]].to_numpy() )
            else:
                for i in range(len(seismicData)):
                    seismicData[i].append(np.nan)


        # #  write out the table
        with h5py.File(f"index{p['index']:06.0f}_msrg.history.h5", 'w') as h5f:
            for col in track.columns:
                h5f.create_dataset(col, data=track[col].to_numpy())

            for it, t in track.iterrows():
                if t['flag_gyre']: 
                    for i in range(len(seismicCols)):
                        h5f.create_dataset('model_number{:0.0f}/{:s}'.format(t['model_number'], seismicCols[i]), data=seismicData[i][it])
                else:
                    continue

