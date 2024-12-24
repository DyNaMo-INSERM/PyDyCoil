#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:33:25 2024

@author: yogehs
"""

k_bT_nm = 4.114
mode_roll = 'reflect'

S_files_arr = ['/Users/yogehs/Downloads/SbcC/S/1/','/Users/yogehs/Downloads/SbcC/S/2/','/Users/yogehs/Downloads/SbcC/S/3/',
 '/Users/yogehs/Downloads/SbcC/S/4/','/Users/yogehs/Downloads/SbcC/S/5/','/Users/yogehs/Downloads/SbcC/S/6/']
trans_files_arr = ['/Users/yogehs/Downloads/SbcC/S_to_C/1/','/Users/yogehs/Downloads/SbcC/V_to_S/1/',
                   ]
C_files_arr = [     '/Users/yogehs/Downloads/SbcC/C/1/','/Users/yogehs/Downloads/SbcC/C/2/','/Users/yogehs/Downloads/SbcC/C/3/',
                 '/Users/yogehs/Downloads/SbcC/C/4/','/Users/yogehs/Downloads/SbcC/C/5/',
                 '/Users/yogehs/Downloads/SbcC/C/6/','/Users/yogehs/Downloads/SbcC/I /1/']

S_long = ['/Users/yogehs/Downloads/SbcC/S/1_long/','/Users/yogehs/Downloads/SbcC/S/2_long/',
             '/Users/yogehs/Downloads/SbcC/S/3_long/']
S_long_new = ['/Users/yogehs/Downloads/SbcC/S/4_long/','/Users/yogehs/Downloads/SbcC/S/5_long/',
             '/Users/yogehs/Downloads/SbcC/S/6_long/']
C_long = ['/Users/yogehs/Downloads/SbcC/C/1_long/',
'/Users/yogehs/Downloads/SbcC/C/2_long/','/Users/yogehs/Downloads/SbcC/C/3_long/','/Users/yogehs/Downloads/SbcC/C/4_long/',
'/Users/yogehs/Downloads/SbcC/C/5_long/','/Users/yogehs/Downloads/SbcC/C/6_long/']




mixed_files = ['/Users/yogehs/Downloads/SbcC/mixed/1_long/','/Users/yogehs/Downloads/SbcC/mixed/2_long/']

C_diff_fs = ['/Users/yogehs/Downloads/SbcC/C/diff_fs_C/0.5 fps/',
           '/Users/yogehs/Downloads/SbcC/C/diff_fs_C/1.3 fps/']

S_2_end_fix = ['/Users/yogehs/Downloads/SbcC/S/S_2_end_fixed/monomer_1/','/Users/yogehs/Downloads/SbcC/S/S_2_end_fixed/monomer_2/']




human_mrn_all = ['/Users/yogehs/Downloads/SbcC/HMRN_H_closed/1/'
                 ,'/Users/yogehs/Downloads/SbcC/HMRN_H_closed/2/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_closed/3/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_closed/4/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_open/1/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_open/2/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_open/3/',
                 '/Users/yogehs/Downloads/SbcC/HMRN_H_open/4/'
                 ]

braided = ['/Users/yogehs/Downloads/SbcC/braided/1/',
           '/Users/yogehs/Downloads/SbcC/braided/2/',
           '/Users/yogehs/Downloads/SbcC/braided/3/']
monomer_long =['/Users/yogehs/Downloads/SbcC/monomer/1/'] 
long_ring_file = ['/Users/yogehs/Downloads/SbcC/ring_noATP/1/',
                  '/Users/yogehs/Downloads/SbcC/ring_noATP/2/',
                  '/Users/yogehs/Downloads/SbcC/ring_noATP/3/',
                  '/Users/yogehs/Downloads/SbcC/ring_ATP/1/'
                  ]

diff_vol_long = ['/Users/yogehs/Downloads/SbcC/vol_diff_S/1/',
                 '/Users/yogehs/Downloads/SbcC/vol_diff_S/2/',
                 '/Users/yogehs/Downloads/SbcC/vol_diff_S/3/',
                 '/Users/yogehs/Downloads/SbcC/vol_diff_S/4/']

long_vids_fps = {"C_1_long":1,"C_2_long":1,"C_3_long":1,"C_4_long":1,"C_5_long":3.33,"C_6_long":0.5,
                 "S_1_long":1,"S_2_long":1,"S_3_long":1,"S_4_long":3.33,"S_5_long":2,"S_6_long":2,
                 'ring_noATP_1':2,'ring_noATP_2':2,'ring_noATP_3':2,'ring_ATP_1':1,
                 'mixed_1_long':0.5,'mixed_2_long':1.3,'HMRN_H_closed_1':1,
                  'HMRN_H_closed_2':1,
                  'HMRN_H_closed_3':1,
                  'HMRN_H_closed_4':1,
                  'HMRN_H_open_1':1,
                  'HMRN_H_open_2':1,
                  'HMRN_H_open_3':1,
                  'HMRN_H_open_4':1,
                  'braided_1':1, 'braided_2':1, 
                  'braided_3':1,'monomer_1':1,'monomer_2':1,
                  'vol_diff_S_1':1,'vol_diff_S_2':1,'vol_diff_S_3':1,'vol_diff_S_4':1}   


all_long_final = S_long+S_long_new+C_long #+mixed_files
color_saha_configs_old ={'S':'red',"C":'darkorange',"closed":'forestgreen',"open":'magenta','b':'dodgerblue','m':'darkcyan','r':'darkviolet','v':'greenyellow'}

color_saha_configs ={'S':'orangered',"C":'peru',"closed":'forestgreen',"open":'magenta','b':'dodgerblue','m':'darkcyan','r':'darkviolet','v':'mediumseagreen'}
cmap_saha_configs ={'S':'Reds',"C":'Oranges',"closed":'Greens','open':'Greys','b':'Blues','m':'Greys','r':'Purples','v':'YlGn'}

final_save_plt_path = '/Users/yogehs/Downloads/SbcC/plots_paper/final_plots_svg_20241108/'