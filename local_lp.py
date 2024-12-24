#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 20:53:44 2024

@author: yogehs
"""


import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15,'lines.markersize':14})

import pandas as pd

import os
import seaborn as sns

import matplotlib.pyplot as plt
 

from class_mech_gui import Filament_mech_gui
from constants import *
from utils import *
import warnings

#%%static LP loop 
plt.rcParams.update({'font.size': 18,'lines.markersize':14,'figure.dpi':100,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6})
  
gui_bool = False
#short files
folder_path_arr = S_files_arr[1:]+C_files_arr[:-1];del_arc_L_arr = [0.10]
#folder_path_arr = diff_vol_long;del_arc_L_arr = [0.10]
#folder_path_arr = S_files_arr[3:5];del_arc_L_arr = [0.10]

#long files
folder_path_arr  =  all_long_final+braided+long_ring_file+diff_vol_long+human_mrn_all+monomer_long;del_arc_L_arr = [0.10]


#folder_path_arr  =  monomer_long;del_arc_L_arr = [0.20]

f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]

df_all = pd.DataFrame()
df_all_1_2 = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])
    #for i in [0]:
        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i],num_segs=10)
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:
            
            for del_l in del_arc_L_arr:
                ex.roll_win_L_p_1_per_1rgn_cont(del_arc_l = del_l,del_arc_l_nm = 10.0)
                df_temp = ex.df_output
                df_temp['config_label'] = f_label_arr[i][0]
                df_temp_12,swap_bool = df_fil1_fil_2_all_L_p_s(df_temp)
                if swap_bool:
                    print(f_label_arr[i])
                    
                    df_temp['swapped_arc_len'] = df_temp['avg_arc_l_win_fil_2'] -0.01
                    df_temp['swapped_arc_len_nm'] =  df_temp['max_arc_l']*(df_temp['avg_arc_l_win_fil_2'] -0.01)

                else:
                    df_temp['swapped_arc_len'] = df_temp['avg_arc_l_win_fil_1'] 
                    df_temp['swapped_arc_len_nm'] =df_temp['max_arc_l']* df_temp['avg_arc_l_win_fil_1'] 

                df_all= pd.concat([df_all,df_temp])
                df_all_1_2= pd.concat([df_all_1_2,df_temp_12])
df_all = df_all.round({'swapped_arc_len': 2, 'avg_arc_l_win_fil_1': 2})
