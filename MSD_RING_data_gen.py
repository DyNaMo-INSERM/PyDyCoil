#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:16:40 2024

@author: yogehs
"""



import chime 

import matplotlib.pyplot as plt

import pandas as pd

import os
import seaborn as sns
 
from class_mech_gui import Filament_mech_gui
from constants import *
from utils import *
from viz_long_files_MSD_dyanmics import *
import warnings
warnings.filterwarnings("ignore")
import sys

plt.rcParams.update({'font.size': 15,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6})

gui_bool = False

del_arc_L_arr = [0.50]
y_min , y_max = 1,100
y_min_norm , y_max_norm = 0.00008,0.8
op_save = '/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/'
#%%lop for S and C long files
#save_path = folder_path+"mech_analysis_gui.csv"

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan']
 

folder_path_arr = all_long_final
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]


df_SC_raw_all = pd.DataFrame()

df_linear_fit = pd.DataFrame()
df_SC_all = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr,df_raw = ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                #c = c_arr[i_file]
      
                #sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_SC_all = pd.concat([df_SC_all,df_msd])
                df_SC_raw_all = pd.concat([df_SC_raw_all,df_raw])



    
df_SC_all.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_SC_all_.csv')
df_SC_raw_all.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_SC_all_.csv')

#%%lop forring files
#save_path = folder_path+"mech_analysis_gui.csv"

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan']
 

folder_path_arr = long_ring_file
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]




df_linear_fit = pd.DataFrame()
df_ring = pd.DataFrame();df_ring_raw = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr,df_raw = ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                c = c_arr[i_file]
      
                sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_ring = pd.concat([df_ring,df_msd])
                df_ring_raw = pd.concat([df_ring_raw,df_raw])



    
df_ring.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_ring_.csv')

df_ring_raw.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_ring_.csv')
#%%C plots  
#R
R = sns.relplot(df_ring,x='deltaT_s',y='MSD',style= 'fil_label',hue='file_label' ,s=100,legend = False)  
R.set(title = 'C ',xlabel = None,xscale="log", yscale="log",box_aspect = 1)
R = sns.relplot(df_ring,x='deltaT_s',y='MSD',style= 'fil_label',hue='file_label',col = 'file_label' ,s=100,legend = False)  
R.set(xlabel = None,xscale="log", yscale="log",box_aspect = 1)

#R norm
R = sns.relplot(df_ring,x='dT_MSD_norm_s',y='MSD_norm',style= 'fil_label',hue='file_label' ,s=100,legend = False) 
#ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
 
R.set(title = 'C  ',xlabel = None,xscale="log", yscale="log",box_aspect = 1,ylim=(y_min_norm,y_max_norm))
R = sns.relplot(df_ring,x='dT_MSD_norm_s',y='MSD_norm',style= 'fil_label',hue='file_label',col = 'file_label' ,s=100,legend = False)  
R.set(xlabel = None,xscale="log", yscale="log",box_aspect = 1)

# tan angle 
R = sns.relplot(df_ring,x='deltaT_s',y='tan_ang_MSD',style= 'fil_label',hue='file_label' ,s=100,legend = False)  
R.set(xlabel = None,title = 'C  ',xscale="log", yscale="log",box_aspect = 1)
R = sns.relplot(df_ring,x='deltaT_s',y='tan_ang_MSD',style= 'fil_label',hue='file_label',col='file_label' ,s=100,legend = False)  
R.set(xlabel = None,xscale="log", yscale="log",box_aspect = 1)
R.fig.subplots_adjust(wspace=0.1, hspace=0)

sns.displot(data=df_ring, x="L_seg", col='file_label')


#%%lop for MRN files
#save_path = folder_path+"mech_analysis_gui.csv"

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan','tab:cyan','tab:cyan','tab:cyan','tab:cyan']
 

folder_path_arr =   human_mrn_all
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]


df_MRN_raw = pd.DataFrame()


df_linear_fit = pd.DataFrame()
df_MRN = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr, df_raw= ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                c = c_arr[i_file]
      
                sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_MRN = pd.concat([df_MRN,df_msd])
            
                df_MRN_raw = pd.concat([df_MRN_raw,df_raw])


    
df_MRN.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_MRN_.csv')
df_MRN_raw.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_MRN_.csv')

#%%lop for braided files
#save_path = folder_path+"mech_analysis_gui.csv"

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan']
 

folder_path_arr = braided
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]



df_braid_raw = pd.DataFrame()

df_linear_fit = pd.DataFrame()
df_braid = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr,df_raw = ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                c = c_arr[i_file]
      
                sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_braid = pd.concat([df_braid,df_msd])
                df_braid_raw= pd.concat([df_braid_raw,df_raw])


df_braid.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_braid_.csv')
df_braid_raw.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_braid_.csv')

#%%lop for monomer files
#save_path = folder_path+"mech_analysis_gui.csv"

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan']
 

folder_path_arr = monomer_long
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]

del_arc_L_arr = [0.98]


df_mono_raw = pd.DataFrame()

df_linear_fit = pd.DataFrame()
df_mono = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr,df_raw = ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                c = c_arr[i_file]
      
                sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_mono = pd.concat([df_mono,df_msd])
            
                df_mono_raw = pd.concat([df_mono_raw,df_raw])
                

df_mono_raw.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_mono_.csv')

df_mono.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_mono_.csv')
#%%lop for diff vols C files

c_arr = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:cyan']
 

folder_path_arr = diff_vol_long
f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]

df_diff_vol_raw = pd.DataFrame()



df_linear_fit = pd.DataFrame()
df_diff_vol = pd.DataFrame()
if __name__ == '__main__': 
    for i in range(len(folder_path_arr)):
        print(folder_path_arr[i])

        ex = Filament_mech_gui(folder_path_arr[i],f_label_arr[i])
        if gui_bool:
            
            app = QApplication(sys.argv)
            ex.load_gui()
            sys.exit(app.exec_())
        else:

            
            for del_l in del_arc_L_arr:
                df_msd,df_msd_arr,df_raw = ex.R_TAN_BEND_MSD(del_l)
                plt.figure(1)
                slope_bend = 1.0/4.0
                slope_R = 3.0/4.0
                X_plot = df_msd['dT_MSD_norm'].to_numpy()
                #Y_bend = df_msd['tan_ang_MSD'].to_numpy()[0]*X_plot**slope_bend 
                Y_R = df_msd['MSD_norm'].to_numpy()[5]*X_plot**slope_R 

                title_temp = df_msd['file_label'].to_numpy()[0]
                
                i_file = i
                c = c_arr[i_file]
      
                sns.scatterplot(df_msd,x='deltaT_s',y='MSD',style= 'fil_label',color =c ,s=100,legend = False)#,ax =ax_R[i_file]) 
                #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
               
                
                df_diff_vol = pd.concat([df_diff_vol,df_msd])
                df_diff_vol_raw = pd.concat([df_diff_vol_raw,df_raw])

df_diff_vol_raw.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/raw_long_diff_vol_.csv')
df_diff_vol.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_diff_vol_.csv')

#%%
df_raw_all = pd.concat(            [df_SC_raw_all,df_ring_raw,df_MRN_raw,df_braid_raw,df_mono_raw,df_diff_vol_raw])

df_raw_all.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/RAW_ALL_MSD.csv')

for name,df_temp in df_raw_all.groupby(by = ['file_label']):
    print(name)
    plt.title(name[0])
    sns.scatterplot(df_temp , x = 'time_s',y = 'end_end_dist_2',hue = 'fil_label')
    plt.legend(frameon = False)
    plt.show()
    