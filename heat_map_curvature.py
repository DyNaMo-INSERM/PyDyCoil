#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:56:25 2024

@author: yogehs

curvature heat maps  
#TODO get the same for the LP
"""


import seaborn as sns

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 50,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6,
                     'font.family':'Times New Roman','axes.linewidth':2.5,
                     'xtick.major.width': 2,'ytick.major.width': 2,'xtick.minor.width': 2,'ytick.minor.width': 2,'xtick.minor.visible':True})

import os 
import pandas as pd 
from class_fil import Filament
import numpy as np 
from constants import *
from utils import *
import warnings
warnings.filterwarnings("ignore")
fit_E_header = ['config_label','avg_arc_L','curv_2_coeff','curv_1_coeff','curv_0_coeff','GOF_Chi'] 
from constants import final_save_plt_path

def drop_0s(arr,curv_arr):
    where_no_0 = np.nonzero(arr)
    del_E = arr[where_no_0]-np.min(arr[where_no_0])

    return(del_E,curv_arr[where_no_0])

def do_energy_landscape(df_all_frame):
    file_label = df_all_frame['config_label'].to_numpy()[0]
    norm = plt.Normalize(-5,0)

    norm_arc_l_final = df_all_frame['norm_arc_L'].to_numpy()
    pt_curv_final = df_all_frame['pt_curv'].to_numpy()
    H, xedges, yedges = np.histogram2d(norm_arc_l_final,pt_curv_final , bins = 20,density=True)
    #energy_landscape = -np.log(H)
    H = H.T
    where_0 = np.where(H == 0)
    where_1 = np.where(H == 1)

    H[where_0] = 1
    energy_landscape = -np.log(H)
    df_fit = pd.DataFrame()
    for i in range(len(energy_landscape)):
        E,curv= drop_0s(energy_landscape[:,i],yedges)
        
        x_normalized = (curv - curv.mean())/curv.std()
        fit_normalized, residuals, _, _, _  = np.polyfit(x_normalized, E, 2,full=True)
        chisq_dof = residuals / (len(x_normalized) - 3)
        f = np.poly1d(fit_normalized)
        avg_arc_L = (xedges[i+1]+xedges[i])*0.5
        print(fit_normalized)
        '''
        plt.title(f"{file_label} norm arc l : {xedges[i]}")
        plt.scatter(x_normalized,E)
        
        plt.plot(x_normalized,f(x_normalized),c ='black' ,label=f'k_{fit_normalized[0]}')
        plt.xlim(-1.75,1.75)        
        plt.legend()plt.show()
        '''
        temp_arr = [file_label,avg_arc_L,fit_normalized[0],fit_normalized[1],fit_normalized[2],chisq_dof]        
        df_dict = pd.DataFrame([dict(zip(fit_E_header, temp_arr))])
        df_fit = pd.concat([df_fit,df_dict], ignore_index=True)

        
    #E_min = np.min(energy_landscape[where_1]);E_max = np.max(energy_landscape[where_1])
    fig = plt.figure(figsize=(20, 20));plt.clf()
    ax = fig.add_subplot(131, title=f'{file_label} curvature dist',aspect='equal')
    plt.imshow(H, interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    #plt.ylim(-0.5,0.5)
    norm = plt.Normalize(-5,0)

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(131, title=f'{file_label}_energy landscaoe')
    
    cax = plt.imshow(energy_landscape, interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap='inferno_r',norm=norm)
    plt.colorbar(cax,shrink=0.3)
    
    yedges = np.round(yedges,3)
    fig = plt.figure(100,figsize=(20, 20))
    ax = fig.add_subplot(131, title=f'{file_label}_energy landscaoe')
    
    ht_sns = sns.heatmap(energy_landscape,cbar_kws={'shrink': 0.3},square = True,norm = norm)
    plt.xticks([]);plt.yticks(ticks = [0,20],labels =[yedges[0],yedges[-1]],rotation=45, fontsize=12);plt.gca().invert_yaxis()
    
    return df_fit
#%% short files 
folder_path_arr = C_files_arr[:-1]#
folder_path_arr = S_files_arr#+C_long
#folder_path_arr  = C_long[:-1]
#folder_path_arr  = S_long

save_path = '/Users/yogehs/Downloads/SbcC/plots_paper/overlays_L_p/'


win_i =0.125

lin_arr = np.linspace(0,1,int(np.ceil(1/(win_i))), False)
l_p_arr = np.ones((2,len(lin_arr)))

f_label_arr = [find_folder_name(os.path.dirname(f))+"_"+find_folder_name(f) for f in folder_path_arr]
max_arc_l_all_files = np.zeros(len(folder_path_arr))
from class_mech_gui import Filament_mech_gui
df_final = pd.DataFrame()
flip_arr = [False,False,False]
if __name__ == '__main__':  
    for k in range(len(folder_path_arr)):
        
        ex = Filament_mech_gui(folder_path_arr[k],f_label_arr[k])#,10,flip_corrd =flip_arr[k])
        temp_dict= {}
        print(ex.num_frames,k)
        N = 40
        df_all_frame= pd.DataFrame()
        for i in range(ex.num_frames):
            ex.i_frame = i
            i_filament = ex.filaments_arr[i]
            temp_dict["pt_curv"] = i_filament.pt_curv_ys

            temp_dict["frame_num"] = i
            temp_dict["config_label"] = f_label_arr[k][0]

            temp_dict["file_label"] = f_label_arr[k]
            temp_dict["norm_arc_L"] = i_filament.arc_l_norm


            df_temp = pd.DataFrame.from_dict(temp_dict)
            df_all_frame = pd.concat([df_all_frame,df_temp])
            df_final = pd.concat([df_final,df_temp]) 
        #do_energy_landscape(df_all_frame)
df_E_fit  = do_energy_landscape(df_final)
plt.show()
sns.scatterplot(df_E_fit,x= 'avg_arc_L',y='curv_2_coeff')
plt.ylim(0,4)
plt.show()
arc_l = df_E_fit['avg_arc_L'].to_numpy()
flex = df_E_fit['curv_2_coeff'].to_numpy()
plt.title(f_label_arr[k][0])
ratio_l1_l2 = []
for i in range(int(len(arc_l)/2)):
    plt.scatter(arc_l[i], flex[i],color = 'tab:blue')
    plt.scatter(1-arc_l[-(i+1)], flex[-(i+1)],color = 'tab:orange')
    temp_r = flex[i]/flex[-(i+1)]
    if temp_r>0:
        
        ratio_l1_l2.append(temp_r)
plt.ylim(0,4)
plt.show()

#%%joint plot 

# Create a scatter plot with KDE heatmap
plt.figure(figsize=(10, 6))
sns.kdeplot(df_final, x='norm_arc_L', y='pt_curv', cmap="viridis", fill=True, thresh=0, levels=100)
#sns.scatterplot(df_final, x='norm_arc_L', y='pt_curv', color='white', s=10, alpha=0.6)
plt.ylim(-0.5,0.5)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Scatter Plot with KDE Heatmap')
plt.show()
#%%
from matplotlib.colors import Normalize
plt.rcParams.update({'font.size': 40,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6,
                     'font.family':'Times New Roman','axes.linewidth':2.5,
                     'xtick.major.width': 2,'ytick.major.width': 2,'xtick.minor.width': 2,'ytick.minor.width': 2,'xtick.minor.visible':True})
color_i = color_saha_configs[f_label_arr[k][0]]
plt.figure(figsize=(20,15))
temp = sns.jointplot(data=df_final, x="norm_arc_L", y="pt_curv", kind="hist"
            ,color = color_i)
cax = temp.fig.add_axes([0.86, 0.2, 0.03, 0.6])  # [left, bottom, width, height]
temp.ax_joint.set_ylim(-0.5, 0.5)
#temp.ylim(-0.5,0.5)
norm = Normalize(vmin=0, vmax=1)
temp.ax_joint.collections[0].set_clim(0,100) 
# Add the color bar
cbar = plt.colorbar(temp.ax_joint.collections[0], cax=cax,norm = norm)


cbar.set_ticks([0, 100])  # Set tick positions
cbar.set_ticklabels(['0', '100'])  # Set custom tick labels

#temp.ax_marg_x.clear()
#temp.ax_marg_x.set_axis_off()
#temp.ax_marg_y.set_yticks([])

# Adjust tick parameters for thicker ticks
temp.ax_joint.tick_params(axis='both', which='major', length=8, width=2,labelsize = 18)  # Main plot ticks
temp.ax_marg_y.tick_params(axis='y', which='major', length=8, width=2,labelsize = 18)    # Y-marginal ticks

# Adjust the linewidth of the main plot and marginal plots
temp.ax_joint.collections[0].set_linewidth(2)  # Adjust linewidth of hexbin edges
for spine in temp.ax_joint.spines.values():
    spine.set_linewidth(2)  # Adjust the linewidth of the plot border

temp.ax_marg_x.spines['bottom'].set_linewidth(2)
temp.ax_marg_y.spines['left'].set_linewidth(2)
temp.fig.suptitle(f_label_arr[k][0])
temp.fig.set_dpi(300)

plt.savefig(f"{final_save_plt_path}curvature_joint_plot_for_{f_label_arr[k][0]}.svg")




#%% 2 d histogram 

#%% 

norm_arc_l_final = df_all_frame['norm_arc_L'].to_numpy()
pt_curv_final = df_all_frame['pt_curv'].to_numpy()
H, xedges, yedges = np.histogram2d(norm_arc_l_final,pt_curv_final , bins = 10,density=True)
#energy_landscape = -np.log(H)
H = H.T
where_0 = np.where(H == 0)
where_1 = np.where(H == 1)

H[where_0] = 1
energy_landscape = -np.log(H)

for i in range(len(energy_landscape)):
    E,curv= drop_0s(energy_landscape[:,i],yedges)
    
    x_normalized = (curv - curv.mean())/curv.std()
    fit_normalized = np.polyfit(x_normalized, E, 2)
    
    f = np.poly1d(fit_normalized)
    plt.title(f"norm arc l : {xedges[i]}")
    plt.scatter(curv,E)
    
    plt.plot(curv,f(x_normalized), c='black')

    plt.xlim(yedges[0],yedges[-1])
    #plt.ylim(-5,-1)
    plt.show()
    

    
#E_min = np.min(energy_landscape[where_1]);E_max = np.max(energy_landscape[where_1])
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(131, title='imshow: square bins',aspect='equal')
plt.imshow(H, interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
#plt.ylim(-0.5,0.5)

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(131, title='energy landscaoe')
sns.heatmap(energy_landscape,cbar_kws={'shrink': 0.3},square = True,norm = norm)
plt.xticks([]);plt.yticks([])
#plt.imshow(energy_landscape, interpolation='nearest', origin='lower')#,extent=[xedges[0], xedges[-1]])
#plt.legend()
  #%%
df_S = df_final[df_final['config_label']=='S']
df_C = df_final[df_final['config_label']=='C']

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Create the first joint plot
g1 = sns.JointGrid(x="norm_arc_L", y="pt_curv", data=df_S)
g1.plot(sns.scatterplot, sns.histplot)
axes[0].set_title('Joint Plot 1: Total Bill vs Tip')

# Create the second joint plot
g2 = sns.JointGrid(x="norm_arc_L", y="pt_curv", data=df_C)
g2.plot(sns.scatterplot, sns.histplot)

# Adjust the layout
plt.tight_layout()
plt.show()
      
