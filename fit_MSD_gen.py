#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 16:25:36 2024

@author: yogehs
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import seaborn as sns 
import pandas as pd 
from constants import color_saha_configs,cmap_saha_configs,final_save_plt_path
kbT = 4.114 #pn nm

y_min , y_max = 1,1000
 
x_renorm_limits = [0.030644186650635584, 4996.008916478943]
y_renorm_limit = [0.06738070130927563, 2.9694807353946584]

x_un_limits = [0.23042308789504504, 78.19536741138113]
y_un_limits = [0.7854736945816706, 148.95038821847717]
y_min_norm , y_max_norm = 0.00008,0.8
fit_l_p_header = ['file_label','config_label','fil_label','avg_arc_L','L_fitted','L_P_fitted'
                  ,'zeta_fitted',"avg_L_p_WLC",'tau_c','r_c'] 
from statistics import geometric_mean
def tau_sum(L,lp,zeta,max_n = 15):
    sum_term = 0
    for n in range(1, max_n + 1):
        tau_n = (zeta / (2 * kbT * lp)) * (L / (n * np.pi)) ** 4
        
        #sum_term += (1 - np.exp(-t / tau_n)) / (n * np.pi) ** 4

# Define the function to be fitted
def delta_R(t, zeta,lp,L , max_n=15):
    sum_term = 0
    for n in range(1, max_n + 1):
        tau_n = (zeta / (2 * kbT * lp)) * (L / (n * np.pi)) ** 4
        sum_term += (1 - np.exp(-t / tau_n)) / (n * np.pi) ** 4
    return (L**4 / lp**2) * sum_term
def delta_R_fit(t,zeta,lp ,  L):
    return delta_R(t, zeta,lp, L)

def fit_func_delta_L(df_all,skip_pt =12):
    '''
    
    this function fits for L_p and zeta for 100 tersm the L is passed as from the iamge  
    
    Parameters
    ----------
    df_all : data frame of file 1 and 2 for a given cofig
    Returns
    -------
    df_fit : TYPE
        DESCRIPTION.

    '''
    df_fit = pd.DataFrame()
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    for i in range(2):
        fl = f"fil{i+1}"
        mask= df_all['fil_label']==fl
        df = df_all[mask]
        if len(df)>1:
    
            t_data = df['deltaT_s'].to_numpy()
            
            MSD_R = df['MSD'].to_numpy()
            
            file_label = df['file_label'].to_numpy()[0]
            fil_label = df['fil_label'].to_numpy()[0]
            i_fil = int(fil_label[-1])-1
            avg_arc_L = df['avg_arc_L'].to_numpy()[0]
            L_seg = df['L_seg'].to_numpy()[0]
            
             
            
            
            mean_l_p = df['mean_l_p'].to_numpy()[0]
            
            L_p_arr_forward = df['l_p_i'].to_numpy()
            L_p_rev_arr = df['l_p_rev_i'].to_numpy()
            L_p_arr = np.concatenate((L_p_arr_forward, L_p_rev_arr))
            L_p_mean = np.mean(L_p_arr)
            marker_arr = ['o','x'];color_arr = ['tab:blue','tab:orange']
            initial_guess = [ 0.01,60]  # Initial guesses for L, lp, zeta, kappa
            
            if skip_pt+3 >len(t_data):
                skip_pt = -2
                
            params, covariance = curve_fit(lambda t,zeta,lp   : delta_R_fit(t,zeta, lp,  L_seg),
                       t_data[:-skip_pt], MSD_R[:-skip_pt], p0=initial_guess)

            zeta_fitted,lp_fitted = params;L_fitted = L_seg
            
            '''
            initial_guess = [ 0.01]  # Initial guesses for L, lp, zeta, kappa
            
            params, covariance = curve_fit(lambda t,zeta: delta_R_fit(t,zeta, L_p_mean,  L_seg),
                       t_data[3:-10], MSD_R[3:-10], p0=initial_guess)

            zeta_fitted = params[0];L_fitted = L_seg;lp_fitted = L_p_mean
            '''
            tau_c_temp= 1.0/ ((4.73/L_fitted)**4 * (kbT*lp_fitted/zeta_fitted))

            r_c_temp = 1.0/ ((90*lp_fitted**2)* (1.0/L_fitted)**4 )


            temp_arr = [file_label,file_label[0],fil_label,avg_arc_L,L_fitted,lp_fitted,zeta_fitted,mean_l_p,tau_c_temp,r_c_temp]
            axs[0].set_title(f"{file_label}_MSD")
            axs[0].scatter(t_data, MSD_R, marker = marker_arr[i_fil], color=color_arr[i])
            axs[0].plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
            axs[0].set_xlabel('Time (s)');axs[0].set_ylabel('MSD  R')
            axs[0].set_xscale('log');axs[0].set_yscale('log')
            axs[0].set_ylim(y_min,y_max);axs[0].set_xlim(0,80);

            
            axs[1].set_title(f"{file_label}_WLC")
            axs[1].hist(L_p_arr,bins = 10, color=color_arr[i],histtype = 'step',density = False,label = f'{len(L_p_arr)}')
            axs[1].axvline(lp_fitted,linestyle = '--',color = color_arr[i])
            axs[1].axvline(np.mean(L_p_arr),linestyle = '-',color = color_arr[i])

            axs[1].set_xlabel('L_P(t)_dist')#;axs[i].set_ylabel('MSD  R')
            axs[1].set_xlim(5,125)
            df_dict = pd.DataFrame([dict(zip(fit_l_p_header, temp_arr))])

            df_fit = pd.concat([df_fit,df_dict], ignore_index=True)

        axs[0].legend(frameon = False);        axs[1].legend(frameon = False)

    return df_fit

def viz_all_MSD(df,label):
    
    R = sns.relplot(df,x='deltaT_s',y='MSD',style= 'fil_label',hue='file_label' ,s=100,legend = False)  
    R.set(title = label,xlabel = None,xscale="log", yscale="log",box_aspect = 1,ylim=(y_min,y_max),xlim =(0,80))
    R = sns.relplot(df,x='deltaT_s',y='MSD',style= 'fil_label',hue='avg_arc_L',col = 'file_label' ,s=100,legend = True)  
    R.set(xscale="log", yscale="log",box_aspect = 1,ylim=(y_min,y_max),xlim =(0,80))

    #R norm
    R = sns.relplot(df,x='dT_MSD_norm_s',y='MSD_norm',style= 'fil_label',hue='file_label' ,s=100,legend = False) 
    #ax_R[i_file].plot(X_plot, Y_R, '--r',label = f'{slope_R:.2f} line')
     
    R.set(title = label,xlabel = None,xscale="log", yscale="log",box_aspect = 1,ylim=(y_min_norm,y_max_norm))
    R = sns.relplot(df,x='dT_MSD_norm_s',y='MSD_norm',style= 'fil_label',hue='file_label',col = 'file_label' ,s=100,legend = False)  
    R.set(xscale="log", yscale="log",box_aspect = 1)

    # tan angle 
    R = sns.relplot(df,x='deltaT_s',y='tan_ang_MSD',style= 'fil_label',hue='file_label' ,s=100,legend = False)  
    R.set(xlabel = None,title =label,xscale="log", yscale="log",box_aspect = 1,ylim = (0.01,0.75),xlim =(0,80))
    R = sns.relplot(df,x='deltaT_s',y='tan_ang_MSD',style= 'fil_label',hue='file_label',col='file_label' ,s=100,legend = False)  
    R.set(xscale="log", yscale="log",box_aspect = 1,ylim = (0.01,0.75),xlim =(0,80))
    R.fig.subplots_adjust(wspace=0.1, hspace=0)
#%%2 segments
#df_all_2 = pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/S_C_ALL_MSD_2_segs_raw.csv')
file_to_ignore = ['C_4_long','S_4_long','mixed_2_long','mixed_1_long']
df_all_2 = pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_SC_all_.csv')
df_ring = pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_ring_.csv')
df_mrn = pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_MRN_.csv')
df_braid= pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_braid_.csv')
df_mono= pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_mono_.csv')
df_diff_vol= pd.read_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/long_diff_vol_.csv')

df_all_2 = pd.concat([df_all_2,df_mrn,df_braid,df_ring,df_mono,df_diff_vol])

df_S_all_2 = df_all_2[df_all_2['config_label']=='S']
df_C_all_2 = df_all_2[df_all_2['config_label']=='C']

df_group_2 = df_all_2.groupby(by = 'file_label')
df_fit_2 = pd.DataFrame(columns = fit_l_p_header)


for name, group in df_group_2:
    print(name)
    if name not in file_to_ignore:

        df_temp = group
        df_fit_2 = pd.concat([df_fit_2,fit_func_delta_L(df_temp,skip_pt=8) ], ignore_index=True)
df_fit_2.to_csv('/Users/yogehs/Downloads/SbcC/plots_paper/MSD/RAW_data/S_C_ALL_MSD_2_segs_L_P_zetta_fitted.csv')
#viz_all_MSD(df_S_all_2,'S');viz_all_MSD(df_C_all_2,'C')
#%%renormalise 

def add_linear_line(x, y, **kwargs):
    ax = plt.gca()  # Get the current axis
    x_vals = np.array(ax.get_xlim())
    y_offset = np.array(ax.get_ylim())
    print(y_offset,x_vals)
    #x_vals = x_vals[10:]
    y_vals =0.0001+ x_vals ** 0.75 # Example: y = 2x (change this to your desired linear equation)
    ax.plot(x_vals, y_vals, '--', color='black')

df_renorm_2 = pd.DataFrame()
file_labels = df_fit_2['file_label'].unique()
fil_labels = ['fil1','fil2']

for fl in file_labels:
    if fl not in file_to_ignore:
        print(fl)
        for fil in fil_labels:
        
            row = df_fit_2.query(f"file_label == '{fl}' and fil_label == '{fil}'")
            mask = ( df_all_2['file_label'] == fl) &( df_all_2['fil_label']==fil)
            L_p= row['L_P_fitted'].to_numpy()[0]
            L = row['L_fitted'].to_numpy()[0]
            zeta =row['zeta_fitted'].to_numpy()[0]
            
            df_temp = df_all_2[mask]
            #fit_func_delta_L(df_temp,skip_pt=12)
            tau_c_temp= 1.0/ ((4.73/L)**4 * (kbT*L_p/zeta))

            tauc_inv =  (4.73/L)**4 * (kbT*L_p)/zeta
            r_c_temp = 1.0/ ((90*L_p**2)* (1.0/L)**4 )

            tau_FR = (zeta/(2*kbT*L_p)) * pow((L/np.pi),4)

            r_c_FR = pow(L,4)/(90*pow(L_p,2))
            print(f"\ntau in loop :{row.iloc[0]['tau_c']},{tau_c_temp},{tauc_inv}")
            print(f"rc in loop :{row.iloc[0]['r_c']},{r_c_temp}")

            df_temp['dt_norm_indi_fit'] = df_temp['deltaT_s']/tau_FR
            
            df_temp['dR_norm_indi_fit'] =df_temp['MSD'] / r_c_temp 
            
            
            df_temp['L_fitted'] =L
            df_temp['L_p_fitted'] =L_p
            df_temp['zeta_fitted'] =zeta
            df_temp['tau_c'] =tau_c_temp
            df_temp['tau_c_inv'] =tauc_inv

            #df_fit_2.iloc[0,'tau_c']=tau_c_temp
    
            df_renorm_2 = pd.concat([df_renorm_2,df_temp ], ignore_index=True)
#%%append local LP valeus to df_fit 
df_all = pd.read_csv(f'{final_save_plt_path}/all_Lp_vs_S_raw_data_10nm_windowsize.csv')

#TODO do it bitch 
columns = df_fit_2.columns.to_list() +['mean_local_lp','median_local_lp','min_local_lp','max_local_lp']
file_labels = df_fit_2['file_label'].unique()
fil_labels = ['fil1','fil2']
df_fit_2_local = pd.DataFrame(columns=columns)
df_fit_2_local = df_fit_2
for i in range(len(df_fit_2)):
    fl,fil = df_fit_2.loc[i,'file_label'],df_fit_2.loc[i,'fil_label']
    print(fl,fil)
    mask = (df_all['fil_label']==fil) & (df_all['file_label']==fl)
    df_local_temp = df_all[mask]
    df_fit_2_local.loc[i,'mean_local_lp'] = df_local_temp['L_p'].mean()
    df_fit_2_local.loc[i,'median_local_lp'] = df_local_temp['L_p'].median()
    df_fit_2_local.loc[i,'min_local_lp'] = df_local_temp['L_p'].min()
    df_fit_2_local.loc[i,'max_local_lp'] = df_local_temp['L_p'].max()
    df_fit_2_local.loc[i,'std_local_lp'] = df_local_temp['L_p'].std()

df_fit_2_local.to_clipboard()
mask = (df_fit_2_local['config_label']=='r')+(df_fit_2_local['config_label']=='H')
df_temp = df_fit_2_local[mask]
df_temp['config_label'] = df_temp['config_label']+'_'+df_temp['file_label'].apply(lambda x: x.split('_')[-2])
df_fit_2_local[mask] = df_temp
#df_fit_2_local['config_label'][[mask]=df_temp['file_label'].apply(lambda x: x.split('_')[-2])
df_fit_2_local['log_zeta'] = np.log(df_fit_2_local['zeta_fitted'] )
#averages for the table 
df_table = df_fit_2_local.groupby(by=['fil_label','config_label'],as_index=False)
df_table =df_table.agg({'log_zeta':['mean','std','sem','median'],
                        'tau_c':['mean','std','sem','median'],
                        'r_c':['mean','std','sem','median'],
                        'mean_local_lp':['mean','std','sem','median'],
                         'median_local_lp':['mean','std','sem','median'],  
                        'L_P_fitted':['mean','std','sem','median'],
                        'L_fitted':['mean','std','sem','median']})



df_fit_2_local.to_csv(f'{final_save_plt_path}/MSD/raw_local_MSD_fit_parameters')
df_table.to_csv(f'{final_save_plt_path}/MSD/stats_local_MSD_fit_parameters')
                        
#%%
R = sns.relplot(df_renorm_2,x='dt_norm_indi_fit',y='dR_norm_indi_fit',col = 'config_label',s=100,legend = True,hue = 'config_label')  

R.set(xlabel = None,xscale="log", yscale="log",box_aspect = 1)#,ylim=(0.01,10),xlim= (0.005,200))

plt.legend()
R = sns.relplot(df_renorm_2,x='deltaT_s',y='MSD',col = 'config_label',s=100,legend = True,hue = 'file_label')  

R.set(xlabel = None,xscale="log", yscale="log",box_aspect = 1)#,ylim=(0.01,10),xlim= (0.005,200))
plt.legend()
x_limits = plt.xlim()  # Get current x-axis limits
y_limits = plt.ylim()
print(x_limits,y_limits)

sns.displot(df_fit_2 ,x= 'L_P_fitted',row = 'config_label')
sns.displot(df_fit_2 ,x= 'zeta_fitted',row = 'config_label',log_scale=True)
sns.displot(df_fit_2 ,x= 'L_fitted',row = 'config_label')
df_group_2 = df_renorm_2.groupby(by = 'config_label')
#%%config by mean values 

import pandas as pd 
df_fit_2 = pd.read_csv('add you path to the results file here')
df_group_2 = df_fit_2.groupby(by = ['config_label','fil_label'])

for name, group in df_group_2:
    print(name,"Stats")
    print(f"mean lp_ {group['L_P_fitted'].mean()}__std_{group['L_P_fitted'].std()}")
    print(f"mean L_ {group['L_fitted'].mean()}__std_{group['L_fitted'].std()}")
    print(f"mean zeta_ {group['zeta_fitted'].mean()}__std_{group['zeta_fitted'].std()}")

    df_temp = group
    
    #plt.show()
#%%
df_group_2 = df_renorm_2.groupby(by = 'config_label')


for name, group in df_group_2:
    print(name)
    df_temp = group
    sns.scatterplot(df_temp, x = 'dt_norm_indi_fit',y = 'dR_norm_indi_fit',hue = 'config_label')
    plt.xscale('log');plt.yscale('log')#;plt.show()
    #plt.show()
    
#%%pred FR method
X_plot = df_renorm_2['dt_norm_indi_fit'].to_numpy()
Y_R = df_renorm_2['dR_norm_indi_fit'].to_numpy()

df_renorm_2 = df_renorm_2.sort_values( by="deltaT_s")

X_unnrom = df_renorm_2['deltaT_s'].to_numpy()
#X_unnrom = np.linspace(0.01,2000,10000)
L_avg_fit = geometric_mean(df_renorm_2['L_fitted'].to_numpy())
L_p_avg_fit = geometric_mean(df_renorm_2['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_renorm_2['zeta_fitted'].to_numpy())

#step 1 pred with delta s and avg values MSD 
MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
#step 2 rescale X_unnrom with avg vales for plotting

X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
#step 3 rescale MSD_pred_avg with avg vales for plotting
MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
#plotting
plt.scatter(X_plot,Y_R,c = 'purple')

plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)


plt.xscale('log');plt.yscale('log')
scale = X_unnrom_rescaled[1:100]

plt.plot(scale, 0.01+scale**0.75,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

plt.legend()


plt.ylabel("F, rescaled MSD In rR");plt.xlabel("t/T_c rescaled delta time")

x_limits = plt.xlim()  # Get current x-axis limits
y_limits = plt.ylim()
print(x_limits,y_limits)
#%%tauc L 

sns.relplot(df_renorm_2 ,x= 'L_fitted',y = 'tau_c',col = 'config_label',hue = 'file_label')
plt.xscale('log');plt.yscale('log');plt.show()


sns.distplot(df_renorm_2,y = 'tau_c')#;plt.yscale('log')
plt.show()
sns.histplot(df_renorm_2,y = 'zeta_fitted')

geometric_mean(df_renorm_2['zeta_fitted'].to_numpy())
#%%fig 2 dynamcis S,c confirmations
mask_other = (df_renorm_2['config_label']=='S') + (df_renorm_2['config_label']=='C') 
df_other = df_renorm_2[mask_other]

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),sharex = True,sharey=True)

axx_name_dict = {'S':axs[0],'C':axs[1]}
df_group_O = df_other.groupby(by=('config_label'))

for name, df_mask in df_group_O:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_mask['dt_norm_indi_fit'].to_numpy()
    MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
    ax_temp.scatter(t_data, MSD_R)
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('t/T_c rescaled delta time');ax_temp.set_ylabel('rescaled MSD In R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);

    
#for inset unscaled
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),sharex = True,sharey=True)
axx_name_dict = {'S':axs[0],'C':axs[1]}

for name, df_mask in df_group_O:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_mask['deltaT_s'].to_numpy()
    MSD_R = df_mask['MSD'].to_numpy()
    ax_temp.scatter(t_data, MSD_R,color = 'black')
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('Time');ax_temp.set_ylabel('MSD in R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);
#%%fig loop on color 
'''
TODO to do on monday
In that case, 
I would just increase the font size of the circles of rescaled MSD ,
just to distinguish from the unscaled ones
'''
plt.rcParams.update({'font.size': 50,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6,
                     'font.family':'Times New Roman','axes.linewidth':2.5,
                     'xtick.major.width': 2,'ytick.major.width': 2,'xtick.minor.width': 2,'ytick.minor.width': 2,'xtick.minor.visible':True})


x_renorm_limits = [0.030644186650635584, 4996.008916478943]
y_renorm_limit = [0.06738070130927563, 2.9694807353946584]

x_un_limits = [0.23042308789504504, 78.19536741138113]
y_un_limits = [0.7854736945816706, 148.95038821847717]

mask_other = (df_renorm_2['config_label']=='b') + (df_renorm_2['config_label']=='m') +(df_renorm_2['config_label']=='r') 
df_other = df_renorm_2[mask_other]


df_group_O = df_renorm_2.groupby(by=('config_label'))

for name, df_mask in df_group_O:
    print(name)
    if name =='H':
        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]

            plt.figure(figsize=(12, 10))
            plt.title(f"{name[0]}_MSD")
            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")
        
            X_unnrom = df_mask['deltaT_s'].to_numpy()
            L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
            L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
            
            #step 1 pred with delta s and avg values MSD 
            MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
            X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
            #step 3 rescale MSD_pred_avg with avg vales for plotting
            MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
            
            
            plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            plt.xscale('log');plt.yscale('log')
            plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
            plt.savefig(f'{final_save_plt_path}/{name}_MSD.svg')

            plt.show()
    else:
        color_i = color_saha_configs[name]

        plt.figure(figsize=(12, 10))
        plt.title(f"{name}_MSD")
        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")
    
        X_unnrom = df_mask['deltaT_s'].to_numpy()
        L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
        L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
        
        #step 1 pred with delta s and avg values MSD 
        MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
        X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
        #step 3 rescale MSD_pred_avg with avg vales for plotting
        MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
        
        
        plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
    
        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        plt.xscale('log');plt.yscale('log')
        plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
        plt.savefig(f'{final_save_plt_path}/{name}_MSD.svg')
        plt.show()
        

#for all rescaled in one plot
fig_h = plt.figure(2,figsize=(12, 10))
ax_h = fig_h.gca()

fig_sb = plt.figure(3,figsize=(12, 10))
ax_sb = fig_sb.gca()

#plt.title(f"all_MSD")

for name, df_mask in df_group_O:
    if name =='H':
        
        print(name,'bitch')

        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        #plt.figure(3,figsize=(12, 10))

        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]

            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")
        
            X_unnrom = df_mask['deltaT_s'].to_numpy()
            L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
            L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
            
            #step 1 pred with delta s and avg values MSD 
            MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
            X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
            #step 3 rescale MSD_pred_avg with avg vales for plotting
            MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
            
            
            ax_h.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            ax_h.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            ax_h.set_xscale('log'); ax_h.set_yscale('log')
            ax_h.set_xlim(x_renorm_limits);ax_h.set_ylim(y_renorm_limit)
        
        #plt.show()

    else:
        print(name,'bitch')
        color_i = color_saha_configs[name]

        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")
    
        X_unnrom = df_mask['deltaT_s'].to_numpy()
        L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
        L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
        
        #step 1 pred with delta s and avg values MSD 
        MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
        X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
        #step 3 rescale MSD_pred_avg with avg vales for plotting
        MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
        
        
        ax_sb.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        ax_sb.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
    
        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        ax_sb.set_xscale('log'); ax_sb.set_yscale('log')
        ax_sb.set_xlim(x_renorm_limits);ax_sb.set_ylim(y_renorm_limit)
fig_sb.set_dpi(300)
fig_h.set_dpi(300)
       
fig_sb.savefig(f'{final_save_plt_path}/all_abcc_MSD.svg')   
fig_h.savefig(f'{final_save_plt_path}/all_human_MSD.svg')   

#plt.savefig(f'{final_save_plt_path}/all_MSD.svg')
#%%fig loop color map

plt.rcParams.update({'font.size': 50,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':20,'ytick.major.size':20,'xtick.minor.size':10,'ytick.minor.size':10,
                     'font.family':'Times New Roman','axes.linewidth':2.5})

X_unnrom = np.logspace(-3.0, 4.5, num=1000)

#X_unnrom = np.linspace(0.05, 5000, 10000)


x_renorm_limits = [0.010644186650635584, 5000.008916478943]
y_renorm_limit = [0.06738070130927563, 2.9694807353946584]

x_un_limits = [0.23042308789504504, 78.19536741138113]
y_un_limits = [0.7854736945816706, 148.95038821847717]

mask_other = (df_renorm_2['config_label']=='b') + (df_renorm_2['config_label']=='m') +(df_renorm_2['config_label']=='r') 
df_other = df_renorm_2[mask_other]


df_group_O = df_renorm_2.groupby(by=('config_label'))

for name, df_mask in df_group_O:
    print(name)
    if name =='H':
        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]
            color_map  =cmap_saha_configs[name[0]]
            
            plt.figure(figsize=(12, 10));plt.title(f"{name[0]}_MSD")
            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")
        
            L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
            L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
            tau_c = df_mask['tau_c'].to_numpy()[0]
            #step 1 pred with delta s and avg values MSD 
            MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
            tau_FR = (zeta_avg_fit/(2*kbT*L_p_avg_fit)) * pow((L_avg_fit/np.pi),4)
            #tau_FR = (zeta_avg_fit/(2*kbT*L_p_avg_fit)) * pow((L_avg_fit/5),4)

            X_unnrom_rescaled = X_unnrom/tau_FR

            #X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
            #step 3 rescale MSD_pred_avg with avg vales for plotting
            MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
            #storing the arrys in the corresponding data frame
            
            #df_mask['MSD_pred_avg_rescaled'] = MSD_pred_avg_rescaled
            i_file=0;cmap = plt.get_cmap(color_map)
        
            for fn, df_file_file in df_mask.groupby(['file_label']):
                print(fn)
                t_data = df_file_file['dt_norm_indi_fit'].to_numpy()
                MSD_R = df_file_file['dR_norm_indi_fit'].to_numpy()
                plt.scatter(t_data, MSD_R,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
                i_file+=1

            #plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = color_i,label = 'predicted FR method',linewidth=5)#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            plt.xscale('log');plt.yscale('log')
            plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
            plt.savefig(f'{final_save_plt_path}/MSD/{name}_MSD.svg')

            plt.show()
    else:
        color_i = color_saha_configs[name]
        color_map  =cmap_saha_configs[name]

        plt.figure(figsize=(12, 10));plt.title(f"{name}_MSD")
        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")
    
        L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
        L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
        tau_c = df_mask['tau_c'].to_numpy()[0]

        #step 1 pred with delta s and avg values MSD 
        MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
        tau_FR = (zeta_avg_fit/(2*kbT*L_p_avg_fit)) * pow((L_avg_fit/np.pi),4)

        X_unnrom_rescaled = X_unnrom/tau_FR
        
        #X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
        #step 3 rescale MSD_pred_avg with avg vales for plotting
        MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
        
        i_file=0;cmap = plt.get_cmap(color_map)
    
        for fn, df_file_file in df_mask.groupby(['file_label']):
            print(fn)
            t_data = df_file_file['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_file_file['dR_norm_indi_fit'].to_numpy()
            
            plt.scatter(t_data, MSD_R,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
            i_file+=1

        #plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        
        plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = color_i,label = 'predicted FR method',linewidth=5)#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        plt.xscale('log');plt.yscale('log')
        plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
        plt.savefig(f'{final_save_plt_path}/MSD/{name}_MSD.svg')
        plt.show()
        

#for all rescaled in one plot
fig_h = plt.figure(2,figsize=(12, 10));ax_h = fig_h.gca()
ax_h.set_title("all human")
fig_sb = plt.figure(3,figsize=(12, 10));ax_sb = fig_sb.gca()
ax_sb.set_title("all sbcc")
from matplotlib.ticker import AutoMinorLocator
ax_sb.xaxis.set_minor_locator(AutoMinorLocator(4))
ax_sb.minorticks_on()

#plt.title(f"all_MSD")

for name, df_mask in df_group_O:
    if name =='H':
        
        print(name,'bitch')

        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        #plt.figure(3,figsize=(12, 10))

        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]
            color_map  =cmap_saha_configs[name[0]]

            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")
        
            L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
            L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
            
            #step 1 pred with delta s and avg values MSD 
            MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
            tau_FR = (zeta_avg_fit/(2*kbT*L_p_avg_fit)) * pow((L_avg_fit/np.pi),4)
            X_unnrom_rescaled = X_unnrom/tau_FR
            
            #X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
            #step 3 rescale MSD_pred_avg with avg vales for plotting
            MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
                        
            i_file=0;cmap = plt.get_cmap(color_map)
        
            for fn, df_file_file in df_mask.groupby(['file_label']):
                print(fn)
                t_data_file = df_file_file['dt_norm_indi_fit'].to_numpy()
                MSD_R_file = df_file_file['dR_norm_indi_fit'].to_numpy()
                ax_h.scatter(t_data_file, MSD_R_file,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
                i_file+=1

            
            #ax_h.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            ax_h.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method',linewidth=5)#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            ax_h.set_xscale('log'); ax_h.set_yscale('log')
            ax_h.set_xlim(x_renorm_limits);ax_h.set_ylim(y_renorm_limit)
        
        #plt.show()

    else:
        print(name,'bitch')
        color_i = color_saha_configs[name]
        color_map  =cmap_saha_configs[name]

        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")
    
        L_avg_fit = geometric_mean(df_mask['L_fitted'].to_numpy())
        L_p_avg_fit = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_mask['zeta_fitted'].to_numpy())
        
        #step 1 pred with delta s and avg values MSD 
        MSD_pred_avg =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
        tau_FR = (zeta_avg_fit/(2*kbT*L_p_avg_fit)) * pow((L_avg_fit/np.pi),4)
        X_unnrom_rescaled = X_unnrom/tau_FR
        
        #X_unnrom_rescaled = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
        #step 3 rescale MSD_pred_avg with avg vales for plotting
        MSD_pred_avg_rescaled =(90*L_p_avg_fit**2)* MSD_pred_avg * (1/L_avg_fit)**4 
        i_file=0;cmap = plt.get_cmap(color_map)
    
        for fn, df_file_file in df_mask.groupby(['file_label']):
            print(fn)
            t_data_file = df_file_file['dt_norm_indi_fit'].to_numpy()
            MSD_R_file = df_file_file['dR_norm_indi_fit'].to_numpy()
            ax_sb.scatter(t_data_file, MSD_R_file,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
            i_file+=1

        #ax_sb.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        ax_sb.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method',linewidth=5)#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
    
        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        ax_sb.set_xscale('log'); ax_sb.set_yscale('log')
        ax_sb.set_xlim(x_renorm_limits);ax_sb.set_ylim(y_renorm_limit)
fig_sb.set_dpi(300)
fig_h.set_dpi(300)
       
fig_sb.savefig(f'{final_save_plt_path}/MSD/all_abcc_MSD.svg')   
fig_h.savefig(f'{final_save_plt_path}/MSD/all_human_MSD.svg')   

#plt.savefig(f'{final_save_plt_path}/all_MSD.svg')

#%%for unscaled in one plots 
from matplotlib import cm
from matplotlib.colors import Normalize
plt.rcParams.update({'font.size': 80})
# Normalize the data so it maps to the colormap properly
#norm = Normalize(vmin=min(np.min(y1), np.min(y2), np.min(y3)), vmax=max(np.max(y1), np.max(y2), np.max(y3)))

# Choose a colormap (e.g., 'viridis')
#colormap = cm.get_cmap('viridis')
for name, df_mask in df_group_O:

    print(name)
    if name =='H':
        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        
        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open

            color_i = color_saha_configs[name[0]]
            color_map  =cmap_saha_configs[name[0]]
            plt.figure(figsize=(12, 10))


            plt.title(f"{name[0]}_MSD")
            i_file=0
            cmap = plt.get_cmap(color_map)
        
            for name, df_mask1 in df_mask.groupby(['file_label']):
                print(name)
                t_data = df_mask1['deltaT_s'].to_numpy()
                MSD_R = df_mask1['MSD'].to_numpy()
                
        
                plt.scatter(t_data, MSD_R,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
                
                #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
                plt.xscale('log');plt.yscale('log')
                ax = plt.gca()
                #ax.spines[['right', 'top']].set_visible(False)
        
                #plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
                plt.xlim(x_un_limits),plt.ylim(y_un_limits)
                i_file+=1
            plt.savefig(f'{final_save_plt_path}/MSD/{name}_unscaled_MSD.svg')

            plt.show()

    else:
        print(name)
        color_i = color_saha_configs[name]
        color_map  =cmap_saha_configs[name]
        plt.figure(figsize=(12, 10))
        #plt.title(f'{name}')

        plt.title(f"{name}_MSD")
        i_file=0
        cmap = plt.get_cmap(color_map)
    
        for name, df_mask1 in df_mask.groupby(['file_label']):
            print(name)
            t_data = df_mask1['deltaT_s'].to_numpy()
            MSD_R = df_mask1['MSD'].to_numpy()
            
    
            plt.scatter(t_data, MSD_R,c=cmap(0.2 + i_file*0.2),edgecolors=color_i)
            
            #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
            plt.xscale('log');plt.yscale('log')
            ax = plt.gca()
            #ax.spines[['right', 'top']].set_visible(False)
    
            #plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
            plt.xlim(x_un_limits),plt.ylim(y_un_limits)
            i_file+=1
        plt.savefig(f'{final_save_plt_path}/MSD/{name}_unscaled_MSD.svg')

        plt.show()

#%%fig 3 dynamcis other confirmations
mask_other = (df_renorm_2['config_label']=='b') + (df_renorm_2['config_label']=='m') +(df_renorm_2['config_label']=='r') 
df_other = df_renorm_2[mask_other]

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 6),sharex = True,sharey=True)

axx_name_dict = {'b':axs[0],'r':axs[1],'m':axs[2]}
df_group_O = df_other.groupby(by=('config_label'))

for name, df_mask in df_group_O:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_mask['dt_norm_indi_fit'].to_numpy()
    MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
    ax_temp.scatter(t_data, MSD_R)
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('t/T_c rescaled delta time');ax_temp.set_ylabel('rescaled MSD In R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);

    
#for inset unscaled
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 18),sharex = True,sharey=True)
axx_name_dict = {'b':axs[0],'r':axs[1],'m':axs[2]}

for name, df_mask in df_group_O:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_mask['deltaT_s'].to_numpy()
    MSD_R = df_mask['MSD'].to_numpy()
    ax_temp.scatter(t_data, MSD_R,color = 'black')
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('Time');ax_temp.set_ylabel('MSD in R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);
    
#%%fig four human 
df_human = df_renorm_2[df_renorm_2['config_label']=='H']
df_human['open_label']=df_human['file_label'].apply(lambda x: x.split('_')[2])

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12),sharex = True,sharey=True)

axx_name_dict = {'closed':axs[0],'open':axs[1]}
df_group_H = df_human.groupby(by=('open_label'))

for name, df_open in df_group_H:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_open['dt_norm_indi_fit'].to_numpy()
    MSD_R = df_open['dR_norm_indi_fit'].to_numpy()
    ax_temp.scatter(t_data, MSD_R)
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('t/T_c rescaled delta time');ax_temp.set_ylabel('rescaled MSD In R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);

    
#for inset unscaled
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12),sharex = True,sharey=True)

axx_name_dict = {'closed':axs[0],'open':axs[1]}
df_group_H = df_human.groupby(by=('open_label'))

for name, df_open in df_group_H:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_open['deltaT_s'].to_numpy()
    MSD_R = df_open['MSD'].to_numpy()
    ax_temp.scatter(t_data, MSD_R,color = 'black')
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('Time');ax_temp.set_ylabel('MSD in R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);
#%%fig loop predict from all configs 

'''
TODO to do on monday
In that case, 
I would just increase the font size of the circles of rescaled MSD ,
just to distinguish from the unscaled ones
'''
plt.rcParams.update({'font.size': 50,'lines.markersize':14,'figure.dpi':200,'lines.color':'black',
                     'xtick.major.size':10,'ytick.major.size':10,'xtick.minor.size':6,'ytick.minor.size':6,
                     'font.family':'Times New Roman','axes.linewidth':2.5,
                     'xtick.major.width': 2,'ytick.major.width': 2,'xtick.minor.width': 2,'ytick.minor.width': 2,'xtick.minor.visible':True})

#X_unnrom = df_renorm_2['deltaT_s'].to_numpy()
X_unnrom = np.linspace(0.1,2000,10000)
L_avg_fit = geometric_mean(df_renorm_2['L_fitted'].to_numpy())
L_p_avg_fit = geometric_mean(df_renorm_2['L_p_fitted'].to_numpy());zeta_avg_fit = geometric_mean(df_renorm_2['zeta_fitted'].to_numpy())

#step 1 pred with delta s and avg values MSD 
MSD_pred_avg_all_config =  delta_R(X_unnrom, zeta_avg_fit,L_p_avg_fit,L_avg_fit )
X_unnrom_rescaled_all_config = X_unnrom * (4.74/L_avg_fit)**4 * (kbT*L_p_avg_fit/zeta_avg_fit) 
#step 3 rescale MSD_pred_avg with avg vales for plotting
MSD_pred_avg_rescaled_all_config =(90*L_p_avg_fit**2)* MSD_pred_avg_all_config * (1/L_avg_fit)**4 



mask_other = (df_renorm_2['config_label']=='b') + (df_renorm_2['config_label']=='m') +(df_renorm_2['config_label']=='r') 
df_other = df_renorm_2[mask_other]


df_group_O = df_renorm_2.groupby(by=('config_label'))

for name, df_mask in df_group_O:
    print(name)
    if name =='H':
        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]
            X_unnrom = np.linspace(0.1,2000,10000)

            plt.figure(figsize=(12, 10))
            plt.title(f"{name[0]}_MSD")
            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")
            L_avg_conf = geometric_mean(df_mask['L_fitted'].to_numpy())
            L_p_avg_conf = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_conf = geometric_mean(df_mask['zeta_fitted'].to_numpy())
            X_unnrom_rescaled = X_unnrom * (4.74/L_avg_conf)**4 * (kbT*L_p_avg_conf/zeta_avg_conf) 

            #step 1 pred with delta s and avg values MSD 
            MSD_pred_conf =  delta_R(X_unnrom, zeta_avg_conf,L_p_avg_conf,L_avg_conf )
            MSD_pred_avg_rescaled =(90*L_p_avg_conf**2)* MSD_pred_conf * (1/L_avg_conf)**4 
            plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = color_i,label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

           
            plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            plt.plot(X_unnrom_rescaled_all_config, MSD_pred_avg_rescaled_all_config,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            plt.xscale('log');plt.yscale('log')
            plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
            plt.savefig(f'{final_save_plt_path}/{name}_MSD.svg')

            plt.show()
    else:
        color_i = color_saha_configs[name]
        X_unnrom = np.linspace(0.1,2000,10000)

        plt.figure(figsize=(12, 10))
        plt.title(f"{name}_MSD")
        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")
        L_avg_conf = geometric_mean(df_mask['L_fitted'].to_numpy())
        L_p_avg_conf = geometric_mean(df_mask['L_p_fitted'].to_numpy());zeta_avg_conf = geometric_mean(df_mask['zeta_fitted'].to_numpy())
        X_unnrom_rescaled = X_unnrom * (4.74/L_avg_conf)**4 * (kbT*L_p_avg_conf/zeta_avg_conf) 

        #step 1 pred with delta s and avg values MSD 
        MSD_pred_conf =  delta_R(X_unnrom, zeta_avg_conf,L_p_avg_conf,L_avg_conf )
        MSD_pred_avg_rescaled =(90*L_p_avg_conf**2)* MSD_pred_conf * (1/L_avg_conf)**4 
        plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = color_i,label = 'predicted FR method',linewidth = 6)#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

        
        plt.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        plt.plot(X_unnrom_rescaled_all_config, MSD_pred_avg_rescaled_all_config,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        #plt.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = color_i,label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)

        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        plt.xscale('log');plt.yscale('log')
        plt.xlim(x_renorm_limits),plt.ylim(y_renorm_limit)
        plt.savefig(f'{final_save_plt_path}/{name}_MSD.svg')
        plt.show()
        

#for all rescaled in two plots human/sbcc

fig_h = plt.figure(2,figsize=(12, 10))
ax_h = fig_h.gca()

fig_sb = plt.figure(3,figsize=(12, 10))
ax_sb = fig_sb.gca()

#plt.title(f"all_MSD")

for name, df_mask in df_group_O:
    if name =='H':
        
        print(name,'bitch')

        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        #plt.figure(3,figsize=(12, 10))

        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]

            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")

            ax_h.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            ax_h.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            ax_h.set_xscale('log'); ax_h.set_yscale('log')
            ax_h.set_xlim(x_renorm_limits);ax_h.set_ylim(y_renorm_limit)
        
        #plt.show()

    else:
        print(name,'bitch')
        color_i = color_saha_configs[name]

        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")

        ax_sb.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        ax_sb.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
    
        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        ax_sb.set_xscale('log'); ax_sb.set_yscale('log')
        ax_sb.set_xlim(x_renorm_limits);ax_sb.set_ylim(y_renorm_limit)
fig_sb.set_dpi(300)
fig_h.set_dpi(300)
       
fig_sb.savefig(f'{final_save_plt_path}/all_abcc_MSD.svg')   
fig_h.savefig(f'{final_save_plt_path}/all_human_MSD.svg')   



#for all rescaled in two plots human/sbcc

fig_h_sb = plt.figure(4,figsize=(12, 10))
ax_h_sb = fig_h_sb.gca()



#plt.title(f"all_MSD")

for name, df_mask in df_group_O:
    if name =='H':
        
        print(name,'bitch')

        df_mask['open_label'] = df_mask['file_label'].apply(lambda x: x.split('_')[2])
        df_group_open = df_mask.groupby(['open_label'])
        #plt.figure(3,figsize=(12, 10))

        for name,df_mask_conf_open in df_group_open:
            df_mask = df_mask_conf_open
            color_i = color_saha_configs[name[0]]

            t_data = df_mask['dt_norm_indi_fit'].to_numpy()
            MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
            df_mask = df_mask.sort_values( by="deltaT_s")

            ax_h_sb.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
            ax_h_sb.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
        
            #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
            ax_h_sb.set_xscale('log'); ax_h_sb.set_yscale('log')
            ax_h_sb.set_xlim(x_renorm_limits);ax_h_sb.set_ylim(y_renorm_limit)
        
        #plt.show()

    else:
        print(name,'bitch')
        color_i = color_saha_configs[name]

        t_data = df_mask['dt_norm_indi_fit'].to_numpy()
        MSD_R = df_mask['dR_norm_indi_fit'].to_numpy()
        df_mask = df_mask.sort_values( by="deltaT_s")

        ax_h_sb.scatter(t_data, MSD_R,edgecolor=color_i,facecolor=color_i,alpha=0.45,linewidth=2 )
        ax_h_sb.plot(X_unnrom_rescaled, MSD_pred_avg_rescaled,color = 'black',label = 'predicted FR method')#, label=f'_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}',color = 'black',zorder = 3)
    
        #plt.xlabel('t/T_c rescaled delta time');plt.ylabel('rescaled MSD In R')
        ax_h_sb.set_xscale('log'); ax_h_sb.set_yscale('log')
        ax_h_sb.set_xlim(x_renorm_limits);ax_h_sb.set_ylim(y_renorm_limit)
fig_h_sb.set_dpi(300)
       
fig_h_sb.savefig(f'{final_save_plt_path}/all_human_sbcc_MSD.svg')   
#%% ring with ad W/O ATP
df_ring = df_renorm_2[df_renorm_2['config_label']=='r']

df_ring_fit_2 = df_fit_2[df_fit_2['config_label']=='r']
df_ring['atp_label']=df_ring['file_label'].apply(lambda x: x.split('_')[-2])
df_ring_fit_2['atp_label']=df_ring_fit_2['file_label'].apply(lambda x: x.split('_')[-2])

df_group_2 = df_ring_fit_2.groupby(by = ['atp_label','fil_label'])

for name, group in df_group_2:
    print(name,"Stats")
    print(f"mean lp_ {group['L_P_fitted'].mean()}__std_{group['L_P_fitted'].std()}")
    print(f"mean L_ {group['L_fitted'].mean()}__std_{group['L_fitted'].std()}")
    print(f"mean zeta_ {group['zeta_fitted'].mean()}__std_{group['zeta_fitted'].std()}")

    df_temp = group


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12),sharex = True,sharey=True)

axx_name_dict = {'ATP':axs[0],'noATP':axs[1]}
df_group_H = df_ring.groupby(by=('atp_label'))

for name, df_open in df_group_H:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_open['dt_norm_indi_fit'].to_numpy()
    MSD_R = df_open['dR_norm_indi_fit'].to_numpy()
    ax_temp.scatter(t_data, MSD_R)
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('t/T_c rescaled delta time');ax_temp.set_ylabel('rescaled MSD In R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    ax_temp.set_xlim(x_renorm_limits);axs[0].set_ylim(y_renorm_limit);

    
#for inset unscaled
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12),sharex = True,sharey=True)

axx_name_dict = {'ATP':axs[0],'noATP':axs[1]}
df_group_H = df_ring.groupby(by=('atp_label'))

for name, df_open in df_group_H:
    print(name)
    ax_temp = axx_name_dict[name]
    ax_temp.set_title(f"{name}_MSD")
    t_data = df_open['deltaT_s'].to_numpy()
    MSD_R = df_open['MSD'].to_numpy()
    ax_temp.scatter(t_data, MSD_R,color = 'black')
    #ax_temp.plot(t_data, delta_R(t_data,  zeta_fitted,lp_fitted,L_fitted), label=f'{i+1}_zeta_{zeta_fitted:.2f}_L_p_{lp_fitted:.2f}', color=color_arr[i])
    ax_temp.set_xlabel('Time');ax_temp.set_ylabel('MSD in R')
    ax_temp.set_xscale('log');ax_temp.set_yscale('log')
    #ax_temp.set_ylim(y_min,y_max);axs[0].set_xlim(0,80);
#%%zeta values 
plt.figure(1,figsize=(18,10))
sns.scatterplot(data=df_fit_2 , x = 'config_label',y='zeta_fitted',hue='config_label')
plt.yscale('log')

plt.figure(1,figsize=(18,10))
sns.boxplot(data=df_fit_2 , x = 'config_label',y='zeta_fitted',hue='config_label')
plt.yscale('log')

plt.figure(2,figsize=(18,10))
sns.scatterplot(data=df_fit_2 , x = 'config_label',y='L_P_fitted',hue='config_label')

plt.figure(3,figsize=(18,10))
sns.boxplot(data=df_fit_2 , x = 'config_label',y='L_P_fitted',hue='config_label')
#%%time trace of the MSD 

file_labels = df_fit_2['file_label'].unique()
fil_labels = ['fil1','fil2']

for fl in file_labels:
    plt.figure(fl);plt.clf()
    for fil in fil_labels:
        
        row = df_fit_2.query(f"file_label == '{fl}' and fil_label == '{fil}'")
        mask = ( df_all_2['file_label'] == fl) &( df_all_2['fil_label']==fil)
        
        df_temp = df_all_2[mask]
        
        df_temp['dt_norm_indi_fit'] = df_temp['deltaT_s']/tau_c_temp
        
        df_temp['dR_norm_indi_fit'] =df_temp['MSD'] / r_c_temp 
        
        
        plt.scatter(df_temp['deltaT_s'],df_temp['MSD'])