#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:35:27 2024

@author: yogehs
"""

import posixpath  # to generate unix paths
from pathlib2 import PurePath, PureWindowsPath, PurePosixPath

import re
from pyqtgraph.Qt import  QtCore
import pyqtgraph as pg
import numpy as np 
import os
from scipy.ndimage.filters import  uniform_filter1d
import matplotlib.pyplot as plt 
from scipy.stats import norm
from collections import defaultdict
import pandas as pd
from constants import *
# Get key mappings from Qt namespace
qt_keys = (
    (getattr(QtCore.Qt, attr), attr[4:])
    for attr in dir(QtCore.Qt)
    if attr.startswith("Key_")
)
keys_mapping = defaultdict(lambda: "unknown", qt_keys)
class KeyPressWindow(pg.GraphicsLayoutWidget):
    sigKeyPress = QtCore.pyqtSignal(object)

    def keyPressEvent(self, ev):
        self.scene().keyPressEvent(ev)
        self.sigKeyPress.emit(ev)
        
def Rollavg(Z, N, mode):      # rolling average on N time steps
 #   Zavg=np.convolve(Z, np.ones((N,))/N, mode=mode)
    Zavg = uniform_filter1d(Z, size=N, mode=mode_roll)
 #   n=int(N/2.); Zavg[range(n)]=Zavg[n+1]; Zavg[range(-1,-n-1,-1)]= Zavg[-n]
 #   if mode=='valid': Z=Z[0: len(Z)-N+1];# Z=Z[0: len(Z)-N+1]
    return Z, Zavg

def errplot(x, y, yerr, **kwargs):
    ax = plt.gca()
    data = kwargs.pop("data")
    data.plot(x=x, y=y, yerr=yerr, kind="scatter", ax=ax, **kwargs)
def violin_plt_arr(df_msd_arr,df_msd_mask):
    arc_l = df_msd_mask['avg_arc_L'].unique()
    fig = plt.figure()
    fil_label = ['fil1','fil2']
    gs = fig.add_gridspec(2, len(arc_l), hspace=0, wspace=0)
    ax_arr = gs.subplots(sharey= True,sharex =True)
    fig.set_size_inches(22, 8)
    
    for i in range(2):
        for j in range(len(arc_l)):
            data_lab_i = f"{fil_label[i]}_{arc_l[j]}"
            mask = (df_msd_arr['data_label'] ==data_lab_i)
            df_mask = df_msd_arr[mask]
            print(df_mask['MSD_arr'])
            sns.violinplot(data=df_mask,x='deltaT',y="L_p_deviation_arr",inner='points',ax=ax_arr[i][j])
def extract_number(folder_name):
    match = re.search(r'\d+', folder_name)
    return int(match.group()) if match else 0           
def find_MSD_metric(df_rdf):
    '''
    input hte dfrdf from the methods class 
    does the computations for only half the Delta T
    return DF_MSD with radial MSD and also ang MSD
    
    '''
    df_msd = pd.DataFrame()    
    df_msd_arr = pd.DataFrame()        

    #df_msd["MSD_arr"] = df_msd['MSD_arr'].astype(object)
    #df_msd["ang_MSD_arr"] = df_msd['ang_MSD_arr'].astype(object)


    uni_arc_L = df_rdf['avg_arc_L'].unique()
    fil_labels = ['fil1','fil2']
    file_label_temp =df_rdf['file_label'][0]
    fps = long_vids_fps[file_label_temp]
    for la in uni_arc_L:
        for fl in fil_labels:
            mask_iso = (df_rdf['avg_arc_L']== la) & (df_rdf['fil_label']== fl)
            df_rdf_iso_arc = df_rdf[mask_iso ]
            end_end_dist_iso = df_rdf_iso_arc['end_end_dist'].to_numpy()
            tan_ang_iso =         df_rdf_iso_arc['tangential_angle_rad'].to_numpy()
            L_p_arr =         df_rdf_iso_arc['L_p'].to_numpy()
            bend_ang_iso =         df_rdf_iso_arc['bend_angle_deg'].to_numpy()
            
            L_p_R_sq_arr =         df_rdf_iso_arc['L_p_Rsq'].to_numpy()
            angle_wrt_r_iso = df_rdf_iso_arc['angle_wrt_r'].to_numpy()
            
            
            N_2 = int(len(L_p_arr)*0.5)

            L_p_arr_N2 = L_p_arr[:N_2];L_p_temp = np.mean(L_p_arr_N2)
            L_arr =  df_rdf_iso_arc['l_seg'].to_numpy()
            mask_good_L_P = L_p_R_sq_arr[:N_2] >=0.99
            
            #print(la,fl,len(L_p_arr),sum(mask_good_L_P),N_2)
            #loop to find the 
            for i in range(1,len(end_end_dist_iso)):

                deviation =          end_end_dist_iso[i:] -end_end_dist_iso[:-i]
                tan_ang_deviation =  tan_ang_iso[i:] - tan_ang_iso[:-i]
                bend_ang_deviation =       bend_ang_iso[i:] - bend_ang_iso[:-i]
                angle_wrt_r_deviation =       angle_wrt_r_iso[i:] - angle_wrt_r_iso[:-i]

                L_p_deviation = L_p_arr[i:] - L_p_arr[:-i]

                deviation_2 = (deviation**2).astype(float)
                deviation_2 = np.round(deviation_2,2)
                MSD = np.mean(deviation_2);std_MSD = (np.std(deviation_2)/np.sqrt(len(deviation_2)-1))

                bend_ang_deviation_2= (bend_ang_deviation**2).astype(float)
                bend_ang_deviation_2 = np.round(bend_ang_deviation_2,4)
                bend_ang_MSD = np.mean(bend_ang_deviation_2);std_bend_ang_MSD = (np.std(bend_ang_deviation_2)/np.sqrt(len(bend_ang_deviation_2)-1))

                
                tan_ang_deviation_2= (tan_ang_deviation**2).astype(float)
                tan_ang_deviation_2 = np.round(tan_ang_deviation_2,4)
                tan_ang_MSD = np.mean(tan_ang_deviation_2);std_tan_ang_MSD = (np.std(tan_ang_deviation_2)/np.sqrt(len(tan_ang_deviation_2)-1))

                angle_wrt_r_deviation_2= (angle_wrt_r_deviation**2).astype(float)
                angle_wrt_r_deviation_2 = np.round(angle_wrt_r_deviation_2,4)
                ang_wrt_r_MSD = np.mean(angle_wrt_r_deviation_2);std_ang_wrt_r_MSD = (np.std(angle_wrt_r_deviation_2)/np.sqrt(len(angle_wrt_r_deviation_2)-1))


                L_p_deviation_2 = (L_p_deviation**2).astype(float)
                L_p_deviation_2 = np.round(L_p_deviation_2,4)
                
                L_p_MSD = np.mean(L_p_deviation_2);std_L_p_MSD = (np.std(L_p_deviation_2)/np.sqrt(len(L_p_deviation_2)-1))

    
                #saving the data frame
                #saving the raw arrays      
                if i+1<0.5*len(end_end_dist_iso):
                    df_temp = pd.DataFrame({'d_r_arr':deviation,
                                            'd_tan_ang_arr':tan_ang_deviation,
                                            'd_bend_ang_arr':bend_ang_deviation,
                                            'd_ang_wrt_r_arr':angle_wrt_r_deviation,

                                            'd_L_p_arr':L_p_deviation,
                                            'd_2_r_arr':deviation_2,
                                            'd_2_tan_ang_arr':tan_ang_deviation_2,
                                            'd_2_bend_ang_arr':bend_ang_deviation_2,
                                            'd_2_ang_wrt_r_arr':angle_wrt_r_deviation_2,

                                            'd_2_L_p_arr':L_p_deviation_2,
                                            'acg_arcl':f'{la}',
                                            'fil_label':f'{fl}',

                                            "data_label":f"{fl}_{la}",
                                            'deltaT':i,
                                            'log_deltaT':np.log10(i),
                                            'file_label':file_label_temp}) 
                    
                    
                    df_msd_arr = pd.concat([df_msd_arr, df_temp], axis=0) 

                L_temp = np.mean(L_arr[:N_2])
                N_df = len(df_msd)
                T_s = i/fps
                df_msd.loc[N_df,'file_label']  = file_label_temp
                df_msd.loc[N_df,'config_label']  = file_label_temp[0]
                df_msd.loc[N_df,'folder_number']  = file_label_temp[1:]


                df_msd.loc[N_df,'data_label'] = f"{fl}_{la}"

                df_msd.loc[N_df,'deltaT']  = i;df_msd.loc[N_df,'log_deltaT']  = np.log(i)
                
                df_msd.loc[N_df,'deltaT_s']  =T_s ;df_msd.loc[N_df,'log_deltaT_s']  = np.log(T_s)

                #radial MSDs
                df_msd.loc[N_df,'MSD']  = MSD;df_msd.loc[N_df,'std_norm_MSD']  = std_MSD
                df_msd.loc[N_df,'std_MSD']  = np.std(deviation_2)
                df_msd.loc[N_df,'MSD_norm']  = (MSD*(90*L_p_temp**2))/(L_temp**4)
                df_msd.loc[N_df,'dT_MSD_norm']  = ((i*L_p_temp)/L_temp**4)
                df_msd.loc[N_df,'dT_MSD_norm_s']  = ((T_s*L_p_temp)/L_temp**4)

                #bending angle MSD 
                df_msd.loc[N_df,'bend_ang_MSD']  = bend_ang_MSD;df_msd.loc[N_df,'std_norm_bend_ang_MSD']  = std_bend_ang_MSD
                df_msd.loc[N_df,'std_bend_ang_MSD']  = np.std(bend_ang_deviation_2)

                #tangential angle MSD 
                df_msd.loc[N_df,'tan_ang_MSD']  = tan_ang_MSD;df_msd.loc[N_df,'std_norm_tan_ang_MSD']  = std_tan_ang_MSD
                df_msd.loc[N_df,'tan_ang_MSD_norm']  = tan_ang_MSD/L_p_temp
                df_msd.loc[N_df,'std_tan_ang_MSD']  = np.std(tan_ang_deviation_2)
                
                 #tangential angle MSD 
                df_msd.loc[N_df,'ang_wrt_r_MSD']  = ang_wrt_r_MSD;df_msd.loc[N_df,'std_norm_ang_wrt_r_MSD']  =     std_ang_wrt_r_MSD
                df_msd.loc[N_df,'std_ang_wrt_r_MSD']  = np.std(angle_wrt_r_deviation_2)
                
                
                #LP MSD 
                df_msd.loc[N_df,'L_p_MSD']  = L_p_MSD;df_msd.loc[N_df,'std_norm_L_p_MSD']  = std_L_p_MSD
                df_msd.loc[N_df,'std_L_p_MSD']  = np.std(L_p_deviation_2)
                

                df_msd.loc[N_df,'mean_l_p']  = L_p_temp
                df_msd.loc[N_df,'l_p_i']  = L_p_arr[i]
                df_msd.loc[N_df,'l_p_rev_i']  = L_p_arr[-i]

                df_msd.loc[N_df,'mean_l_p_good']  = np.mean(L_p_arr_N2[mask_good_L_P])

                df_msd.loc[N_df,'std_l_p']  = np.std(L_p_arr[:N_2])
                df_msd.loc[N_df,'L_seg']  = L_temp
                
                
                

                df_msd.loc[N_df,'avg_arc_L']  = la;df_msd.loc[N_df,'fil_label']  = fl
    return(df_msd,df_msd_arr)

def path2unix(path, nojoin=True, fromwinpath=False):
    """From a path given in any format, converts to posix path format
    fromwinpath=True forces the input path to be recognized as a Windows path (useful on Unix machines to unit test Windows paths)"""
    if not path:
        return path
    if fromwinpath:
        pathparts = list(PureWindowsPath(path).parts)
    else:
        pathparts = list(PurePath(path).parts)
    if nojoin:
        return pathparts
    else:
        return posixpath.join(*pathparts)
def folder_labels_YS_def(file_path ,num_file_labels):
    file_labels_arr= path2unix(file_path)

    return (file_labels_arr[-num_file_labels:])

def df_fil1_fil_2_all_L_p_s(df_l_p):
    '''
    just hte returns the whole L_p  
    labels the fil1 as the one high LP by defualt for L_p vsLP

    Parameters
    ----------
    df_l_p : TYPE
        DESCRIPTION.

    Returns
    -------
    df_fil1_fil2 : TYPE
        DESCRIPTION.

    '''
    df_l_p = df_l_p[["file_label", "file_path",'frame_numer', 'max_arc_l','avg_arc_l_win_fil_1','L_p','config_label']]
    df_l_p['fil_label_L_p'] = "lolol"
    df_l_p['fil_label_arc'] = "lolol"

    
    
    mask_fil1 = df_l_p['avg_arc_l_win_fil_1']<0.5
    df_fil1 = df_l_p[mask_fil1] ; df_fil1['fil_label_arc'] = 'fil1'
    
    df_fil2 = df_l_p[~mask_fil1]; df_fil2['fil_label_arc'] = 'fil2'
    df_fil1['avg_arc_l_']  = df_fil1['avg_arc_l_win_fil_1'] 

    df_fil2['avg_arc_l_']  = 1 - df_fil2['avg_arc_l_win_fil_1'] 
    df_fil2 = df_fil2.sort_values(by=['avg_arc_l_']).reset_index(drop= True) 
    
    lp1_avg = df_fil1.loc[:, 'L_p'].mean()
    lp2_avg = df_fil2.loc[:, 'L_p'].mean()
    
    if lp1_avg>lp2_avg:
        df_fil1['fil_label_L_p'] = 'fil1';df_fil2['fil_label_L_p'] = 'fil2'
        swap_bool = False
    else:
        df_temp_swap = df_fil2
        df_fil2 = df_fil1
        df_fil1 = df_temp_swap
        df_fil1['fil_label_L_p'] = 'fil1';df_fil2['fil_label_L_p'] = 'fil2'
        swap_bool = True

        

    df_fil1['avg_arc_l_int'] = np.round(100*df_fil1['avg_arc_l_'],2).astype(int)
    
    df_fil2['avg_arc_l_int'] = np.round(100*df_fil2['avg_arc_l_'],2).astype(int)
    df_fil1_fil2 = pd.merge(df_fil1, df_fil2,suffixes=("_fil1","_fil2"),on=['avg_arc_l_int','frame_numer'],how = 'inner')
    return df_fil1_fil2,swap_bool
def df_fil1_fil_2_mean_lp_s(df_l_p):
    '''
    just hte mean_l_ps  

    Parameters
    ----------
    df_l_p : TYPE
        DESCRIPTION.

    Returns
    -------
    df_fil1_fil2 : TYPE
        DESCRIPTION.

    '''
    num_label = 2
    file_path = df_l_p['file_path'].to_numpy()[0]

    grouped = df_l_p.groupby("avg_arc_l_win_fil_1", as_index=False)['L_p']
    df_group=grouped.agg([ "mean", "std"]).rename(columns={"avg_arc_l_win_fil_1":"avg_arc_l_" ,"mean": "L_p_mean", "std": "L_p_std"})
    for i in range(num_label):
        
       df_group[f'file_label_{i}'] = folder_labels_YS_def(file_path,num_label)[i]
    df_group['file_label_2'] = df_group['file_label_0']+df_group['file_label_1']
    
    
    df_group['fil_label_arc'] = "hehe"
    df_group['fil_label_L_p'] = "lolol"
    
    
    mask_fil1 = df_group['avg_arc_l_']<0.5
    df_fil1 = df_group[mask_fil1] ; df_fil1['fil_label_arc'] = 'fil1'
    
    df_fil2 = df_group[~mask_fil1]; df_fil2['fil_label_arc'] = 'fil2'
    
    df_fil2['avg_arc_l_']  = 1 - df_fil2['avg_arc_l_'] 
    df_fil2 = df_fil2.sort_values(by=['avg_arc_l_']).reset_index(drop= True) 
    
    lp1_avg = df_fil1.loc[:, 'L_p_mean'].mean()
    lp2_avg = df_fil2.loc[:, 'L_p_mean'].mean()
    
    if lp1_avg>lp2_avg:
        df_fil1['fil_label_L_p'] = 'fil1';df_fil2['fil_label_L_p'] = 'fil2'
    else:
        df_temp_swap = df_fil2
        df_fil2 = df_fil1
        df_fil1 = df_temp_swap
        df_fil1['fil_label_L_p'] = 'fil1';df_fil2['fil_label_L_p'] = 'fil2'

        

    df_fil1['avg_arc_l_int'] = np.round(100*df_fil1['avg_arc_l_'],2).astype(int)
    
    df_fil2['avg_arc_l_int'] = np.round(100*df_fil2['avg_arc_l_'],2).astype(int)
    df_fil1_fil2 = pd.merge(df_fil1, df_fil2,suffixes=("_fil1","_fil2"),on=['avg_arc_l_int','file_label_1','file_label_0','file_label_2'],how = 'inner')
    return df_fil1_fil2
def find_MSD_metric_ori(df_rdf):
    '''
    input hte dfrdf from the methods class 
    does the computations for only half the Delta T
    return DF_MSD with radial MSD and also ang MSD
    
    '''
    df_msd = pd.DataFrame()    
    df_msd_arr = pd.DataFrame()        

    #df_msd["MSD_arr"] = df_msd['MSD_arr'].astype(object)
    #df_msd["ang_MSD_arr"] = df_msd['ang_MSD_arr'].astype(object)


    uni_arc_L = df_rdf['avg_arc_L'].unique()
    fil_labels = ['fil1','fil2']
    file_label_temp =df_rdf['file_label'][0]
    for la in uni_arc_L:
        for fl in fil_labels:
            mask_iso = (df_rdf['avg_arc_L']== la) & (df_rdf['fil_label']== fl)
            df_rdf_iso_arc = df_rdf[mask_iso ]
            end_end_dist_iso = df_rdf_iso_arc['end_end_dist'].to_numpy()
            tan_ang_iso =         df_rdf_iso_arc['tangential_angle_deg'].to_numpy()
            L_p_arr =         df_rdf_iso_arc['L_p'].to_numpy()
            
            L_p_R_sq_arr =         df_rdf_iso_arc['L_p_Rsq'].to_numpy()
            N_2 = int(len(L_p_arr)*0.5)

            L_p_arr_N2 = L_p_arr[:N_2]
            L_arr =  df_rdf_iso_arc['l_seg'].to_numpy()
            mask_good_L_P = L_p_R_sq_arr[:N_2] >=0.99
            
            #print(la,fl,len(L_p_arr),sum(mask_good_L_P),N_2)
            #loop to find the 
            for i in range(1,len(end_end_dist_iso)):

                deviation =          end_end_dist_iso[i:] -end_end_dist_iso[:-i]
                tan_ang_deviation =  tan_ang_iso[i:] - tan_ang_iso[:-i]

                L_p_deviation = L_p_arr[i:] - L_p_arr[:-i]

                deviation_2 = (deviation**2).astype(float)
                deviation_2 = np.round(deviation_2,2)
                MSD = np.mean(deviation_2);std_MSD = (np.std(deviation_2)/np.sqrt(len(deviation_2)-1))


                
                tan_ang_deviation_2= (tan_ang_deviation).astype(float)
                tan_ang_deviation_2 = np.round(tan_ang_deviation_2,4)
                tan_ang_MSD = np.mean(tan_ang_deviation_2);std_tan_ang_MSD = (np.std(tan_ang_deviation_2)/np.sqrt(len(tan_ang_deviation_2)-1))

                L_p_deviation_2 = (L_p_deviation**2).astype(float)
                L_p_deviation_2 = np.round(L_p_deviation_2,4)
                
                L_p_MSD = np.mean(L_p_deviation_2);std_L_p_MSD = (np.std(L_p_deviation_2)/np.sqrt(len(L_p_deviation_2)-1))

    
                #saving the data frame
                #saving the raw arrays      
                if i+1<0.5*len(end_end_dist_iso):
                    df_temp = pd.DataFrame({'MSD_arr':deviation,
                                            'tan_ang_MSD_arr':tan_ang_deviation,

                                            'L_p_deviation_arr':L_p_deviation,
                                            
                                            "data_label":f"{fl}_{la}",
                                            'deltaT':i,
                                            'log_deltaT':np.log10(i)}) 
                    
                    
                    df_msd_arr = pd.concat([df_msd_arr, df_temp], axis=0) 

                L_temp = np.mean(L_arr[:N_2]);L_p_temp = np.mean(L_p_arr_N2)
                N_df = len(df_msd)
                df_msd.loc[N_df,'file_label']  = file_label_temp
                df_msd.loc[N_df,'data_label'] = f"{fl}_{la}"

                df_msd.loc[N_df,'deltaT']  = i;df_msd.loc[N_df,'log_deltaT']  = np.log(i)
                
                #radial MSDs
                df_msd.loc[N_df,'MSD']  = MSD;df_msd.loc[N_df,'std_norm_MSD']  = std_MSD
                df_msd.loc[N_df,'std_MSD']  = np.std(deviation_2)
                df_msd.loc[N_df,'MSD_norm']  = (MSD*(L_p_temp**2))/(L_temp)
                df_msd.loc[N_df,'dT_MSD_norm']  = ((i*L_p_temp)/L_temp**4)
                

                #tangential angle MSD 
                df_msd.loc[N_df,'tan_ang_MSD']  = tan_ang_MSD;df_msd.loc[N_df,'std_norm_tan_ang_MSD']  = std_tan_ang_MSD
                df_msd.loc[N_df,'tan_ang_MSD_norm']  = tan_ang_MSD/L_p_temp
                df_msd.loc[N_df,'std_tan_ang_MSD']  = np.std(tan_ang_deviation_2)
                
                
                #LP MSD 
                df_msd.loc[N_df,'L_p_MSD']  = L_p_MSD;df_msd.loc[N_df,'std_norm_L_p_MSD']  = std_L_p_MSD
                df_msd.loc[N_df,'std_L_p_MSD']  = np.std(L_p_deviation_2)
                

                df_msd.loc[N_df,'mean_l_p']  = L_p_temp
                df_msd.loc[N_df,'mean_l_p_good']  = np.mean(L_p_arr_N2[mask_good_L_P])

                df_msd.loc[N_df,'std_l_p']  = np.std(L_p_arr[:N_2])
                df_msd.loc[N_df,'L_seg']  = L_temp
                
                
                

                df_msd.loc[N_df,'avg_arc_L']  = la;df_msd.loc[N_df,'fil_label']  = fl
    return(df_msd,df_msd_arr)

    
def find_angle(A,B,C):
    '''
    

    Parameters
    ----------
    A : 2 Dvector 
        DESCRIPTION.
    B : 2 Dvector 
        DESCRIPTION.
    C : 2 Dvector 
        DESCRIPTION.

    Returns
    -------
    bend angle between two vecotrs in degrees 
    '''
    BA = B -A
    CB = C - B
    cosine_angle = np.dot(BA, CB) / (np.linalg.norm(BA) * np.linalg.norm(CB))
    angle_bet_ABC = np.degrees(np.arccos(cosine_angle))
    
    return angle_bet_ABC


def findMiddle_arr(input_list):
    
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        return input_list[int(middle)]
    
def find_folder_name(path):
    return os.path.basename(os.path.dirname(path))
def get_arrays_rgn(N=11):
    ar1 = np.linspace(0,0.5,N)
    ar2 = np.linspace(0.5,1,N)
    rgn1 = [];rgn2 = []
    
    for i in range(N-1):
        rgn1.append([ar1[i],ar1[i+1]])
        rgn2.append([ar2[i],ar2[i+1]])
    rgn1 = np.array(rgn1)
    rgn2 = np.array(rgn2)

    return(rgn1,rgn2)