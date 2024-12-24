#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:32:57 2024

@author: yogehs
"""
import numpy as np 
import pandas as pd
from scipy.optimize import curve_fit

from scipy.ndimage.filters import gaussian_filter1d, uniform_filter1d

from scipy.stats import norm

from constants import *
from utils import *
class Filament():
    '''
    defining stuff for a single filament, at a given time
    
    i.e. one csv file output 
    
    pass the file path for the csv file 
    
    '''
    def __init__(self,file_path,num_segs=10,flip_corrd = False):
        self.file_path = file_path
        #self.slice_num = slice_num
        
        #resolution loaded from tiff file, px/nm
        #self.resolution = img_res
        self.num_segs = num_segs
        #unit in nm 
        self.X =[];self.Y =[]
        self.X_c =[];self.Y_c =[]
        self.flip_corrd =flip_corrd 
        #array of (x,y) positions 
        self.vector_arr =[]
        
        #point curvature 
        self.pt_curv = [];self.smooth_pt_curv =[]
        self.pt_curv_ys=[] ;self.smooth_pt_curv_ys=[]        

        self.mean_curv_seg = np.zeros(num_segs);self.std_curv_seg = np.zeros(num_segs)
        self.mean_curv_ys_seg = np.zeros(num_segs);self.std_curv_ys_seg = np.zeros(num_segs)
        
        #the bending angle from scaling the curvature with arc len in degress 
        #w.r.t x axis and in degress
        self.bending_ang_seg =np.zeros(num_segs)
        self.bending_ang_fil = 0
        
        #tangential angle w.r.t x axis from the derivative array
        
        self.tan_angle_arr = 0
        #arc lenght 
        self.arc_l = [0];
        
        self.arc_l_norm = [] #normalised within a frame 
        self.arc_l_norm_folder = [] #normalised withing a folder 
        
        self.max_arc_l = 0
        self.max_max_arc_l_folder = 0 # maximum of of the max arc l in all folders
        self.mean_arc_l_seg = np.zeros(num_segs)
        self.arc_l_seg = np.zeros(num_segs)
        #square of the dist 
        self.end_to_end_dist_2_raw =[0]
        self.end_to_end_dist_2_seg =np.zeros(num_segs)
        
        
        #window for moving arc lenght (normalised)
        self.arc_win_A = 0 ; self.arc_win_B = 1
        self.mask_win = []
        #binnin the filament 
        self.bin_edges = []
        #persistence_len
        self.L_p = 0;self.L_p_std = 0
        self.L_p_R_sq = 0;self.L_p_dist = 0
        #bending energy in units of KbT
        self.E_bend_arr = np.zeros(num_segs);self.E_bend_fil = 0
        #o_p dataframe
        self.df_output = pd.DataFrame()
    def init_frame_params(self):
      
        #array of (x,y) positions 
        self.vector_arr =[]
        
        #point curvature 
        
        #arc lenght 
        self.arc_l = [0];
    
        self.arc_l_norm = [] #normalised within a frame 
        self.arc_l_norm_folder = [] #normalised withing a folder 
        
        self.end_to_end_dist_2_raw =[]
  

        
    def read_data(self):
        df_kapp = pd.read_csv(self.file_path)
        
        self.df_output = pd.DataFrame(columns=df_kapp.columns)
        #self.find_bin_edges()

        self.pt_curv=df_kapp["Point Curvature (um-1)"].multiply( df_kapp["Point Curvature Sign"])
        self.pt_curv *= 0.001
        self.pt_curv = np.array(self.pt_curv)
        temp ,self.smooth_pt_curv = Rollavg(self.pt_curv,len(self.pt_curv),mode_roll)


        self.X = df_kapp['X-Coordinate (um)']
        self.Y = df_kapp['Y-Coordinate (um)']
        self.X *= 1000 ;  self.Y *=1000
        self.X= np.array(self.X);self.Y= np.array(self.Y)
        
        if self.flip_corrd :
            print("flipping cord")
            self.X = self.X[::-1]
            self.Y = self.Y[::-1]

        self.mask_win = np.array([True]*len(self.X))
        self.centering_X_Y()
    def centering_X_Y(self):
        #make edge 1 at 0,0
        
        xx = self.X- self.X[0]
        yy = self.Y- self.Y[0]
        #vectorise the Edge 1 and edge 2 
        
        V = np.array([xx[-1]-xx[0],yy[-1]-yy[0]])
        #find theta 
        self.theta_orient =np.angle(V[0]+1j *V[1])
        #define the roation matrix for this theta 
        rotMatrix = np.array([[np.cos(self.theta_orient), np.sin(self.theta_orient)], 
                                 [-np.sin(self.theta_orient),  np.cos(self.theta_orient)]])
        for i in range(len(xx)):
            temp = rotMatrix.dot(np.array([xx[i],yy[i]]))
            self.X_c.append(temp[0]);self.Y_c.append(temp[1])
        self.X_c = np.array(self.X_c);self.Y_c = np.array(self.Y_c)
    def load_data(self):
        self.read_data()
        self.prep_data()
        self.custom_curvature()
    def custom_curvature(self):
    
        dx_dt = np.gradient(self.X_c)
        dy_dt = np.gradient(self.Y_c)
        
        dx_dt_raw = np.gradient(self.X)
        dy_dt_raw = np.gradient(self.Y)
        
        self.tan_angle_arr = np.arctan2(dy_dt,dx_dt)
        ds_dt = np.sqrt(dx_dt * dx_dt + dy_dt * dy_dt)
        d2s_dt2 = np.gradient(ds_dt)
        d2x_dt2 = np.gradient(dx_dt)
        d2y_dt2 = np.gradient(dy_dt)
    
        self.pt_curv_ys = (d2x_dt2 * dy_dt - dx_dt * d2y_dt2)/ (dx_dt**2 + dy_dt**2)**1.5

    def prep_data(self):
        '''vetorise hte coordinates
        find the arc length 
        find the bin edges 
        find end to end distance raw'''
        self.vector_arr =[]
        for i in range(len(self.X)):
            #self.vector_arr.append(np.array([self.X[i],self.Y[i]]))
            self.vector_arr.append(np.array([self.X_c[i],self.Y_c[i]]))

        self.find_arc_L()        
        self.find_bin_edges()

        self.find_end_end_dist_raw()
        self.fit_L_p()
        self.dist_L_P()
        self.find_avg_cruvature()
    def mask_data(self):
        self.X= self.X[self.mask_win];self.Y= self.Y[self.mask_win]
        self.X_c= self.X_c[self.mask_win];self.Y_c= self.Y_c[self.mask_win]

        self.pt_curv = self.pt_curv[self.mask_win]
        
        self.vector_arr =[]

        for i in range(len(self.pt_curv)):
            #self.vector_arr.append(np.array([self.X[i],self.Y[i]]))
            self.vector_arr.append(np.array([self.X_c[i],self.Y_c[i]]))

        self.arc_l = self.arc_l[self.mask_win] ;self.arc_l_norm = self.arc_l_norm [self.mask_win]
        #self.max_arc_l = self.arc_l[-1]-self.arc_l[0]
        #self.arc_l = self.arc_l[self.mask_win] -self.arc_l[self.mask_win][0]
        #self.arc_l_norm = self.arc_l_norm [self.mask_win]-self.arc_l_norm [self.mask_win][0]

        #self.find_arc_L()
        self.find_bin_edges()
        self.find_end_end_dist_raw()
        #self.fit_L_p()
        self.find_avg_cruvature()


    def find_bin_edges(self):
        self.bin_edges =[]
        #print(min(self.arc_l_norm),max(self.arc_l_norm))
        bins_temp = np.linspace(min(self.arc_l_norm),max(self.arc_l_norm),self.num_segs+1)

        for i in range(self.num_segs):
            self.bin_edges.append((bins_temp[i],bins_temp[i+1]))
            
            
    def fit_L_p(self):
        #arc_l_win = self.max_arc_l*(self.arc_l_norm-self.arc_l_norm[0])
        
        arc_l_win_seg = self.arc_l_seg
        #pars, cov = curve_fit(self.wlc_LP_fit,arc_l_win ,self.end_to_end_dist_2_raw,p0 =[20])
        pars, cov = curve_fit(self.wlc_LP_fit,arc_l_win_seg ,self.end_to_end_dist_2_seg,p0 =[2])

        stdevs = np.sqrt(np.diag(cov))
        #print("fit parameters : ",pars,"  ",stdevs)
        self.L_p =np.round(pars[0],3) ;self.L_p_std = stdevs
        end_end_dist_2_fit = self.wlc_LP_fit(arc_l_win_seg,self.L_p)

        # residual sum of squares
        ss_res = np.sum((self.end_to_end_dist_2_seg - end_end_dist_2_fit) ** 2)
        
        # total sum of squares
        ss_tot = np.sum((self.end_to_end_dist_2_seg - np.mean(self.end_to_end_dist_2_seg)) ** 2)
        
        # r-squared
        self.L_p_R_sq =np.round( 1 - (ss_res / ss_tot),3)
        self.dist_L_P()
        #print(avg_arcl_1_arr)
        
    def dist_L_P(self):
        '''
        in this module I will find the LP from the distribution of hte curvature 
        '''
        mean_fit,std_fit=norm.fit(self.pt_curv-np.mean(self.pt_curv))
        self.L_p_dist =np.round(1.0/(std_fit**2),3)

    def wlc_LP_fit(self,L,L_p):
        end_end_dist_2 = 4*L*L_p  - 8*(L_p**2)*(1-np.exp(-L/(2*L_p)))
        return end_end_dist_2
    
    
    def find_E_bend(self):
        '''
        calculated in the unites of KbT, felxural rigidity = L_p * KBT 
        E_bend = k_f  * segment length * avg_curvatrue**2/(2)
        '''
        segment_len= np.array([0]+list(self.arc_l_seg))
        segment_len = segment_len[1:] - segment_len[:-1]
        for i in range(self.num_segs):

            self.E_bend_arr[i] = 0.5* self.L_p* segment_len[i] * (self.mean_curv_seg[i]**2)
        self.E_bend_fil = np.round(np.sum(self.E_bend_arr),4)
        #self.E_bend_fil = 0.5* self.L_p* (self.arc_l_seg[-1]-self.arc_l_seg[0]) * (np.mean(self.pt_curv_ys)**2)

    def find_end_end_dist_raw(self):

        self.end_to_end_dist_2_raw =[]
        pt_a = self.vector_arr[0]

        N = len(self.vector_arr)
        for i in range(N):
            pt_b= self.vector_arr[i]
            end_to_end_dist_a_b = np.linalg.norm(pt_a-pt_b)

            self.end_to_end_dist_2_raw.append(end_to_end_dist_a_b**2)
        self.end_to_end_dist_2_raw = np.array(self.end_to_end_dist_2_raw)


    def find_avg_cruvature(self):
        '''
        finds and assings the follwoing 
        - avg curvature for all  segment in a filament
        - STD curvature for all  segment in a filament 
        - average arc lenght for all  segment in a filament
        -finds the overall bend angle of the filament w.r.t x axis 
        - bend angle of each segment w.r.t x axis 
        '''
        dL = self.arc_l[-1]- self.arc_l[0]
        self.bending_ang_fil = np.rad2deg(dL*np.mean(self.pt_curv))
        
        le_0 ,he_0 = self.bin_edges[0]
        for i in range(self.num_segs):

           le,he = self.bin_edges[i]
           bin_mask = (le<=self.arc_l_norm) & (self.arc_l_norm<he)
           temp_curv_1 = self.pt_curv[bin_mask] 
           temp_end_end_dist = self.end_to_end_dist_2_raw[bin_mask]
           #print(temp_end_end_dist)
           self.mean_arc_l_seg[i] = (le+he)*0.5
           self.arc_l_seg[i] = (he-le_0)*self.max_arc_l

           self.end_to_end_dist_2_seg[i] = temp_end_end_dist[-1]
           self.mean_curv_seg[i] = np.mean(temp_curv_1) 
           self.std_curv_seg[i] = np.std(temp_curv_1) 
           self.bending_ang_seg[i] = np.rad2deg(self.arc_l_seg[i]*np.mean(temp_curv_1) )

        #self.find_E_bend()
    def find_arc_L(self):
        self.arc_l = [];self.arc_l_norm = []

        arc_l_pre_sum = []
        for i in range(len(self.vector_arr)-1):
            arc_l_pre_sum.append(np.linalg.norm(self.vector_arr[i+1]-self.vector_arr[i]))
            
        for j in range(len(arc_l_pre_sum)):
            self.arc_l.append(np.sum(arc_l_pre_sum[:j]))
        self.arc_l.append(np.sum(arc_l_pre_sum[:]))
            
        self.arc_l = np.array(self.arc_l)
        self.arc_l_norm  = np.array(self.arc_l/np.max(self.arc_l))
        self.max_arc_l = np.max(self.arc_l)
     
    def cal_flex_rigid(self,segment_len):
        
        self.flex_rigid = (2*k_bT_nm)/(segment_len*self.mean_curv_seg**2)