#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:39:19 2024

@author: yogehs
"""
import shutil

from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QPushButton
import pyqtgraph as pg
import os
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.optimize import curve_fit


from scipy.stats import norm

#from class_fil import Filament
from utils import *
from constants import *
from PIL import Image
from PIL.TiffTags import TAGS
import copy
from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QPushButton,QFileDialog


pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

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
from pyqtgraph.Qt import QtWidgets, QtCore
class Filament_mech_gui():
    """ takes a list of file paths to make an filament onject
    """
    def __init__(self,folder_path,file_label,num_segs=5 ,del_bool =False,flip_corrd =False):
        super().__init__()
        self.app = pg.mkQApp("Filament mech") 
        
        self.mw = QtWidgets.QMainWindow()
        self.view = pg.GraphicsLayoutWidget()
        self.mw.setCentralWidget(self.view)

        #frame stuff
        temp_file_arr = sorted([f for f in os.listdir(folder_path) if f.endswith('.csv')])
        temp_file_arr = sorted(temp_file_arr, key=extract_number)

        self.file_frame_arr = [folder_path+f for f in temp_file_arr]

        self.num_frames = len(self.file_frame_arr)
        self.frame_num_arr = np.arange(1,self.num_frames+1,1)

        self.i_frame = 0
        self.filaments_arr = []
        self.resolution = 1
        #window for moving arc lenght (normalised)
        self.arc_win_A = 0 ; self.arc_win_B = 0.5
        self.mask_win_gui = []
        self.masked_filament = Filament(self.file_frame_arr[0])
        self.df_output = pd.DataFrame(columns=['arc_L_segment','end_end_dist_2_segment'],dtype=object)
        #self.df_output['arc_L_segment'] = self.df_output['arc_L_segment'].astype('object')
        #self.df_output['end_end_dist_2_segment'] = self.df_output['end_end_dist_2_segment'].astype('object')
        self.file_label = file_label
        self.init_all_frames(num_segs,flip_corrd)
        #os.makedirs(folder_path+"/results",exist_ok = True)
        self.save_path = folder_path+"/results/"
        if del_bool:shutil.rmtree(self.save_path,ignore_errors=True )
        os.makedirs(folder_path+"/results",exist_ok = True)

        

        
    def load_gui(self):
        self.create_app()
        self.load_frame()
        
    def init_resolution(self):
        '''
        this func will find the resolution of image from the tiff file 
        '''
        tiff_files = sorted([f for f in os.listdir(self.folder_path) if f.endswith('.tif')])
        img = Image.open(self.folder_path+'/'+tiff_files[0])
        meta_dict = {TAGS[key] : img.tag[key] for key in img.tag_v2}
        #print(meta_dict['XResolution'][0])
        self.resolution = meta_dict['XResolution'][0][0]/meta_dict['XResolution'][0][1]
        
    def add_button(self,btn):
        proxy = QtWidgets.QGraphicsProxyWidget()
        proxy.setWidget(btn)
    
        return(proxy)
        #win.addItem(proxy)
        
    
    def btn_next_frame(self):
        self.i_frame +=1
        #self.save_sess_log()

        if self.i_frame==self.num_frames:
            self.file_label_itm.setText("C'est fini")
            self.AFM_plt.clear();self.L_p_plt.clear();self.X_Y_plt.clear()

        else:
            #self.onButtonClicked_t_rup_redo()
            self.load_frame()
    def btn_prev_frame(self):
        self.i_frame -=1
        if self.i_frame<0:
            self.file_label_itm.setText("Index 0, no more file post")
            self.AFM_plt.clear();self.L_p_plt.clear();self.X_Y_plt.clear()

        else:
            #self.onButtonClicked_t_rup_redo()
            #self.pop_sess_log()
            self.load_frame()
    def init_all_frames(self,num_segs,flip_corrd):
        #define the list of objects
        for i in range(self.num_frames):
            i_filament = Filament(self.file_frame_arr[i],num_segs,flip_corrd)
            i_filament.load_data()

            self.filaments_arr.append(i_filament)
        
    def load_frame(self):
        self.X_Y_plt.clear();        self.AFM_plt.clear()

        self.X_Y_plt.addItem(self.arc_len_win)

        self.X_Y_plt.addItem(self.arc_len_win2)
        self.moved_arc_len_win()

        self.file_label_itm.setText("frame number : "+str(self.frame_num_arr[self.i_frame]))

        i_filament = self.filaments_arr[self.i_frame]
        X = i_filament.X ; Y = i_filament.Y  
        arc_l = i_filament.arc_l_norm ; curv = i_filament.pt_curv
        curv_sm = i_filament.smooth_pt_curv
        
        #adding the X,Y pos
        pen_prop = pg.mkPen('y',cosmetic = False, width=3 )#style=QtCore.Qt.DashLine)
        self.X_Y_plt.enableAutoRange(axis='x')
        self.X_Y_plt.setAutoVisible(x=True)
        #self.X_Y_plt.plot(X,Y)
        scatter_plt = pg.ScatterPlotItem()
        scatter_plt.setData(x = arc_l, y = curv,pen = pen_prop)
        frame_label= 'frame num'+str(self.frame_num_arr[self.i_frame])+"_resolution(px/nm)_"+str(self.resolution)
        self.X_Y_plt.addItem(scatter_plt)

        #self.X_Y_plt.plot(arc_l,curv,pen = pen_prop,connect='finite',symbol = 'star',name=frame_label)
        #self.X_Y_plt.setPen('y', width=3, style=QtCore.Qt.DashLine)          ## Make a dashed yellow line 2px wide
        self.X_Y_plt.addLegend(frame=False,QSizeF = 10)

        #self.arc_len_win.setRegion()
        
    def find_MSD(self):
        '''
        find the mean squre flutations in 
        - position x_Y
        - end to end dist 
        -curvature 
        all along the arc lenght acroos all the frames, to capture the dyanmics 
        '''
        
    def prep_frames(self):
        
        #prepare teh relavante data from all frames 
        print("lolo")    
    def onButtonClicked_L_p(self):
        #compute the L_p
        
        
        print('lol')
        
        
    def key_pressed_YS(self,ev):
        if keys_mapping[ev.key()]=='D':
            self.btn_next_frame()
        elif keys_mapping[ev.key()]=='A':
            self.btn_prev_frame()

        elif keys_mapping[ev.key()]=='S':
            self.save_sess()

    def create_app(self):
        
        #initial setup 
        self.mw.resize(1000,800)
        self.view = KeyPressWindow(show=True)
        
        self.mw.setCentralWidget(self.view)
        
        #self.view.sigKeyPress.connect(lambda event: print(keys_mapping[event.key()]))
        self.view.sigKeyPress.connect(lambda event: self.key_pressed_YS(event))

        self.mw.show()
        self.mw.setWindowTitle('Filament Mech')
        #creating buttons  
        self.compute_L_p_btn = QtWidgets.QPushButton("get L_p")
        self.prev_frame_btn = QtWidgets.QPushButton("Prev frame")
        self.next_frame_btn = QtWidgets.QPushButton("Next frame")
        self.output_label_yell = pg.LabelItem("Persitence len: ")
        self.output_label_green = pg.LabelItem("Persitence len: ")
        self.file_label_itm = pg.LabelItem("Persitence len: ")

        #self.output_label_yell.setStyleSheet("background-color: yellow")
        #self.output_label_green.setStyleSheet("background-color: green")

        self.save_btn = QtWidgets.QPushButton("Save value")

        #connecting the buttons
        self.prev_frame_btn.clicked.connect(self.btn_prev_frame)
        self.next_frame_btn.clicked.connect(self.btn_next_frame)
        self.compute_L_p_btn.clicked.connect(self.onButtonClicked_L_p)
        self.save_btn.clicked.connect(self.save_sess)
        #adding plots 
        self.AFM_plt =self.view.addPlot() 

        self.X_Y_plt =self.view.addPlot() 
        self.L_p_plt =self.view.addPlot() 

        self.L_p_plt.setLabel("bottom", "Arc length nm")
        self.L_p_plt.setLabel("left", "End -end distance - square  nm*2")
        self.X_Y_plt.setLabel("bottom", "normalised Arc length nm")
        self.X_Y_plt.setLabel("left", "Pt curvature  nm-1")
        self.X_Y_plt.setXRange(0,1)

        #self.L_p_plt.setYLink(self.X_Y_plt)
        #self.L_p_plt.setXLink(self.X_Y_plt)


        ## create fview
        self.view.nextRow()

        self.b_pf = self.view.addItem(self.add_button(self.prev_frame_btn))
        self.view.addItem(self.file_label_itm) 

        self.b_nf = self.view.addItem(self.add_button(self.next_frame_btn))
        self.view.nextRow()

        self.view.addItem(self.AFM_plt)
        self.view.addItem(self.X_Y_plt)
        self.view.addItem(self.L_p_plt)

        self.view.nextRow()

        self.view.addItem(self.output_label_yell) 
        self.view.addItem(self.output_label_green) 

        self.btn_save = self.view.addItem(self.add_button(self.save_btn));


        
        #adding arc_len win
        
        self.arc_len_win = pg.LinearRegionItem(values=[0,0.5],pen = 'yellow',brush = 'yellow')
        self.arc_len_win.setBounds([0,1])
        #this is done to stop the lines from moving and changing the windowsize
        self.arc_len_win.lines[0].setMovable = False
        #self.otr_label = pg.InfLineLabel(self.otr_win.lines[0], "OTR WIN ",color='black', position=0.95, rotateAxis=(0,-1), anchor=(1, 1))

        self.arc_len_win.sigRegionChanged.connect(self.moved_arc_len_win)
        
        self.arc_len_win2 = pg.LinearRegionItem(values=[0.5,1],pen = 'green',brush = 'green')
        self.arc_len_win2.setBounds([0,1])
        #this is done to stop the lines from moving and changing the windowsize
        self.arc_len_win2.lines[0].setMovable = False
        #self.otr_label = pg.InfLineLabel(self.otr_win.lines[0], "OTR WIN ",color='black', position=0.95, rotateAxis=(0,-1), anchor=(1, 1))

        self.arc_len_win2.sigRegionChanged.connect(self.moved_arc_len_win)
    def update_arc_len_win(self,rgn):
        ''''ipdates and stores the legs of win 
        the mask is for {win_a,win_b)
                         '''
                         
        i_filament = self.filaments_arr[self.i_frame]

        #i_filament.load_data()
        self.arc_win_A = rgn[0] ; self.arc_win_B = rgn[1]
        self.mask_win_gui = (i_filament.arc_l_norm >= self.arc_win_A)& (i_filament.arc_l_norm < self.arc_win_B)
        self.masked_filament = copy.deepcopy(i_filament)
        self.masked_filament.mask_win = self.mask_win_gui
        self.masked_filament.arc_win_A = rgn[0] ; self.masked_filament.arc_win_B = rgn[1]
        self.masked_filament.mask_data()
        self.masked_filament.fit_L_p()
        self.masked_filament.dist_L_P()

        self.masked_filament.find_E_bend()
        
    def update_afm_plot(self,pen_prop):
        i_filament =  self.filaments_arr[self.i_frame]
     

        xx = i_filament.X_c ; yy = i_filament.Y_c

        #self.AFM_plt.plot()
        scatter_plt = pg.ScatterPlotItem()
        scatter_plt1 = pg.ScatterPlotItem()

        scatter_plt.setData(x = xx[self.mask_win_gui], y = yy[self.mask_win_gui],pen = pen_prop)
        scatter_plt1.setData(x = xx[~self.mask_win_gui], y = yy[~self.mask_win_gui])
        pen_prop_neg =  pg.mkPen('gray' , width=3)

        scatter_plt1.setPen(pen_prop_neg)
        scatter_plt1.setOpacity(0.2)
        self.AFM_plt.addItem(scatter_plt)
        self.AFM_plt.addItem(scatter_plt1)


        
    def moved_arc_len_win(self):
        self.L_p_plt.clear()
        rgn = self.arc_len_win.getRegion()
        rgn2 = self.arc_len_win2.getRegion()
        pen_prop =  pg.mkPen('yellow' , width=3)

        self.update_arc_len_win(rgn)
        self.update_afm_plot(pen_prop)

        i_filament =  self.masked_filament      
        L_p_yellow = i_filament.L_p;R_sq_yell = i_filament.L_p_R_sq
        E_end_avg_yel = i_filament.E_bend_fil
        arc_l = i_filament.arc_l_norm ; curv = i_filament.pt_curv
        curv_sm = i_filament.smooth_pt_curv
        end_end_dist_sq = i_filament.end_to_end_dist_2_raw
        arc_l_seg = i_filament.arc_l_seg;end_end_seg = i_filament.end_to_end_dist_2_seg
        
        arc_l_win = i_filament.max_arc_l*(i_filament.arc_l_norm-i_filament.arc_l_norm[0])

        fitted_end_end_dist_sq = i_filament.wlc_LP_fit(arc_l_seg,L_p_yellow)

        #self.L_p_plt.plot(arc_l,curv_sm,pen = pen_prop,connect='finite',symbol = 'star')
        scatter_plt = pg.ScatterPlotItem()
        scatter_plt.setData(x = arc_l_seg, y = end_end_seg,pen = pen_prop)
        pen_yellow_fit = pg.mkPen('black', width=3, style=QtCore.Qt.DashLine)  
        self.L_p_plt.plot(arc_l_seg,fitted_end_end_dist_sq,pen = pen_yellow_fit,label=f'Lp_yel{L_p_yellow}')
        
        self.L_p_plt.addItem(scatter_plt)
        self.output_label_yell.setText(f"L_p_yel ={L_p_yellow} , R_sq_yel = {R_sq_yell}, E_B ={E_end_avg_yel}",color = 'orange')


        #regeion 2
        pen_prop =  pg.mkPen('green',cosmetic = True, width=3, style=QtCore.Qt.DashLine)

        self.update_arc_len_win(rgn2)
        self.update_afm_plot(pen_prop)

        i_filament =  self.masked_filament  
        L_p_green = i_filament.L_p;R_sq_gr = i_filament.L_p_R_sq
        E_end_avg_green = i_filament.E_bend_fil

        arc_l = i_filament.arc_l_norm ; curv = i_filament.pt_curv
        curv_sm = i_filament.smooth_pt_curv
        end_end_dist_sq = i_filament.end_to_end_dist_2_raw
        arc_l_seg = i_filament.arc_l_seg;end_end_seg = i_filament.end_to_end_dist_2_seg
        
        
        fitted_end_end_dist_sq = i_filament.wlc_LP_fit(arc_l_seg,L_p_green)
        
        pen_green_fit = pg.mkPen('green', width=3, style=QtCore.Qt.DashLine)  
        scatter_plt_2 = pg.ScatterPlotItem()

        scatter_plt_2.setData(x = arc_l_seg, y = end_end_seg,pen = pen_prop)
        self.L_p_plt.plot(arc_l_seg,fitted_end_end_dist_sq,pen = pen_green_fit,label=f'Lp_gr{L_p_green}')
        
        self.L_p_plt.addItem(scatter_plt_2)

        self.L_p_plt.addLegend()
        self.output_label_green.setText(f"L_p_gr ={L_p_green} , R_sq_gr = {R_sq_gr}, E_B ={E_end_avg_green}",color = 'green')


    def save_sess(self):
        rgn = self.arc_len_win.getRegion()
        rgn2 = self.arc_len_win2.getRegion()

        self.update_arc_len_win(rgn);self.save_rgn()

        self.update_arc_len_win(rgn2);self.save_rgn()
        chime.success()

        
    def save_rgn(self,csv_file_name = "mech_analysis_gui.csv",bool_export_fig=True):
        N_output = len(self.df_output)
        temp_avg_arc_len = np.round(0.5*(self.masked_filament.arc_win_B +self.masked_filament.arc_win_A),3)
        if temp_avg_arc_len<=0.5:temp_fil_label = 'fil1'
        else:temp_fil_label = 'fil2'
        
        #print(N_trouble_sess,self.log_afs_df.loc[self.i_file])
        self.df_output.loc[N_output,'file_label']  = self.file_label
        self.df_output.loc[N_output,'file_path'] = self.file_frame_arr[self.i_frame]
        self.df_output.loc[N_output,'frame_numer'] = self.i_frame
        self.df_output.loc[N_output,'max_arc_l'] = self.filaments_arr[self.i_frame].max_arc_l

        self.df_output.loc[N_output,'Arc_win_rgn_A'] = self.masked_filament.arc_win_A
        self.df_output.loc[N_output,'Arc_win_rgn_B'] = self.masked_filament.arc_win_B


        self.df_output.loc[N_output,'avg_arc_l_win_fil_1'] = temp_avg_arc_len
        self.df_output.loc[N_output,'avg_arc_l_win_fil_2'] =1- temp_avg_arc_len
        
        self.df_output.loc[N_output,'fil_label'] =temp_fil_label

        self.df_output.loc[N_output,'Arc_win_size'] = np.round(self.masked_filament.arc_win_B -self.masked_filament.arc_win_A,3)
        self.df_output['max_arc_L_segment'] = np.max(self.masked_filament.arc_l_seg) 

        self.df_output['arc_L_segment'] = self.df_output['arc_L_segment'].astype(object) 
        self.df_output['end_end_dist_2_segment'] = self.df_output['end_end_dist_2_segment'].astype(object) 

        self.df_output.loc[N_output,'arc_L_segment'] = np.array2string(self.masked_filament.arc_l_seg, formatter={'float_kind':lambda x: "%.2f" % x},separator=',') 
        self.df_output.loc[N_output,'end_end_dist_2_segment'] = np.array2string(self.masked_filament.end_to_end_dist_2_seg, formatter={'float_kind':lambda x: "%.2f" % x},separator=',') 
        self.df_output.loc[N_output,'E_bend_arr_kBT'] = np.array2string(self.masked_filament.E_bend_arr, formatter={'float_kind':lambda x: "%.2f" % x},separator=',') 
        self.df_output.loc[N_output,'Mean_Curvature_segment'] = np.array2string(self.masked_filament.mean_curv_seg, formatter={'float_kind':lambda x: "%.2f" % x},separator=',') 
        self.df_output.loc[N_output,'Std_Curvature_segment'] = np.array2string(self.masked_filament.std_curv_seg, formatter={'float_kind':lambda x: "%.2f" % x},separator=',') 

        
        
        self.df_output.loc[N_output,'E_bend_filament_kBT'] = self.masked_filament.E_bend_fil

        #self.df_output.loc[N_output,'end_end_dist_2_segment'] = self.masked_filament.end_to_end_dist_2_seg
       
        self.df_output.loc[N_output,'L_p'] = self.masked_filament.L_p
        self.df_output.loc[N_output,'L_p_std'] = self.masked_filament.L_p_std
        self.df_output.loc[N_output,'R_sq_fit'] = self.masked_filament.L_p_R_sq
        self.df_output.loc[N_output,'L_P_dist'] = self.masked_filament.L_p_dist
        self.df_output.loc[N_output,'mean_pt_curv_for_dist'] = np.mean(self.masked_filament.pt_curv)

        self.df_output.to_csv(self.save_path+csv_file_name)
        if bool_export_fig:
            exporter_L_p = pg.exporters.ImageExporter(self.L_p_plt)
            exporter_AFM= pg.exporters.ImageExporter(self.AFM_plt)
            exporter_curv= pg.exporters.ImageExporter(self.X_Y_plt)
    
            exporter_L_p.export(self.save_path+f'{N_output}_L_P_plot_in_df_output.png')
            exporter_AFM.export(self.save_path+f'{N_output}_AFM_plot_in_df_output.png')
            exporter_curv.export(self.save_path+f'{N_output}_Curvat_plot_in_df_output_fr_num_{self.i_frame}.png')
    def roll_win_L_p_1_px(self,del_arc_l = 0.1):
        '''
        compute and save LP from wlc and dist  
        in two region simulateneously, x,0.1+x and 1-(0.1+x), 1-x
        fil1 and fil2 flipped  
        doest this over a rollwing window moved at 1% distance 
        
        saves the fil
        '''
        color_arr = ['yellow','green']
        for i in range(self.num_frames):
            self.i_frame = i
            #rolls over the arc len in a pixel length
            roll_win_size = 1.0/self.filaments_arr[i].max_arc_l
            temp_win_a = 0.0 
            temp_win_b = temp_win_a+del_arc_l
            win_num = 0
            while temp_win_a+del_arc_l<=0.5:
                fig, (ax1, ax2) = plt.subplots(1, 2)

                rgn1 = [temp_win_a,temp_win_b]
                rgn2 = [1- temp_win_b,1- temp_win_a]
                rgn_arr = [rgn1,rgn2]
                #plt.clf()
                for i_rgn in range(2):
                    
                    temp_rgn = rgn_arr[i_rgn]

                    self.update_arc_len_win(temp_rgn)
                    
                    i_fil_mask =  self.masked_filament  
                    L_p_green = i_fil_mask.L_p;R_sq_gr = i_fil_mask.L_p_R_sq
                    E_B = i_fil_mask.E_bend_fil;pt_curv_0mean = i_fil_mask.pt_curv - np.mean(i_fil_mask.pt_curv )
                    
                    arc_l_seg = i_fil_mask.arc_l_seg;end_end_seg = i_fil_mask.end_to_end_dist_2_seg
                    
                    
                    fitted_end_end_dist_sq = i_fil_mask.wlc_LP_fit(arc_l_seg,L_p_green)
                    ax1.scatter(arc_l_seg,end_end_seg,c = color_arr[i_rgn], label = f"[{temp_rgn[0]},{temp_rgn[1]}]")
                    ax1.plot(arc_l_seg,fitted_end_end_dist_sq,color = color_arr[i_rgn],linestyle ='dashed',label = f"L_p_{L_p_green}_R_sq_{R_sq_gr}_E_b_{E_B}")
                    
                    
                    ax2.hist(pt_curv_0mean,density =True,color = color_arr[i_rgn],histtype = "step")
                    xmin, xmax = plt.xlim()
                    x = np.linspace(xmin, xmax, 100)
                    mean_fit,std_fit=norm.fit(pt_curv_0mean)
                    
                    y = norm.pdf(x, mean_fit, std_fit)
                    ax2.plot(x, y,'--',color = color_arr[i_rgn],linewidth = 4,label = f"L_p_dist{i_fil_mask.L_p_dist} ")
                    
                    
                    self.save_rgn(f"roll_win_auto_1_px.csv",False)


                fig.suptitle(f"Frame_number_{i}_win_number_{win_num}")
                ax1.set_title("WLC");ax1.legend()
                ax2.set_title("LP from dist");ax2.legend()

                
                ax1.set_xlabel("Arc length nm");ax1.set_ylabel("End to end distance square nm^2")
                ax2.set_xlabel("mean corrected kappa")
                #fig.show()
                fig.savefig(self.save_path+f"frame_{i}_win_num_{win_num}_win_size_{del_arc_l}_L_P_fit_pixel.png")
                temp_win_a +=roll_win_size
                temp_win_b = temp_win_a+del_arc_l;win_num+=1
    def L_p_plot_poster(self,temp_win_a,ax,i_c,del_arc_l = 0.1):
        color_arr = ['deepskyblue','hotpink'];marker_arr= ['v','s']
        config_label = ["S","C"]
        self.i_frame = 0
        temp_win_b = temp_win_a+del_arc_l
        rgn1 = [temp_win_a,temp_win_b]

        self.update_arc_len_win(rgn1)
        
        i_fil_mask =  self.masked_filament  
        L_p_green = i_fil_mask.L_p;R_sq_gr = i_fil_mask.L_p_R_sq
        E_B = i_fil_mask.E_bend_fil;pt_curv_0mean = i_fil_mask.pt_curv - np.mean(i_fil_mask.pt_curv )
        
        arc_l_seg = i_fil_mask.arc_l_seg;end_end_seg = i_fil_mask.end_to_end_dist_2_seg
        
        
        fitted_end_end_dist_sq = i_fil_mask.wlc_LP_fit(arc_l_seg,L_p_green)
        ax.scatter(arc_l_seg,end_end_seg,c = color_arr[i_c],marker = marker_arr[i_c],s = 250)
        ax.plot(arc_l_seg,fitted_end_end_dist_sq,color = color_arr[i_c],linestyle ='dashed',label = f"{config_label[i_c]} L_p {L_p_green:.2f}",linewidth = 3.5,alpha =0.8)
            
    def segment_diag_plot(self,ax,del_arc_l = 0.1):
        color_arr = ['deepskyblue','hotpink'];marker_arr= ['v','s']
        config_label = ["S","C"]
        self.i_frame = 0
        
        lin_arr = np.linspace(0,1,int(np.ceil(1/(del_arc_l))), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0
           
        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+del_arc_l

            self.update_arc_len_win([temp_win_a,temp_win_b])
            
            i_fil_mask =  self.masked_filament 
            xx = i_fil_mask.X
            yy = i_fil_mask.Y
            
            
            ax.scatter(xx,yy,s = 50)
                
                       
              
    def roll_win_L_p_1_per(self,del_arc_l = 0.1,roll_size = 0.1):
        '''
        compute and save LP from wlc and dist  
        in two region simulateneously, x,0.1+x and 1-(0.1+x), 1-x
        fil1 and fil2 flipped  
        doest this over a rollwing window moved at 10% of win size (del_arc_l)
        
        saves the fil
        '''
        color_arr = ['deepskyblue','hotpink'];marker_arr= ['v','s']
        for i in range(self.num_frames):
            self.i_frame = i
            #rolls over the arc len in a percent
            roll_win_size = roll_size*del_arc_l
            temp_win_a = 0.0 
            temp_win_b = temp_win_a+del_arc_l
            win_num = 0
            while temp_win_a+del_arc_l<=0.5:
                #fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5));plt.cla()
                plt.figure(i,figsize=(7,9));plt.clf()
                rgn1 = [temp_win_a,temp_win_b]
                rgn2 = [1- temp_win_b,1- temp_win_a]
                rgn_arr = [rgn1,rgn2]
                #plt.clf()
                for i_rgn in range(2):
                    temp_rgn = rgn_arr[i_rgn]

                    self.update_arc_len_win(temp_rgn)
                    
                    i_fil_mask =  self.masked_filament  
                    L_p_green = i_fil_mask.L_p;R_sq_gr = i_fil_mask.L_p_R_sq
                    E_B = i_fil_mask.E_bend_fil;pt_curv_0mean = i_fil_mask.pt_curv - np.mean(i_fil_mask.pt_curv )
                    
                    arc_l_seg = i_fil_mask.arc_l_seg;end_end_seg = i_fil_mask.end_to_end_dist_2_seg
                    
                    
                    fitted_end_end_dist_sq = i_fil_mask.wlc_LP_fit(arc_l_seg,L_p_green)
                    plt.scatter(arc_l_seg,end_end_seg,c = color_arr[i_rgn],marker = marker_arr[i_rgn],s = 250)
                    plt.plot(arc_l_seg,fitted_end_end_dist_sq,color = color_arr[i_rgn],linestyle ='dashed',label = f"Segemnt Dimer_{i_rgn+1} L_p {L_p_green:.1f}",linewidth = 3.5,alpha =0.8)
                    
                    '''
                    ax2.hist(pt_curv_0mean,density =True,color = color_arr[i_rgn],histtype = "step")
                    xmin, xmax = plt.xlim()
                    x = np.linspace(xmin, xmax, 100)
                    mean_fit,std_fit=norm.fit(pt_curv_0mean)
                    
                    y = norm.pdf(x, mean_fit, std_fit)
                    ax2.plot(x, y,'--',color = color_arr[i_rgn],linewidth = 4,label = f"L_p_dist{i_fil_mask.L_p_dist} ")
                    '''
                    
                    self.save_rgn(f"roll_win_auto_{100*roll_size:.2f}_percent.csv",False)


                #fig.suptitle(f"Frame_number_{i}_win_number_{win_num}")
                #ax1.set_title("WLC");ax1.legend()
                #ax2.set_title("LP from dist");ax2.legend()

                
                #ax1.set_xlabel("Arc length nm");ax1.set_ylabel("End to end distance square nm^2")
                #ax2.set_xlabel("mean corrected kappa")
                #fig.show()
                #plt.legend(frameon= False);plt.show()
                
                #plt.savefig(self.save_path+f"poster_plot.png");break
                temp_win_a +=roll_win_size
                temp_win_b = temp_win_a+del_arc_l;win_num+=1     
    def roll_win_L_p_1_per_1rgn_cont(self,del_arc_l = 0.1,del_arc_l_nm = None,df_renorm = None):
        '''
        del_arc_l_nm : window size in nm
        del_arc_l: window size in normalised arc length
        compute and save LP from wlc and dist  
        in one region simulateneously, x,0.1+x 
        fil1 
        doest this over a rollwing window moved at 10% of win size (del_arc_l)
        
        saves the fil
        '''

        roll_size = 0.1
        color_arr = ['deepskyblue','hotpink'];marker_arr= ['v','s']
        for i in range(self.num_frames):
            self.i_frame = i
            if df_renorm is None:
                start_AL = 0;end_AL =1
            else:
                start_AL = df_renorm['norm_l_p_start'].to_numpy()[i]
                end_AL = df_renorm['norm_l_p_end'].to_numpy()[i]
            if del_arc_l_nm is not None:
                del_arc_l = float(del_arc_l_nm/self.filaments_arr[i].max_arc_l)
                roll_size = 0.1
                print(self.filaments_arr[i].max_arc_l,del_arc_l,roll_size)
            #rolls over the arc len in a percent
            roll_win_size = roll_size*del_arc_l
            temp_win_a = start_AL
            temp_win_b = temp_win_a+del_arc_l
            win_num = 0
            while temp_win_a+del_arc_l<=end_AL:
                #fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5));plt.cla()
                #plt.figure(i,figsize=(7,9));plt.clf()
                rgn1 = [temp_win_a,temp_win_b]
                rgn_arr = [rgn1];i_rgn = 0
                temp_rgn = rgn_arr[i_rgn]

                self.update_arc_len_win(temp_rgn)
                
                i_fil_mask =  self.masked_filament  
                L_p_green = i_fil_mask.L_p;R_sq_gr = i_fil_mask.L_p_R_sq
                E_B = i_fil_mask.E_bend_fil;pt_curv_0mean = i_fil_mask.pt_curv - np.mean(i_fil_mask.pt_curv )
                
                arc_l_seg = i_fil_mask.arc_l_seg;end_end_seg = i_fil_mask.end_to_end_dist_2_seg
                
                
                fitted_end_end_dist_sq = i_fil_mask.wlc_LP_fit(arc_l_seg,L_p_green)
                self.save_rgn(f"roll_win_auto_{100*roll_size:.2f}_percent.csv",False)

                temp_win_a +=roll_win_size
                temp_win_b = temp_win_a+del_arc_l;win_num+=1
                
                # plt.title(f"Frame_number_{i}_win_number_{win_num}")
                # plt.scatter(arc_l_seg,end_end_seg,c = color_arr[i_rgn],marker = marker_arr[i_rgn],s = 250)
                # plt.plot(arc_l_seg,fitted_end_end_dist_sq,color = color_arr[i_rgn],linestyle ='dashed',label = f"Segemnt Dimer_{i_rgn+1} L_p {L_p_green:.1f}",linewidth = 3.5,alpha =0.8)
                
               
                '''
                
                ax2.hist(pt_curv_0mean,density =True,color = color_arr[i_rgn],histtype = "step")
                xmin, xmax = plt.xlim()
                x = np.linspace(xmin, xmax, 100)
                mean_fit,std_fit=norm.fit(pt_curv_0mean)
                
                y = norm.pdf(x, mean_fit, std_fit)
                ax2.plot(x, y,'--',color = color_arr[i_rgn],linewidth = 4,label = f"L_p_dist{i_fil_mask.L_p_dist} ")
                '''
                


                #fig.suptitle(f"Frame_number_{i}_win_number_{win_num}")
                #ax1.set_title("WLC");ax1.legend()
                #ax2.set_title("LP from dist");ax2.legend()
                #ax1.set_xlabel("Arc length nm");ax1.set_ylabel("End to end distance square nm^2")
                #ax2.set_xlabel("mean corrected kappa")
                #fig.show()
                #plt.legend(frameon= False);plt.show()
                
                #plt.savefig(self.save_path+f"Frame_number_{i}_win_number_{win_num}.png");break
     
    def curv_viz_disp_plot(self,win_i):
        '''
        2D histogram of the curvature at a particular segment over all segments 
        
        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/win_i)), False)
        fig, ax = plt.subplots(nrows=1, ncols=len(lin_arr),sharey=True)
        df_test = pd.DataFrame( columns=['avg_arc_L','pt_curvatures','frame_number'])
        max_len = 0
        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i

            if la <0.5:
                fil_label = "fil1"
                avg_arc_l = 0.5*(temp_win_a+temp_win_b)

            else:
                fil_label = "fil2"
                avg_arc_l = 1- 0.5*(temp_win_a+temp_win_b)

            
            for j in range(self.num_frames):
                self.i_frame = j
                self.update_arc_len_win([temp_win_a,temp_win_b])
                i_fil_mask =  self.masked_filament
                if max_len<len(i_fil_mask.pt_curv):
                    max_len = len(i_fil_mask.pt_curv)
                for k in range(len(i_fil_mask.pt_curv)):
                    N_test = len(df_test)
                    df_test.loc[N_test,'file_label']  = self.file_label
                    df_test.loc[N_test,'avg_arc_L']  = np.round(avg_arc_l,3)
                    df_test.loc[N_test,'pt_curvatures']  =i_fil_mask.pt_curv[k]
                    df_test.loc[N_test,'frame_number']  = j
                    df_test.loc[N_test,'fil_label']  = fil_label
                
        ax = sns.displot(df_test,
                    x = 'frame_number',y='pt_curvatures', binwidth=(1, .005),
                    cbar=False,row = 'fil_label',col= 'avg_arc_L',hue = 'fil_label')
        
        ax.fig.suptitle(self.file_label)
        ax.fig.savefig(self.save_path+f'/{self.file_label}_curv_dist_overtime.png') 
        ax.fig.savefig(f"/Users/yogehs/Downloads/SbcC/{self.file_label}_curv_dist_overtime.png") 

           
        return df_test,max_len
    def tangential_angle_L_p(self,win_i):
        '''
        dynamics of mid point of each segment, using the X_c and Y_c

        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/win_i)), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0

        for i in range(len(lin_arr)):
            la = lin_arr[i]

            if la <0.5:
                temp_win_a =0;temp_win_b= la+win_i

                fil_label = "fil1"
                avg_arc_l = 0.5*(temp_win_a+temp_win_b)
                total_arl = temp_win_b

            else:
                temp_win_a =la;temp_win_b= 1
                fil_label = "fil2"
                avg_arc_l = 1- 0.5*(temp_win_a+temp_win_b)
                total_arl = 1-temp_win_a

            
            for j in range(self.num_frames):
                self.i_frame = j
                self.update_arc_len_win([temp_win_a,temp_win_b])
                
                i_filament =  self.filaments_arr[self.i_frame]
                max_len = i_filament.max_arc_l
                i_fil_mask =  self.masked_filament


                tan_angle_mask =  i_filament.tan_angle_arr[self.mask_win_gui]
                avg_tan_angle = np.average(i_filament.pt_curv_ys[self.mask_win_gui])*total_arl*max_len
                #avg_tan_angle = np.average(tan_angle_mask)

                
                if fil_label =='fil1':
                    tan_ang_0 = tan_angle_mask[1]
                    tan_ang_1 = tan_angle_mask[-1]
                    tan_ang_diff = tan_ang_1-tan_ang_0
                    delta_tan_ang =tan_ang_diff - avg_tan_angle
                    print("fil1 ",delta_tan_ang,avg_tan_angle,la,j)
                else:
                    tan_ang_0 = tan_angle_mask[-1]
                    tan_ang_1 = tan_angle_mask[1]
                    tan_ang_diff = tan_ang_1-tan_ang_0
                    
                    
                    #tan_ang_0 = tan_angle_mask[1]
                    #tan_ang_1 = tan_angle_mask[-1]
                    #tan_ang_diff = tan_ang_1-tan_ang_0
                    delta_tan_ang = avg_tan_angle-tan_ang_diff
                    print("fil1 ",delta_tan_ang,avg_tan_angle)


                N_test = len(df_test)
                df_test.loc[N_test,'file_label']  = self.file_label
                df_test.loc[N_test,'config_label']  = self.file_label[0]

                df_test.loc[N_test,'avg_arc_L']  = np.round(avg_arc_l,3)
                df_test.loc[N_test,'total_arc_L']  = total_arl

                df_test.loc[N_test,'avg_arc_L_whole_seg']  = np.round(0.5*(temp_win_a+temp_win_b),3)
                df_test.loc[N_test,'tan_angle_rad_1_0']  =tan_ang_diff
                df_test.loc[N_test,'tan_angle_deg_1_0']  = np.degrees(tan_ang_diff)

                df_test.loc[N_test,'cosof_tan_angle_deg_1_0']  =np.cos(tan_ang_diff)
                df_test.loc[N_test,'log_cosof_tan_angle_deg_1_0']  = np.log(np.cos(tan_ang_diff))


                df_test.loc[N_test,'cosof_delta_ang']  =np.cos(delta_tan_ang)
                df_test.loc[N_test,'log_cosof_delta_ang']  = np.log(np.cos(delta_tan_ang))





                df_test.loc[N_test,'end_end_dist_2']  =i_fil_mask.end_to_end_dist_2_raw[-1]
                df_test.loc[N_test,'end_end_dist']  =np.sqrt(i_fil_mask.end_to_end_dist_2_raw[-1])
                df_test.loc[N_test,'bend_angle_deg']  =i_fil_mask.bending_ang_fil
                df_test.loc[N_test,'L_p']  =i_fil_mask.L_p
                df_test.loc[N_test,'L_p_Rsq']  =i_fil_mask.L_p_R_sq

                df_test.loc[N_test,'contour_L']  =i_fil_mask.max_arc_l
                df_test.loc[N_test,'contour_s']  =i_fil_mask.max_arc_l *total_arl

                #df_test.loc[N_test,f'{fil_label}_end_end_dist_{avg_arc_l}']  =i_fil_mask.end_to_end_dist_2_raw[-1]

                df_test.loc[N_test,'frame_number']  = j
                df_test.loc[N_test,'fil_label']  = fil_label



        df_test.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_cosine_corr_L_p.csv")
        #df_test.to_csv(f"{self.save_path}{self.file_label}_lateral_fluct.csv")

        return(df_test)
    def lateral_fluct_dist(self,win_i):
        '''
        dynamics of mid point of each segment, using the X_c and Y_c

        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/win_i)), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0

        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i

            if la <0.5:
                fil_label = "fil1"
                avg_arc_l = 0.5*(temp_win_a+temp_win_b)

            else:
                fil_label = "fil2"
                avg_arc_l = 1- 0.5*(temp_win_a+temp_win_b)
                
            temp_list = []
            
            for j in range(self.num_frames):
                self.i_frame = j
                
                self.update_arc_len_win([temp_win_a,temp_win_b])
                
                i_filament =  self.filaments_arr[self.i_frame]
             

                x_c_win = i_filament.X_c[self.mask_win_gui] ;y_c_win = i_filament.Y_c[self.mask_win_gui]
                x_win = i_filament.X[self.mask_win_gui] ;y_win = i_filament.Y[self.mask_win_gui]
                i_fil_mask =  self.masked_filament
                
                temp_list += [i_fil_mask.end_to_end_dist_2_raw[-1]]
                if j==0:
                    disp_a_b = 0 
                else:
                    temp_vec_b_i  = np.array(findMiddle_arr(x_c_win),findMiddle_arr(y_c_win)) 
                    disp_a_b  = np.linalg.norm(temp_a_i_1-temp_vec_b_i)
 
                N_test = len(df_test)
                df_test.loc[N_test,'file_label']  = self.file_label
                df_test.loc[N_test,'avg_arc_L']  = np.round(avg_arc_l,3)
                df_test.loc[N_test,'avg_arc_L_whole_seg']  = np.round(0.5*(temp_win_a+temp_win_b),3)

                df_test.loc[N_test,'X_C']  =  findMiddle_arr(x_c_win)
                df_test.loc[N_test,'Y_C']  = findMiddle_arr(y_c_win)
                
                df_test.loc[N_test,'X_raw']  = findMiddle_arr(x_win)
                df_test.loc[N_test,'Y_raw']  = findMiddle_arr(y_win)


                df_test.loc[N_test,'X_0']  = i_filament.X[0]
                df_test.loc[N_test,'Y_0']  = i_filament.Y[0]
                # this is the distance between two consecutive frames
                
                df_test.loc[N_test,'cons_frame_disp']  = disp_a_b

                #this is the angle between the end points of the protein image in a frame 
                df_test.loc[N_test,'thetha_deg']  = np.degrees(i_filament.theta_orient)





                df_test.loc[N_test,'end_end_dist_2']  =i_fil_mask.end_to_end_dist_2_raw[-1]
                df_test.loc[N_test,'end_end_dist']  =np.sqrt(i_fil_mask.end_to_end_dist_2_raw[-1])
                df_test.loc[N_test,'bend_angle_deg']  =i_fil_mask.bending_ang_fil
                df_test.loc[N_test,'L_p']  =i_fil_mask.L_p
                df_test.loc[N_test,'L_p_Rsq']  =i_fil_mask.L_p_R_sq

                df_test.loc[N_test,'l_seg']  =i_fil_mask.max_arc_l

                #df_test.loc[N_test,f'{fil_label}_end_end_dist_{avg_arc_l}']  =i_fil_mask.end_to_end_dist_2_raw[-1]

                df_test.loc[N_test,'frame_number']  = j
                df_test.loc[N_test,'fil_label']  = fil_label
                temp_a_i_1 = np.array(findMiddle_arr(x_c_win),findMiddle_arr(y_c_win)) 



        df_test.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_lateral_fluct.csv")
        #df_test.to_csv(f"{self.save_path}{self.file_label}_lateral_fluct.csv")

        return(df_test)


    def radial_dist_func(self,win_i):
        '''
        dynamics of end end distr segment wise
        #TODO find MSD for different time intervals
        #TODO what is the avaerage angle 
        
        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/win_i)), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0

        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i

            if la <0.5:
                fil_label = "fil1"
                avg_arc_l = 0.5*(temp_win_a+temp_win_b)

            else:
                fil_label = "fil2"
                avg_arc_l = 1- 0.5*(temp_win_a+temp_win_b)
            temp_list = []
            
            for j in range(self.num_frames):
                self.i_frame = j
                self.update_arc_len_win([temp_win_a,temp_win_b])
                i_fil_mask =  self.masked_filament
                temp_list += [i_fil_mask.end_to_end_dist_2_raw[-1]]

                N_test = len(df_test)
                df_test.loc[N_test,'file_label']  = self.file_label
                df_test.loc[N_test,'avg_arc_L']  = np.round(avg_arc_l,3)
                df_test.loc[N_test,'avg_arc_L_whole_seg']  = np.round(0.5*(temp_win_a+temp_win_b),3)

                df_test.loc[N_test,'end_end_dist_2']  =i_fil_mask.end_to_end_dist_2_raw[-1]
                df_test.loc[N_test,'end_end_dist']  =np.sqrt(i_fil_mask.end_to_end_dist_2_raw[-1])
                df_test.loc[N_test,'tangential_angle_deg']  =i_fil_mask.bending_ang_fil
                df_test.loc[N_test,'L_p']  =i_fil_mask.L_p
                print(avg_arc_l,i_fil_mask.L_p,i_fil_mask.L_p_R_sq)
                df_test.loc[N_test,'L_p_Rsq']  =i_fil_mask.L_p_R_sq

                df_test.loc[N_test,'l_seg']  =i_fil_mask.max_arc_l

                #df_test.loc[N_test,f'{fil_label}_end_end_dist_{avg_arc_l}']  =i_fil_mask.end_to_end_dist_2_raw[-1]

                df_test.loc[N_test,'frame_number']  = j
                df_test.loc[N_test,'fil_label']  = fil_label

        df_msd,df_msd_arr =find_MSD_metric_ori(df_test)
        mask = df_msd['deltaT']<=(df_msd['deltaT'].max()+1)*0.5
        df_msd = df_msd[mask] 
        '''

        fg =sns.FacetGrid(df_msd, col="avg_arc_L",row = "fil_label",hue = 'fil_label',height = 6) 
        fg.map_dataframe(errplot, "deltaT", "MSD", "std_MSD")
        fg.set(xscale = 'log');fg.set(yscale= 'log')

        fg.fig.suptitle(df_msd['file_label'][0]+"MSD_radial")
        fg.add_legend()

        #fg.set(ylim=(-0.2,1),xlim=(0,20))
        fg.fig.savefig(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/radial/{self.file_label}_MSD_rad_overtime.png")
        fg =sns.FacetGrid(df_msd, col="avg_arc_L",hue = "fil_label",row = "fil_label",height=6)  
        fg.map_dataframe(errplot, "deltaT", "ang_MSD", "std_ang_MSD")
        fg.set(xscale = 'log');fg.set(yscale= 'log')

        fg.fig.suptitle(df_msd['file_label'][0]+"MSD_ang")
        #fg.set(yscale='log',xscale='log')

        #fg.set(ylim=(-0.2,2))
        #fg.set(xlim=(0,20))

        fg.add_legend()
        fg.fig.savefig(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/angle/{self.file_label}_MSD_angle_overtime.png")
        '''
        df_msd.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_MSD.csv")
        df_msd.to_csv(f"{self.save_path}{self.file_label}_MSD.csv")
        df_test.to_csv(f"{self.save_path}{self.file_label}_raw_preMSD.csv")

        return(df_msd,df_msd_arr)
    def R_TAN_BEND_MSD(self,win_i):
        '''
        dynamics of end end distr segment wise
        #TODO find MSD for different time intervals
        #TODO what is the avaerage angle 
        
        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/(win_i))), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0    
        fps = long_vids_fps[self.file_label]


        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i
            
            temp_win_a_2 =la+win_i;temp_win_b_2= la+2*win_i

            if la <0.5:
                fil_label = "fil1"
                avg_arc_l = np.round(0.5*(temp_win_a+temp_win_b),3)

            else:
                fil_label = "fil2"
                avg_arc_l = np.round(1- 0.5*(temp_win_a+temp_win_b),3)
            temp_list = []
            
            for j in range(self.num_frames):
                
                if temp_win_b_2<=1:
                    if fil_label =='fil1':
                        self.i_frame = j
                        print(temp_win_a_2,temp_win_b_2)
                        self.update_arc_len_win([temp_win_a_2,temp_win_b_2])
                        i_fil_mask =  self.masked_filament
                        pt_C = np.array([i_fil_mask.X_c[-1],i_fil_mask.Y_c[-1]])
        
                        self.update_arc_len_win([temp_win_a,temp_win_b])
                        i_fil_mask =  self.masked_filament
                        
                        #finding the bend angle 
                        pt_A = np.array([i_fil_mask.X_c[0],i_fil_mask.Y_c[0]])
                        pt_B = np.array([i_fil_mask.X_c[-1],i_fil_mask.Y_c[-1]])
                        
                        bend_angle =find_angle(pt_A, pt_B, pt_C)
                        
                        
                    else:
                        
                        #as the sense of filament reading is changing,
                        #the def of pta and ptb and ptc needs to change 
                        self.i_frame = j

                        self.update_arc_len_win([temp_win_a_2,temp_win_b_2])
                        i_fil_mask =  self.masked_filament
                        pt_A = np.array([i_fil_mask.X_c[-1],i_fil_mask.Y_c[-1]])
                        pt_B = np.array([i_fil_mask.X_c[0],i_fil_mask.Y_c[0]])
    
            
                        self.update_arc_len_win([temp_win_a,temp_win_b])
                        i_fil_mask =  self.masked_filament
                        pt_C = np.array([i_fil_mask.X_c[0],i_fil_mask.Y_c[0]])
                        
                        bend_angle =find_angle(pt_A, pt_B, pt_C)

                else:
                    self.i_frame = j
                    self.update_arc_len_win([temp_win_a,temp_win_b])

                    i_fil_mask =  self.masked_filament

                    bend_angle = 0

                i_filament =  self.filaments_arr[self.i_frame]

                if fil_label =='fil1':
                    tan_angle_mask =  i_filament.tan_angle_arr[self.mask_win_gui]
                    tan_angle = tan_angle_mask[1]
                    pt_A = np.array([i_fil_mask.X_c[0],i_fil_mask.Y_c[0]])
                    pt_B = np.array([i_fil_mask.X_c[-1],i_fil_mask.Y_c[-1]])
                    
                    pt_mid = np.array([findMiddle_arr(i_fil_mask.X_c),findMiddle_arr(i_fil_mask.Y_c)])
                    
                    angle_wrt_r = find_angle(pt_A, pt_mid,pt_B )
                else:
                    tan_angle_mask =  i_filament.tan_angle_arr[self.mask_win_gui]
                    tan_angle = tan_angle_mask[-1]

                    pt_A = np.array([i_fil_mask.X_c[-1],i_fil_mask.Y_c[-1]])
                    pt_B = np.array([i_fil_mask.X_c[0],i_fil_mask.Y_c[0]])
                    
                    pt_mid = np.array([findMiddle_arr(i_fil_mask.X_c),findMiddle_arr(i_fil_mask.Y_c)])
                    
                    angle_wrt_r = find_angle(pt_A, pt_mid,pt_B )

                N_test = len(df_test)
                df_test.loc[N_test,'file_label']  = self.file_label
                df_test.loc[N_test,'config_label']  = self.file_label[0]
                df_test.loc[N_test,'folder_number']  = self.file_label[1:]

                df_test.loc[N_test,'frame_number']  = j
                df_test.loc[N_test,'time_s']  = j/fps


                df_test.loc[N_test,'avg_arc_L']  = avg_arc_l
                df_test.loc[N_test,'avg_arc_L_whole_seg']  = np.round(0.5*(temp_win_a+temp_win_b),3)

                df_test.loc[N_test,'end_end_dist_2']  =i_fil_mask.end_to_end_dist_2_raw[-1]
                df_test.loc[N_test,'end_end_dist']  =np.sqrt(i_fil_mask.end_to_end_dist_2_raw[-1])
                df_test.loc[N_test,'tangential_angle_rad']  = tan_angle
                df_test.loc[N_test,'bend_angle_deg']  = bend_angle
                df_test.loc[N_test,'angle_wrt_r']  = angle_wrt_r

                
                df_test.loc[N_test,'L_p']  =i_fil_mask.L_p
                df_test.loc[N_test,'L_p_Rsq']  =i_fil_mask.L_p_R_sq

                df_test.loc[N_test,'l_seg']  =i_fil_mask.arc_l[-1]-i_fil_mask.arc_l[0]

                #df_test.loc[N_test,f'{fil_label}_end_end_dist_{avg_arc_l}']  =i_fil_mask.end_to_end_dist_2_raw[-1]

                df_test.loc[N_test,'fil_label']  = fil_label

        df_msd,df_msd_arr =find_MSD_metric(df_test)
        
        #masking the first half of DT 
        mask = df_msd['deltaT']<=(df_msd['deltaT'].max()+1)*0.5
        df_msd = df_msd[mask] 
        
        mask_fil1 = df_msd['fil_label']=='fil1'
        df_msd['avg_arc_L_bend_angle'] = 0
        df_msd['avg_arc_L_bend_angle'][mask_fil1] = np.round(df_msd['avg_arc_L'][mask_fil1],3)
        df_msd['avg_arc_L_bend_angle'][~mask_fil1 ]= np.round(df_msd['avg_arc_L'][~mask_fil1] - 0.1,3)
        
        
        '''

        fg =sns.FacetGrid(df_msd, col="avg_arc_L",row = "fil_label",hue = 'fil_label',height = 6) 
        fg.map_dataframe(errplot, "deltaT", "MSD", "std_MSD")
        fg.set(xscale = 'log');fg.set(yscale= 'log')

        fg.fig.suptitle(df_msd['file_label'][0]+"MSD_radial")
        fg.add_legend()

        #fg.set(ylim=(-0.2,1),xlim=(0,20))
        fg.fig.savefig(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/radial/{self.file_label}_MSD_rad_overtime.png")
        fg =sns.FacetGrid(df_msd, col="avg_arc_L",hue = "fil_label",row = "fil_label",height=6)  
        fg.map_dataframe(errplot, "deltaT", "ang_MSD", "std_ang_MSD")
        fg.set(xscale = 'log');fg.set(yscale= 'log')

        fg.fig.suptitle(df_msd['file_label'][0]+"MSD_ang")
        #fg.set(yscale='log',xscale='log')

        #fg.set(ylim=(-0.2,2))
        #fg.set(xlim=(0,20))

        fg.add_legend()
        fg.fig.savefig(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/angle/{self.file_label}_MSD_angle_overtime.png")
        df_msd.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_MSD.csv")
        df_msd.to_csv(f"{self.save_path}{self.file_label}_MSD.csv")
        df_test.to_csv(f"{self.save_path}{self.file_label}_raw_preMSD.csv")

        return(df_msd,df_msd_arr)    
        '''

        #df_msd.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_{win_i}_MSD_R_TAN_BEND.csv")
        #df_msd.to_csv(f"{self.save_path}{self.file_label}_{win_i}_MSD_R_TAN_BEND.csv")
        #df_msd_arr.to_csv(f"{self.save_path}{self.file_label}_{win_i}_MSD_distributions.csv")
        df_msd_arr.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_{win_i}_MSD_distributions.csv")
        #df_test.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_{win_i}_pre_MSD.csv")

        return(df_msd,df_msd_arr,df_test)
    def avg_curv_viz(self,win_i):
        '''

        
        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/(win_i))), False)
        df_test = pd.DataFrame( columns=['avg_arc_L','end_end_dist','frame_number'])
        max_len = 0

        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i
            avg_arc_l = np.round(0.5*(temp_win_a+temp_win_b),3)

            temp_list=[];temp_curv_all_frames=[]
            for j in range(self.num_frames):

                self.i_frame = j
                self.update_arc_len_win([temp_win_a,temp_win_b])
                i_fil_mask =  self.masked_filament

                temp_curv_all_frames+=i_fil_mask.pt_curv.tolist()
                
            
            N_test = len(df_test)
            df_test.loc[N_test,'file_label']  = self.file_label
            df_test.loc[N_test,'config_label']  = self.file_label[0]
            df_test.loc[N_test,'folder_number']  = self.file_label[1:]



            df_test.loc[N_test,'avg_arc_L']  = avg_arc_l

            

            df_test.loc[N_test,'l_seg']  =i_fil_mask.max_arc_l


            df_test.loc[N_test,'avg_pt_curv_ys']  = np.mean(temp_curv_all_frames)
            df_test.loc[N_test,'std_pt_curv_ys']  = np.std(temp_curv_all_frames)

        #df_test.to_csv(f"/Users/yogehs/Downloads/SbcC/MSD_scatters/{self.file_label}_{win_i}_pre_MSD.csv")

        return(df_test)       
    def curv_viz_joint_plot(self,win_i):
        '''
        2D histogram of the curvature at a particular segment over all segments 
        
        '''
        lin_arr = np.linspace(0,1,int(np.ceil(1/win_i)), False)
        fig, ax = plt.subplots(nrows=1, ncols=len(lin_arr),sharey=True)
        df_test = pd.DataFrame()

        for i in range(len(lin_arr)):
            la = lin_arr[i]
            temp_win_a =la;temp_win_b= la+win_i
            if la <0.5:
                fil_label = "1"
                avg_arc_l = 0.5*(temp_win_a+temp_win_b)

            else:
                fil_label = "2"
                avg_arc_l = 1- 0.5*(temp_win_a+temp_win_b)
            avg_arc_l=np.round(avg_arc_l,3)
            temp_list=[];temp_frame=[]
            for j in range(self.num_frames):

                self.i_frame = j
                self.update_arc_len_win([temp_win_a,temp_win_b])
                i_fil_mask =  self.masked_filament
                temp_list += i_fil_mask.pt_curv.tolist()
                temp_frame += [j]*len(i_fil_mask.pt_curv.tolist())
            df_temp = pd.DataFrame({f'{fil_label}_pt_curvature_{avg_arc_l}':temp_list,f'{fil_label}_temp_frame_{avg_arc_l}':temp_frame})

            df_test = pd.concat([df_test,df_temp], ignore_index=False, axis=1)

        return df_test
 