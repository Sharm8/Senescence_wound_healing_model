from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
import random
import math
import Parameters_WH as param
import csv
import glob
import os
import pandas as pd
import datetime
import time
import shutil

# 1 pixel = 3 microns
# 1 mcs = 27 secs 
    
 
class MechanismsSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)  


    def start(self):
        return 
                      

    def step(self,mcs):
       
        field_Proteinase = self.field.Proteinase
        field_Inflammation_signal = self.field.Inflammation_signal
        field_PDGF = self.field.PDGF
        field_CSF = self.field.CSF
        lx = self.dim.x 
        ly = self.dim.y 
                
        if mcs % 50 == 0:
            
            # Fibroblast and Myofibroblast proliferation in response to PDGF
            for cell in self.cell_list_by_type(self.FIBROBLAST, self.MYOFIBROBLAST):
                PDGF_conc = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                # activated fibroblasts only
                if PDGF_conc > param.actF_col3:   
                    neighbor_list = self.get_cell_neighbor_data_list(cell)
                    surface_area_with_medium = neighbor_list.common_surface_area_with_cell_types(cell_type_list=[self.MEDIUM])
                    surface = cell.surface
                    # contact inhibition 
                    rs = surface_area_with_medium/surface 
                    threshold = 0.05
                    if rs < threshold:
                        cell.targetVolume = cell.volume
                        cell.lambdaVolume = 2.0
                        
                    elif rs >= threshold: 
                        # cell growth
                        Gmax_Fib = ((param.Gmax_Fib * param.Fib_TarVol * param.mcs_to_sec) / (24*60*60*1))  
                        PDGF_conc_halfmaxgrowth = param.PDGF_conc_halfmaxgrowth 
                        cell.targetVolume += (Gmax_Fib * PDGF_conc / (PDGF_conc_halfmaxgrowth + PDGF_conc)) * 50
                        cell.targetSurface = math.sqrt(cell.targetVolume)*4
                        PDGF_conc = max(0, (PDGF_conc - param.PDGF_conc_halfmaxgrowth))
                           
                            
            
            # Macrophage proliferation in response to CSF 
            for cell in self.cell_list_by_type(self.MACROPHAGE):
                neighbor_list = self.get_cell_neighbor_data_list(cell)
                surface_area_with_medium = neighbor_list.common_surface_area_with_cell_types(cell_type_list=[self.MEDIUM])
                surface = cell.surface
                # contact inhibition
                rs = surface_area_with_medium/surface 
                threshold = 0.05
                CSF_conc = field_CSF[cell.xCOM, cell.yCOM, cell.zCOM]
                if rs < threshold:
                    cell.targetVolume = cell.volume
                    cell.lambdaVolume = 2.0
                    
                elif rs >= threshold:
                   # cell growth  
                    Gmax_Mac = ((param.Gmax_Mac * param.Mac_TarVol * param.mcs_to_sec) / (24*60*60*1)) 
                    CSF_conc_halfmaxgrowth = param.CSF_conc_halfmaxgrowth   
                    cell.targetVolume += (Gmax_Mac * CSF_conc / (CSF_conc_halfmaxgrowth + CSF_conc)) * 50
                    cell.targetSurface = math.sqrt(cell.targetVolume)*4
                    CSF_conc = max(0, (CSF_conc - param.CSF_conc_halfmaxgrowth))     
                
         
           
        
        if mcs >= 5:
            
            if mcs % 50 == 0:
                
                # ECM production by activated fibroblasts
                for cell in self.cell_list_by_type(self.FIBROBLAST): 
                    
                    PDGF_value = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                    
                    if PDGF_value > param.actF_col3:
                        PDGF_value = max(0, (PDGF_value - param.actF_col3))
                        # new ecm
                        cell_size = int(param.Ecm_TarSur/4)
                        x = int(random.randrange(0,200-cell_size,1))
                        y = int(random.randrange(0,200-cell_size,1))
                        cell = self.cell_field[x, y, 0]
                        if not cell:
                            ecm = self.new_cell(self.ECM)
                            self.cell_field[x:x + cell_size - 1, y:y + cell_size - 1, 0] = ecm
                            ecm.targetVolume = param.Ecm_TarVol 
                            ecm.lambdaVolume = 2.0
                            ecm.targetSurface = param.Ecm_TarSur 
                            ecm.lambdaSurface = 2.0
                            
                    
            
            
            if mcs % 200 == 0:
                
                
                               
                # Activated fibroblasts differentiate into myofibroblasts
                fib_len = len(self.cell_list_by_type(self.FIBROBLAST))
                myofibroblast_percentage = ((np.random.binomial(n = fib_len, p = param.Fib_actF_myoF_chance, size = 1))/fib_len)*100
                for cell in self.cell_list_by_type(self.FIBROBLAST):
                    
                    if cell.dict['deactivated_myoF'] != 1:
                        
                        if 50<cell.xCOM<150:
                            if 50<cell.yCOM<200:
                                chance = random.randint(0, 100)
                                if chance <= myofibroblast_percentage:
                                    PDGF_value = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                                    neighbour_list = self.get_cell_neighbor_data_list(cell) 
                                    common_area_with_ECM = neighbour_list.common_surface_area_with_cell_types(cell_type_list=[self.ECM])
                                    # mechanical stimulation through contact with ECM and PDGF exposure
                                    if (common_area_with_ECM > param.Col_myoF_deact and PDGF_value > param.actF_col3):
                                        cell.type = self.MYOFIBROBLAST
                                        cell.targetVolume = param.Myofib_TarVol 
                                        cell.lambdaVolume = 2.0
                                        cell.targetSurface = param.Myofib_TarSur 
                                        cell.lambdaSurface = 2.0
                                        cell.dict['myo_start'] = mcs
                                        PDGF_value = max(0, (PDGF_value - param.actF_col3))
                            
                # ECM production by myofibroblasts
                for cell in self.cell_list_by_type(self.MYOFIBROBLAST):   
                    PDGF_value = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                    if PDGF_value > param.actF_col3:
                        PDGF_value = max(0, (PDGF_value - param.actF_col3))
                        
                        n = 2
                        for i in range(n):
                            cell_size = int(param.Ecm_TarSur/4)
                            x = int(random.randrange(0,200-cell_size,1))
                            y = int(random.randrange(0,200-cell_size,1))
                            #x = min(int(200-cell_size), int(random.randrange(int(cell.xCOM), int(cell.xCOM)+40, 1)))     
                            #y = min(int(200-cell_size), int(random.randrange(int(cell.yCOM), int(cell.yCOM)+40, 1)))
                            cell = self.cell_field[x, y, 0]
                            if not cell:
                                ecm = self.new_cell(self.ECM)
                                self.cell_field[x:x + cell_size - 1, y:y + cell_size - 1, 0] = ecm
                                ecm.targetVolume = param.Ecm_TarVol 
                                ecm.lambdaVolume = 2.0
                                ecm.targetSurface = param.Ecm_TarSur 
                                ecm.lambdaSurface = 2.0
                        
                     
                # Myofibroblast death, de-differentiation or senescence 
                myofib_len = len(self.cell_list_by_type(self.MYOFIBROBLAST))              
                senescence_percentage_ccn1 = ((np.random.binomial(n = myofib_len, p = param.senescence_percentage_param, size = 1))/myofib_len)*100
                for cell in self.cell_list_by_type(self.MYOFIBROBLAST):
                    if mcs > cell.dict['myo_start'] + ((24*60*60*3)/param.mcs_to_sec):
                        neighbour_list = self.get_cell_neighbor_data_list(cell) 
                        common_area_with_ECM = neighbour_list.common_surface_area_with_cell_types(cell_type_list=[self.ECM])
                        chance = random.randint(0, 100)
                        if chance >= senescence_percentage_ccn1:
                            # exposure to inflammation or release of mechanical tension
                            if common_area_with_ECM < (param.Col_myoF_deact) or field_Inflammation_signal[cell.xCOM, cell.yCOM, cell.zCOM] > param.IL6:
                                chance = random.randint(0, 100)
                                # de-differentiation to fobriblast state
                                if chance <= 50:
                                   cell.type = self.FIBROBLAST
                                   cell.dict['deactivated_myoF'] = 1
                                   
                                else:
                                   # myofibroblast death
                                   self.delete_cell(cell)
                                    
                        else:
                            # During the proliferative phase of wound healing
                            if mcs < 40000:                                                                                                                                            
                                neighbour_list = self.get_cell_neighbor_data_list(cell) 
                                common_area_with_ECM = neighbour_list.common_surface_area_with_cell_types(cell_type_list=[self.ECM])
                                # CCN1-mediated myofibroblast differentiation
                                if common_area_with_ECM > param.Col_myoF_deact: 
                                    cell.type = self.SNCMYOF
                                    cell.targetVolume = param.Sncmyof_TarVol 
                                    cell.lambdaVolume = 2.0 
                                    cell.targetSurface = param.Sncmyof_TarSur 
                                    cell.lambdaSurface = 2.0
                                    cell.dict['snc_start'] = mcs  
                                    cell.dict['snc_mech'] = 'CCN1'
                
                # Paracrine and juxtacrine secondary senescence in myofibroblasts
                fib_len = len(self.cell_list_by_type(self.FIBROBLAST))                
                senescence_percentage = ((np.random.binomial(n = fib_len, p = param.senescence_percentage_param, size = 1))/fib_len)*100
                for cell in self.cell_list_by_type(self.MYOFIBROBLAST):
                    if mcs > cell.dict['myo_start'] + ((24*60*60*3)/param.mcs_to_sec):
                        chance = random.randint(0, 100)
                        if chance <= senescence_percentage:
                            PDGF_value = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                            Inflammation_value = field_Inflammation_signal[cell.xCOM, cell.yCOM, cell.zCOM]
                            for neighbor_tuple in self.get_cell_neighbor_data_list(cell):
                                if neighbor_tuple[0] and neighbor_tuple[0].type == self.SNCMYOF:
                                    # Juxtacrine senescence during the first 3 days
                                    if mcs < neighbor_tuple[0].dict['snc_start'] + ((24*60*60*3)/param.mcs_to_sec):
                                        if neighbor_tuple[1] > 0 and PDGF_value > param.myoF_to_sncF:         
                                            cell.type = self.SNCMYOF
                                            cell.targetVolume = param.Sncmyof_TarVol 
                                            cell.lambdaVolume = 2.0 
                                            cell.targetSurface = param.Sncmyof_TarSur 
                                            cell.lambdaSurface = 2.0
                                            PDGF_value = max(0, (PDGF_value - param.myoF_to_sncF))
                                            cell.dict['snc_start'] = mcs
                                            cell.dict['snc_mech'] = 'JUX'    
                                            
                                    # paracrine senescence during the following days
                                    elif mcs >= neighbor_tuple[0].dict['snc_start'] + ((24*60*60*3)/param.mcs_to_sec):
                                        if Inflammation_value > 0:
                                            cell.type = self.SNCMYOF
                                            cell.targetVolume = param.Sncmyof_TarVol 
                                            cell.lambdaVolume = 2.0 
                                            cell.targetSurface = param.Sncmyof_TarSur 
                                            cell.lambdaSurface = 2.0
                                            Inflammation_value = max(0, (Inflammation_value - param.IL6))
                                            cell.dict['snc_start'] = mcs                    
                                            cell.dict['snc_mech'] = 'PARA'    
                
                
                # cell speed in simulation ###
                
                #for cell in self.cell_list_by_type(self.FIBROBLAST): # or cell list by type, depending on what you want
                    #if 50<cell.xCOM<150:
                        #if 50<cell.xCOM<150:
                            #vx = cell.xCOM - cell.xCOMPrev 
                            #vy = cell.yCOM - cell.yCOMPrev
                            #vz = cell.zCOM - cell.zCOMPrev
                    #cell.dict['prev_x'] = cell.xCOM
                    
                            #print(vx, vy, vz)
                
    #def finish(self):
        

    def on_stop(self):
        
        return

        
class PlotsSteppable(SteppableBasePy):     
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.cell_plot = None
        self.ECM_plot = None
        self.Field_plot = None
        self.area_plot = None
        

    def start(self):
        

        self.timestamp = '{:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.datetime.now())
        self.startTime = time.time()
        
        self.Fib_plot = self.add_new_plot_window(title='Fibroblast Population Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Cell Count',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        self.Fib_plot.add_plot("FibroblastCount", style='Lines', color='Blue', size=3)
        self.Fib_plot.add_plot("FibroblastTargetCount", style='Dots', color='Red', size=3)
        
        
        
        self.MyoFib_plot = self.add_new_plot_window(title='Myofibroblast Population Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Cell Count',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        self.MyoFib_plot.add_plot("Myofibroblast", style='Lines', color='Blue', size=3)
        self.MyoFib_plot.add_plot("MyofibroblastTargetCount", style='Dots', color='Red', size=3)
        
        
        
        
        self.SncmyoFib_plot = self.add_new_plot_window(title='Senescent Myofibroblast Population Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Cell Count',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        self.SncmyoFib_plot.add_plot("Senescent Myofibroblast", style='Lines', color='Blue', size=3)
        self.SncmyoFib_plot.add_plot("SenescentMyofibroblastTargetCount", style='Dots', color='Red', size=3)
        
        self.SncmyoFib_type_plot = self.add_new_plot_window(title='Senescent Myofibroblast Population Data - Mechanism',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Cell Count',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        self.SncmyoFib_type_plot.add_plot("CCN1", style='Lines', color='Blue', size=3)
        self.SncmyoFib_type_plot.add_plot("Juxtacrine", style='Lines', color='Red', size=3)
        self.SncmyoFib_type_plot.add_plot("Paracrine", style='Lines', color='Green', size=3)
        
        
        self.Mac_plot = self.add_new_plot_window(title='Macrophage Population Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Cell Count',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        self.Mac_plot.add_plot("Macrophage", style='Lines', color='Blue', size=3)
        self.Mac_plot.add_plot("MacrophageTargetCount", style='Dots', color='Red', size=3)
        
        
        
        self.ECM_plot = self.add_new_plot_window(title='ECM Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='ECM',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
                                                     
        self.ECM_plot.add_plot("ECM", style='Lines', color='Blue', size=3)
        self.ECM_plot.add_plot("ECMTarget", style='Dots', color='Red', size=3)
        
        
        self.Field_plot = self.add_new_plot_window(title='Chemical Field Data',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Concentration',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
                                                     
        
        self.Field_plot.add_plot("PDGF", style='Lines', color='Green', size=2)
        self.Field_plot.add_plot("Proteinase", style='Lines', color='Purple', size=2)
        self.Field_plot.add_plot("CSF", style='Lines', color='Yellow', size=2)
        self.Field_plot.add_plot("Inflammation", style='Lines', color='Blue', size=2)
        

        self.Area_plot = self.add_new_plot_window(title='Wound area plot',
                                                     x_axis_title='Time (Monte Carlo Step)',
                                                     y_axis_title='Wound area',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False,
                                                     config_options={'legend': True})
                                                     
        self.Area_plot.add_plot("Area", style='Lines', color='Blue', size=3)
        self.Area_plot.add_plot("AreaTarget", style='Dots', color='Red', size=3)
        
        
        
        
        if self.output_dir is not None:
            Fib_data_file_name = "\\Simulation\\fibroblast.csv"  
            #Fib_data_file_name = "Fib_pop_data.csv"  
            Fib_input_path = self.output_dir + Fib_data_file_name
            Fib_df = pd.read_csv(Fib_input_path, usecols=[0,1], header = None)            
            self.Fib_dict = Fib_df.set_index(0)[1].to_dict()
            for day in self.Fib_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.Fib_plot.add_data_point("FibroblastTargetCount", mcs, self.Fib_dict[day])
               
             
        
            MyoFib_data_file_name = "\\Simulation\\myofibroblast.csv"  
            #MyoFib_data_file_name = "myof_pop_data.csv"  
            MyoFib_input_path = self.output_dir + MyoFib_data_file_name
            MyoFib_df = pd.read_csv(MyoFib_input_path, usecols=[0,1], header = None)            
            self.MyoFib_dict =  MyoFib_df.set_index(0)[1].to_dict()
            for day in self.MyoFib_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.MyoFib_plot.add_data_point("MyofibroblastTargetCount", mcs, self.MyoFib_dict[day])
        
            Mac_data_file_name = "\\Simulation\\macrophage.csv"  
            #Mac_data_file_name = "Mac_pop_data.csv"  
            Mac_input_path = self.output_dir + Mac_data_file_name
            Mac_df = pd.read_csv(Mac_input_path, usecols=[0,1], header = None)            
            self.Mac_dict =  Mac_df.set_index(0)[1].to_dict()
            for day in self.Mac_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.Mac_plot.add_data_point("MacrophageTargetCount", mcs, self.Mac_dict[day])
        
        
            SncmyoFib_data_file_name = "\\Simulation\\senescentmyofibroblast.csv"  
            #SncmyoFib_data_file_name = "sncmyof_pop_data.csv"  
            SncmyoFib_input_path = self.output_dir + SncmyoFib_data_file_name
            SncmyoFib_df = pd.read_csv(SncmyoFib_input_path, usecols=[0,1], header = None)            
            self.SncmyoFib_dict = SncmyoFib_df.set_index(0)[1].to_dict()
            for day in self.SncmyoFib_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.SncmyoFib_plot.add_data_point("SenescentMyofibroblastTargetCount", mcs, self.SncmyoFib_dict[day])
        
        
            ECM_data_file_name = "\\Simulation\\ECM.csv"  
            #ECM_data_file_name = "ECM_data.csv"  
            ECM_input_path = self.output_dir + ECM_data_file_name
            ECM_df = pd.read_csv(ECM_input_path, usecols=[0,1], header = None)            
            self.ECM_dict = ECM_df.set_index(0)[1].to_dict()
            for day in self.ECM_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.ECM_plot.add_data_point("ECMTarget", mcs, self.ECM_dict[day])
            
        
            
        
            Area_data_file_name = "\\Simulation\\wound_closure_rate.csv"  
            #Area_data_file_name = "wound_closure_rate.csv"  
            Area_input_path = self.output_dir + Area_data_file_name
            Area_df = pd.read_csv(Area_input_path, usecols=[0,1], header = None)            
            self.Area_dict = Area_df.set_index(0)[1].to_dict()
            for day in self.Area_dict:
                mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                self.Area_plot.add_data_point("AreaTarget", mcs, self.Area_dict[day])
            
            
        
                    
       
        self.Fib_AbsErr_sum = 0
        self.Fib_AbsErr_list = []
        self.Fib_list = []
        self.Fib_list_wound = []
        
        
        self.MyoFib_AbsErr_sum = 0
        self.MyoFib_AbsErr_list = []
        self.MyoFib_list = []

        self.SncmyoFib_AbsErr_sum = 0
        self.SncmyoFib_AbsErr_list = []
        self.SncmyoFib_list = []
        
        self.Mac_AbsErr_sum = 0
        self.Mac_AbsErr_list = []
        self.Mac_list = []
        
        self.ECM_AbsErr_sum = 0
        self.ECM_AbsErr_list = []
        self.ECM_list = []
        
        self.Area_AbsErr_sum = 0
        self.Area_AbsErr_list = []
        self.Area_list = []
        
        self.mcs_list = []
        self.ccn1_list = []
        self.jux_list = []
        self.para_list = []
        self.PDGF_list = []
        self.CSF_list = []
        self.Inflammation_list = []
        self.Proteinase_list = []
        
    

    def step(self, mcs):
        
        
        if mcs % 100 == 0:
            self.mcs_list.append(mcs)
           
            ccn1_list = []
            jux_list = []
            para_list = []
            num_wound_fib_list = []
            num_ECM_list = []
            num_pix = 0
            
            
            
            for cell in self.cell_list:
                if 50<cell.xCOM<150 and 50<cell.yCOM<200:
                    if cell.type == self.FIBROBLAST:
                        num_wound_fib_list.append(cell)
                        
                    elif cell.type == self.ECM:
                        num_ECM_list.append(cell)
               
            
            for x,y,z in self.every_pixel():
                if 50<x<150:
                    if 50<y<200:
                       cell_at_pixel = self.cell_field[x,y,z]  
                       if cell_at_pixel:
                           if cell_at_pixel.type != self.ECM:
                                num_pix += 1
                          
            
            close_rate = (num_pix/((150-50)*(200-50))) * 100
            
            for cell in self.cell_list_by_type(self.SNCMYOF):
                
                if cell.dict['snc_mech'] == 'CCN1':
                    ccn1_list.append(cell)
                    #print('ccn1')
                elif cell.dict['snc_mech'] == 'JUX':
                    jux_list.append(cell)
                    #print('jux')
                elif cell.dict['snc_mech'] == 'PARA':
                    para_list.append(cell)
                    #print('para')
                
                
            num_ccn1 = len(ccn1_list) 
            num_jux = len(jux_list)
            num_para = len(para_list)
            num_Fibroblast = len(num_wound_fib_list)
            num_Myofibroblast = len(self.cell_list_by_type(self.MYOFIBROBLAST))    
            num_Macrophage = len(self.cell_list_by_type(self.MACROPHAGE))
            num_Sncmyof = len(self.cell_list_by_type(self.SNCMYOF))
            #num_wound_fib = len(num_wound_fib_list) 
            num_ECM = len(num_ECM_list)
            
            # Make sure Secretion plugin is loaded
            # make sure this field is defined in one of the PDE solvers
            secretor_PDGF = self.get_field_secretor("PDGF")
            secretor_Proteinase = self.get_field_secretor("Proteinase")
            secretor_Inflammation_signal = self.get_field_secretor("Inflammation_signal")
            secretor_CSF = self.get_field_secretor("CSF")
            
            num_PDGF = secretor_PDGF.totalFieldIntegral()
            num_Proteinase = secretor_Proteinase.totalFieldIntegral()
            num_CSF = secretor_CSF.totalFieldIntegral()
            num_Inflammation = secretor_Inflammation_signal.totalFieldIntegral()
            
            
            self.Fib_list.append(num_Fibroblast)
            self.MyoFib_list.append(num_Myofibroblast)
            self.SncmyoFib_list.append(num_Sncmyof)
            self.Mac_list.append(num_Macrophage)
            self.ECM_list.append(num_ECM)
            self.Area_list.append(close_rate)
            self.ccn1_list.append(num_ccn1)
            self.jux_list.append(num_jux)
            self.para_list.append(num_para)
            self.PDGF_list.append(num_PDGF)
            self.CSF_list.append(num_CSF)
            self.Inflammation_list.append(num_Inflammation)
            self.Proteinase_list.append(num_Proteinase)
            
            
            self.Area_plot.add_data_point("Area", mcs, close_rate)
            self.Fib_plot.add_data_point("FibroblastCount", mcs, num_Fibroblast)
            self.MyoFib_plot.add_data_point("Myofibroblast", mcs, num_Myofibroblast)
            self.Mac_plot.add_data_point("Macrophage", mcs, num_Macrophage)
            self.SncmyoFib_plot.add_data_point("Senescent Myofibroblast", mcs, num_Sncmyof)
            self.ECM_plot.add_data_point("ECM", mcs, num_ECM)
            self.Field_plot.add_data_point("PDGF", mcs, num_PDGF)
            self.Field_plot.add_data_point("Proteinase", mcs, num_Proteinase)
            self.Field_plot.add_data_point("CSF", mcs, num_CSF)
            self.Field_plot.add_data_point("Inflammation", mcs, num_Inflammation)
            self.SncmyoFib_type_plot.add_data_point("CCN1", mcs, num_ccn1)
            self.SncmyoFib_type_plot.add_data_point("Juxtacrine", mcs, num_jux)
            self.SncmyoFib_type_plot.add_data_point("Paracrine", mcs, num_para)
            #self.Fib_wound_plot.add_data_point("WoundFibroblastCount", mcs, num_wound_fib)
            

        if mcs > 101:
            for day in self.Fib_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    Fib_AbsErr = abs(num_Fibroblast - self.Fib_dict[day]) / (max(self.Fib_dict, key = self.Fib_dict.get))
                    self.Fib_AbsErr_sum += Fib_AbsErr
                    self.Fib_AbsErr_list.append(Fib_AbsErr)
                    
            for day in self.MyoFib_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    MyoFib_AbsErr = abs(num_Myofibroblast - self.MyoFib_dict[day]) / (max(self.MyoFib_dict, key = self.MyoFib_dict.get))
                    self.MyoFib_AbsErr_sum += MyoFib_AbsErr
                    self.MyoFib_AbsErr_list.append(MyoFib_AbsErr)
                           
                    
            for day in self.SncmyoFib_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    SncmyoFib_AbsErr = abs(num_Sncmyof - self.SncmyoFib_dict[day]) / (max(self.SncmyoFib_dict, key = self.SncmyoFib_dict.get))
                    self.SncmyoFib_AbsErr_sum += SncmyoFib_AbsErr
                    self.SncmyoFib_AbsErr_list.append(SncmyoFib_AbsErr)
                           
                    
            for day in self.Mac_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    Mac_AbsErr = abs(num_Macrophage - self.Mac_dict[day]) / (max(self.Mac_dict, key = self.Mac_dict.get))
                    self.Mac_AbsErr_sum += Mac_AbsErr
                    self.Mac_AbsErr_list.append(Mac_AbsErr)
                       
                 
             
            for day in self.Area_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    Area_AbsErr = abs(close_rate - self.Area_dict[day]) / (max(self.Area_dict, key = self.Area_dict.get))
                    self.Area_AbsErr_sum += Area_AbsErr
                    self.Area_AbsErr_list.append(Area_AbsErr)
                    
            for day in self.ECM_dict:
                day_to_mcs = ((24*60*60*int(day))/param.mcs_to_sec)
                #print(day_to_mcs, 'this is day to mcs')
                if day_to_mcs == mcs:
                    ECM_AbsErr = abs(num_ECM - self.ECM_dict[day]) / (max(self.ECM_dict, key = self.ECM_dict.get))
                    self.ECM_AbsErr_sum += ECM_AbsErr
                    self.ECM_AbsErr_list.append(ECM_AbsErr)
      

    def finish(self):
        
        out_path = Path(self.output_dir).joinpath('AreaPlot.png')
        out_path_fib = Path(self.output_dir).joinpath('FibroblastPlot.png')
        out_path_myofib = Path(self.output_dir).joinpath('MyofibroblastPlot.png')
        out_path_sncmyofib = Path(self.output_dir).joinpath('SenescentMyofibroblastPlot.png')
        out_path_mac = Path(self.output_dir).joinpath('MacrophagePlot.png')
        out_path_ECM = Path(self.output_dir).joinpath('ECMPlot.png')    
        out_path_Field = Path(self.output_dir).joinpath('FieldPlot.png')
        out_path_snc_type = Path(self.output_dir).joinpath('SncType.png')
        
        
        self.Area_plot.save_plot_as_png(str(out_path), 2000, 2000)
        self.Fib_plot.save_plot_as_png(str(out_path_fib), 2000, 2000)
        self.MyoFib_plot.save_plot_as_png(str(out_path_myofib), 2000, 2000)
        self.SncmyoFib_plot.save_plot_as_png(str(out_path_sncmyofib), 2000, 2000)
        self.Mac_plot.save_plot_as_png(str(out_path_mac), 2000, 2000)
        self.ECM_plot.save_plot_as_png(str(out_path_ECM), 2000, 2000)    
        self.Field_plot.save_plot_as_png(str(out_path_Field), 2000, 2000)
        self.SncmyoFib_type_plot.save_plot_as_png(str(out_path_snc_type), 2000, 2000)
        
        
        
        plot_data_dict = {'MCS' : self.mcs_list, 'Fibroblast' : self.Fib_list, 'Myofibroblast': self.MyoFib_list, 
        'Senescent_Fibroblast' : self.SncmyoFib_list, 'Macrophage' : self.Mac_list, 'ECM' : self.ECM_list, 'Area' : self.Area_list, 
        'CCN1_mediated_senescence' : self.ccn1_list, 'Juxtacrine_senescence' : self.jux_list, 'Paracrine_senescence' : self.para_list, 
        'PDGF' : self.PDGF_list, 'CSF' : self.CSF_list, 'Inflammation' : self.Inflammation_list, 'Proteinase' : self.Proteinase_list}
        
        
        #pd.set_option('max_columns', 200)
        #pd.set_option('max_rows', 200)
        pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.width', 1000) 
        df_pop = pd.DataFrame(plot_data_dict) 
        #df.to_csv('Population_data_WH.csv')
        
        all_AbsErr_list = {'Fib_AbsErr_sum' : self.Fib_AbsErr_list, 'MyoFib_AbsErr_sum' : self.MyoFib_AbsErr_list,
        'SncmyoFib_AbsErr_sum' : self.SncmyoFib_AbsErr_list, 'Mac_AbsErr_sum' : self.Mac_AbsErr_list, 'ECM_AbsErr_sum' : self.ECM_AbsErr_list, 'Area_AbsErr_sum' : self.Area_AbsErr_list}
        
                
        AbsErr_Ave = {}
        for i in all_AbsErr_list:
            if len(all_AbsErr_list[i]) != 0:
                AbsErr_Ave[i] = sum(all_AbsErr_list[i]) / len(all_AbsErr_list[i])
        #print(AbsErr_Ave)
        print(len(AbsErr_Ave), 'this is len of abserr list')
        if len(all_AbsErr_list) != 0:
            QualVal = sum(AbsErr_Ave.values()) 
            
        print("\n\n\t\t final total abs error=",QualVal)
              
        #time stamp
        #self.timestamp = '{:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.datetime.now())
        #         self.startTime = datetime.datetime.now()
        runTime = time.time() - self.startTime
                
        headers  = "Quality_measure" + "\t"
        Line     = str(QualVal) + "\t"
        headers += "RunTime\t"
        Line    += str(int(runTime)) + "\t"
        
        for name in dir(param):
           if not name.startswith("__") and name != "math": 
               headers += name + "\t"
               Line += str(getattr(param,name)) + "\t"
        
        for key, value in AbsErr_Ave.items():
            headers += key + "\t"
            Line += str(value) + "\t"
        
              
        # Writing quality measure and population data to files ##################
        #self.file_obj_qual.write(str(headers)+"\n"+str(theLine))
        #self.file_obj.write(str(QualVal), str(all_AbsErr_sum_list))
        #self.file_obj_qual.close()
        
        output_dir = self.output_dir
        if output_dir is not None:
            #"mcs" not available in finish   
            #output_path_pop = Path(output_dir).joinpath('PopData_step_' + str(mcs).zfill(3) + '.txt')
            #output_path_qual = Path(output_dir).joinpath('QualData_step_' + str(mcs).zfill(3) + '.txt')
            output_path_pop = Path(output_dir).joinpath('PopData_final.txt')    
            output_path_qual = Path(output_dir).joinpath('QualData_final.txt')   
            with open(output_path_pop, 'w') as fout_pop:
                #print(plot_data_dict, file = fout_pop)
               fout_pop.write(str(df_pop))
            fout_pop.close()  
               
            with open(output_path_qual, 'w') as fout_qual:
            #   #fout_qual.write(str(headers)+"\n"+str(theLine))   
               fout_qual.write(str(headers)+"\n"+str(Line))   
            fout_qual.close() 
        
        '''
        #################################################################################################    
        ###### code to move the output file to the top level of the home directory for this project, 
        ###### and making it "append", and adding a time stamp to start of line.   
        SummaryFileName = "Quality_data.txt" 
        fileDir = os.path.dirname(os.path.abspath(__file__))
        SummFilePathName = fileDir+"/"+SummaryFileName
      ##fileDirUp = os.path.abspath(os.path.join(fileDir,".."))  # go up one directory
      ##SummFilePathName = fileDirUp+"/"+SummaryFileName
        print("\n\t\t SummFilePathName:",SummFilePathName,"\n")
        try:                
            self.SummaryFile_obj = open(SummFilePathName, mode='a') # append
        except IOError:
            print("Could not open file. \n",SummFilePathName,"\n in the project's directory.")
            self.stop_simulation()
        
        print(str(headers)+"\n"+str(Line)+"\n", file=self.SummaryFile_obj)  
        self.SummaryFile_obj.close()
        '''
        
        
    def on_stop(self):
        
        return


class ClearSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        return
        

    def step(self, mcs):
        if mcs >= 5:
            if mcs == 1:
                for cell in self.cell_list:
                    if cell.type == self.MACROPHAGE:
                        # macrophage chemotaxis towards inflammation field              
                        cd_Inflammation_signal = self.chemotaxisPlugin.addChemotaxisData(cell, "Inflammation_signal")
                        cd_Inflammation_signal.setLambda(120)
                        #cd_Inflammation_signal.assignChemotactTowardsVectorTypes([self.SNCMYOF])
            
            if mcs % 100 == 0:  
                
                field_Proteinase = self.field.Proteinase
                field_Inflammation_signal = self.field.Inflammation_signal
                field_PDGF = self.field.PDGF
                field_CSF = self.field.CSF
                                              
                
                for cell in self.cell_list:
                    # Proteinase-mediated ECM breakdown 
                    if cell.type == self.ECM:
                        if field_Proteinase[cell.xCOM,cell.yCOM,cell.zCOM] >= param.proteinase_thr:
                            self.delete_cell(cell)
                            field_Proteinase[cell.xCOM,cell.yCOM,cell.zCOM] = max(0, (field_Proteinase[cell.xCOM,cell.yCOM,cell.zCOM] - param.proteinase_thr))
                    
                    # Senescent cell clearance
                    if cell.type == self.SNCMYOF or cell.type == self.SNCF:
                        neighbor_list = self.get_cell_neighbor_data_list(cell)
                        common_area_with_mac = neighbor_list.common_surface_area_with_cell_types(cell_type_list=[self.MACROPHAGE])
                        chance = random.randint(1, 100)
                        if chance <= 50:
                            # during the fibrolytic phase of senescence
                            if mcs > cell.dict['snc_start'] + ((24*60*60*3)/param.mcs_to_sec): 
                                # macrophage mediated clearance
                                if common_area_with_mac > 0:
                                    self.delete_cell(cell) 
                                    
                        else:
                            # other clearance mechanisms 
                            if field_Inflammation_signal[cell.xCOM,cell.yCOM,cell.zCOM] > param.IL6*2:
                                self.delete_cell(cell) 
                    
                # macrophage clearance during late-stage healing                
                if mcs > 55000:
                    for cell in self.cell_list_by_type(self.MACROPHAGE):
                        if field_CSF[cell.xCOM, cell.yCOM, cell.zCOM] < param.CSF_thr:
                            self.delete_cell(cell)                
                            break   
                                       
            
            
            
            # Rate of macrophage and fibroblast removal 
            if mcs > round(3.3*24*param.hour_to_mcs) == 0:
                if len(self.cell_list_by_type(self.MACROPHAGE)) > 1:
                    random_mac = random.choice(list(self.cell_list_by_type(self.MACROPHAGE)))
                    self.delete_cell(random_mac)
                
                if len(self.cell_list_by_type(self.FIBROBLAST)) > 1:
                    random_fib = random.choice(list(self.cell_list_by_type(self.FIBROBLAST)))
                    self.delete_cell(random_fib)
                
                
                
                
    
    
    def finish(self):
        
        return

    def on_stop(self):
        
        return
   
   
   
class WoundSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        return
        
    def step(self, mcs):
        
        
        if mcs == 2:
            
            for cell in self.cell_list:
                if cell.type == self.FIBROBLAST:
                    cell.targetVolume = param.Fib_TarVol 
                    cell.lambdaVolume = 1.0
                    cell.targetSurface = param.Fib_TarSur 
                    cell.lambdaSurface = 1.0
                    
                elif cell.type == self.MYOFIBROBLAST:
                    cell.targetVolume = param.Myofib_TarVol 
                    cell.lambdaVolume = 2.0
                    cell.targetSurface = param.Myofib_TarSur 
                    cell.lambdaSurface = 2.0
                    
                elif cell.type == self.MACROPHAGE:
                    cell.targetVolume = param.Mac_TarVol 
                    cell.lambdaVolume = 2.0
                    cell.targetSurface = param.Mac_TarSur 
                    cell.lambdaSurface = 2.0
                           
                elif cell.type == self.SNCMYOF:
                    cell.targetVolume = param.Sncmyof_TarVol 
                    cell.lambdaVolume = 2.0
                    cell.targetSurface = param.Sncmyof_TarSur 
                    cell.lambdaSurface = 2.0
                
                elif cell.type == self.ECM:
                    cell.targetVolume = param.Ecm_TarVol 
                    cell.lambdaVolume = 2.0
                    cell.targetSurface = param.Ecm_TarSur 
                    cell.lambdaSurface = 2.0
                
                else:
                    cell.targetVolume = cell.volume
                    cell.lambdaVolume = 1.0
                    cell.targetSurface = cell.surface
                    cell.lambdaSurface = 1.0

            # chemotaxis
            for cell in self.cell_list_by_type(self.FIBROBLAST, self.MYOFIBROBLAST):
                
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "PDGF")
                cd.setLambda(60.0) 
                #cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
                
                if cell.type == self.FIBROBLAST:
                    cell.dict['deactivated_myoF'] = 0
            
            
        if mcs > 0:    
            if mcs % 100 == 0:  
                PDGF = self.field.PDGF
                # fibroblast and myofibroblast travel along ECM
                for cell in self.cell_list_by_type(self.FIBROBLAST, self.MYOFIBROBLAST):
                    neighbour_list = self.get_cell_neighbor_data_list(cell) 
                    common_area_with_ECM = neighbour_list.common_surface_area_with_cell_types(cell_type_list=[self.ECM])
                    
                    if common_area_with_ECM > 0:    
                        cd = self.chemotaxisPlugin.getChemotaxisData(cell, "PDGF")
                        cd.setLambda(120)
                        #cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
            
                
                lowLambda = 60
                for cell in self.cell_list_by_type(self.FIBROBLAST, self.MYOFIBROBLAST):
                    neighbour_data = self.get_cell_neighbor_data_list(cell)
                    for neighbour, common_area in neighbour_data:
                        if neighbour:
                            if neighbour.type == self.MACROPHAGE:
                                cd = self.chemotaxisPlugin.getChemotaxisData(cell, "PDGF")
                                cd.setLambda = lowLambda
 
                       
                        
                        
                    
    def finish(self):
        
        return

    def on_stop(self):
        
        return



class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        # parent child position will be randomised between mitosis events
        self.set_parent_child_position_flag(0) 
        
    
    def step(self, mcs):
        if mcs % 20 == 0:
            cells_to_divide = []
            mac_to_divide = []
            myof_to_divide = []
            
            # cell division 
            for cell in self.cell_list_by_type(self.FIBROBLAST):
               if cell.volume >= ((14/param.vox_to_um)**2) * 2:
                    cells_to_divide.append(cell)
                  
            for cell in self.cell_list_by_type(self.MACROPHAGE):
               if cell.volume >= ((18/param.vox_to_um)**2) * 2:
                    mac_to_divide.append(cell)
                
            for cell in self.cell_list_by_type(self.MYOFIBROBLAST):
               if cell.volume >= ((18/param.vox_to_um)**2) * 2:
                    myof_to_divide.append(cell)
                    cell.dict['myo_start'] = mcs
                    
                    
            for cell in cells_to_divide:
                self.divide_cell_along_minor_axis(cell)
                           
            for cell in mac_to_divide:
                self.divide_cell_random_orientation(cell)
        
            for cell in myof_to_divide:
                self.divide_cell_random_orientation(cell)
        
    def update_attributes(self):
        # reducing parent target volume before cloning
        self.parent_cell.targetVolume /= 2.0                
        self.parent_cell.targetSurface /= 2.0
        #self.clone_parent_2_child()            
        self.child_cell.targetVolume  = self.parent_cell.targetVolume                 
        self.child_cell.targetSurface = self.parent_cell.targetSurface                  
        self.child_cell.lambdaVolume  = self.parent_cell.lambdaVolume                 
        self.child_cell.lambdaSurface = self.parent_cell.lambdaSurface
        # for more control of what gets copied from parent to child use cloneAttributes function
        self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[]) 
   
class SecretionSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SecretionBasePy.__init__(self, frequency)
    
    
    def step(self, mcs):
          
        if mcs % 50 == 0:
            Proteinase_Field_secretor = self.get_field_secretor("Proteinase")
            PDGF_secretor = self.get_field_secretor("PDGF")
            CSF_Field_secretor = self.get_field_secretor("CSF")
            Inflammation_Field_secretor = self.get_field_secretor("Inflammation_signal")
            
            field_PDGF = self.field.PDGF
            field_CSF = self.field.CSF
            Inflammation_signal_field = self.field.Inflammation_signal
            proteinase = self.field.Proteinase
            
            for cell in self.cell_list:
                if cell.type == self.MACROPHAGE:
                    Inflammation_value = Inflammation_signal_field[cell.xCOM, cell.yCOM, cell.zCOM]
                    if Inflammation_value < param.IL6/2:
                        PDGF_secretor.secreteOutsideCellAtBoundary(cell, param.PDGF_secretion_coeff)
                        
                    neighbour_list = self.get_cell_neighbor_data_list(cell) 
                    common_area_with_ECM = neighbour_list.common_surface_area_with_cell_types(cell_type_list=[self.ECM])        
            
                    if common_area_with_ECM > param.Mac_TarSur or Inflammation_value > param.IL6/2:
                        Proteinase_Field_secretor.secreteOutsideCellAtBoundary(cell, param.Proteinase_secretion_coeff)
            
                elif cell.type == self.FIBROBLAST:
                    
                    PDGF_conc = field_PDGF[cell.xCOM, cell.yCOM, cell.zCOM]
                    if PDGF_conc > param.actF_col3:
                        CSF_Field_secretor.secreteOutsideCellAtBoundaryOnContactWith(cell, param.CSF_secretion_coeff, [self.MACROPHAGE])
                                
                elif cell.type == self.MYOFIBROBLAST:
                    
                    CSF_Field_secretor.secreteOutsideCellAtBoundaryOnContactWith(cell, param.CSF_secretion_coeff, [self.MACROPHAGE])
                    PDGF_secretor.secreteOutsideCellAtBoundary(cell, param.PDGF_secretion_coeff)
                    proteinase[cell.xCOM, cell.yCOM, cell.zCOM] = max(0, (proteinase[cell.xCOM, cell.yCOM, cell.zCOM] - param.proteinase_thr)) # TIMP
               
             
                elif cell.type == self.SNCMYOF:
                    # fibrogenic phase
                    if mcs < cell.dict['snc_start'] + ((24*60*60*3)/param.mcs_to_sec):
                        PDGF_secretor.secreteOutsideCellAtBoundary(cell, param.PDGF_secretion_coeff)
                        cell_size = int(param.Ecm_TarSur/4)
                        x = int(random.randrange(0,200-cell_size,1))
                        y = int(random.randrange(0,200-cell_size,1))
                        cell = self.cell_field[x, y, 0]
                        if not cell:
                            ecm = self.new_cell(self.ECM)
                            self.cell_field[x:x + cell_size - 1, y:y + cell_size - 1, 0] = ecm
                            ecm.targetVolume = param.Ecm_TarVol 
                            ecm.lambdaVolume = 2.0
                            ecm.targetSurface = param.Ecm_TarSur 
                            ecm.lambdaSurface = 2.0
                        
                    # fibrolytic phase        
                    elif mcs >= cell.dict['snc_start'] + ((24*60*60*3)/param.mcs_to_sec):    
                        Inflammation_Field_secretor.secreteOutsideCellAtBoundary(cell, param.Inflammation_Field_secretion_coeff)
                        Proteinase_Field_secretor.secreteOutsideCellAtBoundary(cell, param.Proteinase_secretion_coeff)
                        CSF_Field_secretor.secreteOutsideCellAtBoundary(cell, param.CSF_secretion_coeff)
                
                elif cell.type == self.SNCF:
                    Inflammation_Field_secretor.secreteOutsideCellAtBoundary(cell, param.Inflammation_Field_secretion_coeff)
                    Proteinase_Field_secretor.secreteOutsideCellAtBoundary(cell, param.Proteinase_secretion_coeff)
                    CSF_Field_secretor.secreteOutsideCellAtBoundary(cell, param.CSF_secretion_coeff)
                                
                            
                            

