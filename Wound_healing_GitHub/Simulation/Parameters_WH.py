
# 1 voxel (side) = 3 microns
# 1 MCS = 27 seconds  # wound healing process occurs over ~10 days 

# PSO
actF_col3              = 8.1044e-02 # min:0.05 max:2.0## 
myoF_to_sncF           = 4.2876e+00 # min:1.0 max:5.0## 
proteinase_thr         = 2.0926e-02 # min:0.02 max:2.0 ***** ## 
Col_myoF_deact         = 8.0205e+00 # min:5.0 max:11.0## 
IL6                    = 2.3734e+00 # min:1.0 max:3.0## 
CSF_thr                = 3.2847e+00 # min:2.0 max:5.0##
#


Fib_actF_myoF_chance = .20
senescence_percentage_param = .15


vox_to_um = 3
mcs_to_sec = 27 
Fib_TarVol = (14/vox_to_um)**2
Fib_TarSur = (14/vox_to_um)*4
Myofib_TarVol = (18/vox_to_um)**2
Myofib_TarSur = (18/vox_to_um)*4
Mac_TarVol = (18/vox_to_um)**2
Mac_TarSur = (18/vox_to_um)*4
Kera_TarVol = (13/vox_to_um)**2
Kera_TarSur = (13/vox_to_um)*4
Sncmyof_TarVol = (14/vox_to_um)**2
Sncmyof_TarSur = (14/vox_to_um)*4
Ecm_TarVol = (7/vox_to_um)**2
Ecm_TarSur = (7/vox_to_um)*4


hour_to_mcs = 3600/mcs_to_sec
min_to_mcs = 60/mcs_to_sec
#mod_step = 10

PDGF_secretion_coeff = 1 
Proteinase_secretion_coeff = 1 
CSF_secretion_coeff = 1 
Inflammation_Field_secretion_coeff = 1 


PDGF_conc_halfmaxgrowth = 4.6*(10**-10) 
CSF_conc_halfmaxgrowth = 1.3*(10**-9) 
Gmax_Fib = 0.9 
Gmax_Mac = 0.8 





