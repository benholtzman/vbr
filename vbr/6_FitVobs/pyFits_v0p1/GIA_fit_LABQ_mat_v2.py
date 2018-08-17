# FIT the measured LAB and Vs to find reasonable models
# need to rewrite this as a class..., after working out the best criterion for LAB finding.
# python GIA_fit_LABQ_mat_v2.py 'path/to/box' 'matlab_box.mat'
#
# e.g., python GIA_fit_LABQ_mat_v2.py ~/Dropbox/Shared_Release/00_Test_Boxes/VBRcloset/2018-07-02-BigBox Box_2018-07-02-BigBox_VBR_py.ma

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from time import time
from imp import reload  # this lets us reload a module without restarting the kernel !
import pandas as pd
import sys

# our functions:
from pyBox import pyBox
import fit_LAB_fcns as flab
import fit_LAB_plt_fcns as fplt
reload(flab)

# ================================================
# set pyBox attribute dictionary
attrs={}
if len(sys.argv) > 1:
    attrs['path'] = sys.argv[1]
else:
    attrs['path']='/Users/ben/0_vbr_projects/VBRcloset/2018-07-02-BigBox'
if len(sys.argv) > 2:
    attrs['boxobj'] = sys.argv[2]
else:
    attrs['boxobj']='Box_2018-07-02-BigBox_VBR_py.mat'


# settings for fitting
per_bw_max = 30. # max period to use for fitting  [s]
per_bw_min = 10. # min period to use for fitting [s]

# observations to fit to
# zLAB: seismically identified plate thickness for plot [km]
# Vs_adavg: average Vs in adiabatic portion [km/s]
obs={'hot':{'zLAB':80,'Vs_adavg':4.6},
    'cold':{'zLAB':180,'Vs_adavg':4.7},
    'Q_LAB':800,'eta_LAB':1e24}

# ================================================
# READ PRIOR MODELS IN FROM THE MATLAB STRUCTURES:
# initialize pyBox class and load the box in
times={}
times['t0'] = time()
BoxObj=pyBox(attrs)
BoxObj.loadBox()
box=BoxObj.box

# pull some vectors to make processing easier
zPlate_vec = BoxObj.var2['range']
Tpot_vec = BoxObj.var1['range']
f_band = BoxObj.f_band
print("\n")
print('Tpot_vec = ', str(Tpot_vec) )
print('zPlate_vec = ', str(zPlate_vec) )
print('full frequency band:')
print(f_band[0],f_band[-1])

# ======================================================================
# find the values and indeces for the frequency range of interest
fmin,fmax,i_fmin,i_fmax=BoxObj.f_band_limits(per_bw_min,per_bw_max,'period')
times['t1'] = time()

# ================================================
# FIND BEST MODELS FOR LAB DEPTH
# ================================================
fits={'hot':{},'cold':{}}
for temp in ['hot','cold']:
    print("\n"+temp+' plate best fit:')
    fits[temp]['Res_LAB'],fits[temp]['ind'],fits[temp]['LAB'] =flab.find_LAB_Q_Res(box,obs['Q_LAB'],obs[temp]['zLAB'],i_fmin,i_fmax)
    fits[temp]['Rmin']=fits[temp]['Res_LAB'].min()
    print(fits[temp]['Rmin'])
    ij_best_all=np.where(fits[temp]['Res_LAB']==fits[temp]['Rmin'])
    fits[temp]['ij_best']=ij_best_all[0][0],ij_best_all[1][0]

    # find average Vs best fitting
    fits[temp]['Res_Vs_adavg_mat'],fits[temp]['Vs_adavg_mat']= flab.find_Vs_adavg_Res(box,obs[temp]['Vs_adavg'],i_fmin,i_fmax)

    # FIND THE JOINT BEST FITTING MODEL !
    # normalize and multiply pointwise the residual matrixes?
    maxval=np.max(fits[temp]['Res_Vs_adavg_mat'])
    fits[temp]['Res_Vs_N']=fits[temp]['Res_Vs_adavg_mat']/maxval
    maxval=np.max(fits[temp]['Res_LAB'])
    fits[temp]['Res_zLAB_N']=fits[temp]['Res_LAB']/maxval

    # see manual... maybe not implemented correctly. Menke book, Ch 11 !
    fits[temp]['P_Vs'] = (2*np.pi*fits[temp]['Res_Vs_N'])**(-0.5)*np.exp(-0.5*fits[temp]['Res_Vs_N'])
    fits[temp]['P_zPlate'] = (2*np.pi*fits[temp]['Res_zLAB_N'])**(-0.5)*np.exp(-0.5*fits[temp]['Res_zLAB_N'])
    fits[temp]['P_JOINT']=(fits[temp]['P_Vs']*fits[temp]['P_zPlate'])**2

    Pmax = fits[temp]['P_JOINT'].max()
    ij_best_all = np.where(fits[temp]['P_JOINT']==Pmax)
    fits[temp]['ij_best']=ij_best_all[0][0],ij_best_all[1][0]

# ADD IN HERE THE JOINT FITTING WITH Vs !
ij_best_COLD=fits['cold']['ij_best']
ij_best_HOT=fits['hot']['ij_best']
LAB_bestfits = pd.DataFrame(
     {"COLDplate" : ij_best_COLD, 'cold_indLAB': fits['cold']['ind'][ij_best_COLD],
      "HOTplate" : ij_best_HOT, 'hot_indLAB': fits['hot']['ind'][ij_best_HOT]})
LAB_bestfits.to_pickle('LAB_bestfits.pkl')

times['t2'] = time()
print(str(times['t2']-times['t1'])+' seconds to do all the finding')

# NOW HERE RUN THE AVERAGE Vs at freq band of Vs surface waves (use kernels to WEIGHT freq to depth !
# ( a job for emily this !) !

# ======================
# PLOT THE CURVES !
# ======================

Layout=fplt.buildLayout() # Define locations, widths & heights for all plots
plt.figure(1,figsize=(9,7))
# Q(Z) PLOT
ax = fplt.plot_Q_profiles(Layout['Qz'],BoxObj,fits,obs,i_fmin,i_fmax)
# Vs(Z) PLOT
ax2 = fplt.plot_Vs_profiles(Layout['Vsz'],BoxObj,fits,obs,i_fmin,i_fmax)
# Residuals TOP row: HOT
ax3,ax4,ax5=fplt.residPlots(Layout,'rh',Tpot_vec,zPlate_vec,fits,'hot','red')
# Residuals BOTTOM row: COLD
ax6,ax7,ax8=fplt.residPlots(Layout,'rb',Tpot_vec,zPlate_vec,fits,'cold','blue')

figname = 'fitLAB_YTmaxwell_v2.png'
#figname = 'fitLAB_YTtempdep.png'
#figname = 'fitLAB_JFandrade.png'
#figname = 'fitLAB_JFeburgers.png'
plt.savefig(figname)
plt.show()

# ============
# sys.exit()
# ============
