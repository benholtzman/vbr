# PLOT THE GIA predictions...

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from time import time
from imp import reload  # this lets us reload a module without restarting the kernel ! very convenient now...
import pandas as pd
import sys

# our functions:
import fit_LAB_fcns as flab
reload(flab)

# ================================================
# READ PRIOR MODELS IN DIRECTLY FROM THE MATLAB STRUCTURES:

#SNAfit = True
#if SNAfit==True:
    # not yet:
    # pathtodata = './SNA_fits_Qfit/VBR_Tp1325_Z200_An/'
    # maybe this one? or the dropbox 00_Boxes ?
    #path = '/Users/ben/0_vbr_git/VBRcloset/y161210_SNA_fit/'
#path = '/Users/ben/Dropbox/0_VBR_WORK/0_y17_Projects/Boxes/2017-07-20-SNA_forGIA/'
#matobj = 'Box_2017-07-20-SNA_forGIA_VBR_py.mat'
path = '/Users/ben/0_vbr_git/VBRcloset/2018-07-02-BigBox/'
matobj = 'Box_2018-07-02-BigBox_VBR_py.mat'
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)
zLAB_obs_km_COLD = 180 # plot this line as the seismically identified plate thickness..

#elif SNAfit==False:
    #pathtodata = './TNA_fits_Qfit/VBR_Tp1425_Z100_An/'
zLAB_obs_km_HOT = 80 # plot this line as the seismically identified plate thickness..
# b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)

# -----------------
# read in the Box

t0 = time()
# when you do stuff that takes time, see how long it takes !
hoozit = sio.whosmat(path+matobj) # get a little info without loading it..
print(hoozit)
# load the mat file:
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)
box = b['Box']
print('box is a ' + str(type(box)) + ', with shape:')
print(box.shape)
n_var1 = box.shape[0]
n_var2 = box.shape[1]
print(box[0,2]._fieldnames)
print(box[0,2].info._fieldnames)
print(box[0,2].run_info._fieldnames)

zPlate_vec = box[0,0].info.var2range # this flipped from var1range
Tpot_vec = box[0,0].info.var1range
print('Tpot_vec = ', str(Tpot_vec) )
print('zPlate_vec = ', str(zPlate_vec) )

frame_test = box[0,1].Frames[-1]
# frequency band:
f_band = frame_test.VBR.input.SV.f
print('full frequency band:')
print(f_band)
# ======================================================================
# FIND THE MODELS THAT BEST FIT THE LAB DEPTH, using the right frequency band:

#f_band of interest:
# find the index corresponding to the minimum and max freq
per_bw_max = 200.
per_bw_min = 20.
f_bw_min = 1/per_bw_max
f_bw_max = 1/per_bw_min
print('subset frequency band:')
print(f_bw_min,f_bw_max)

# find the index of frequency value closest to the min and max:
# surely there is a more pythonic way to do this, but whatever.

fmin = f_band[f_band>=f_bw_min][0]
fmax = f_band[f_band<=f_bw_max][-1]
i_fmax = int(np.where(f_band==fmax)[0])
i_fmin = int(np.where(f_band==fmin)[0])
print('actual frequency band and indexes:')
print(fmin,fmax,i_fmin,i_fmax)
# np.where : https://docs.scipy.org/doc/numpy/reference/generated/numpy.where.html

# depth vector (different for each model!):
# Z_km = box[1,1].run_info.Z_km

# for debugging:
# frame = box[15,15].Frames[-1]
# frame.VBR.out.anelastic.AndradePsP.Qa


# # ================================================
# FIND BEST MODELS FOR LAB DEPTH
# # ================================================

Q_LAB = 800 ; # if it breaks, this might have to be higher--
# add a try-catch
eta_LAB = 1e24
# or can try a gradient? dQ/dz_LAB ?

# COLD PLATE:
Res_lab_Q_mat_COLD, ind_zLAB_Q_mat_COLD = flab.find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km_COLD,i_fmin,i_fmax)
Rmin = Res_lab_Q_mat_COLD.min()
print(Rmin)
ij_best_all = np.where(Res_lab_Q_mat_COLD==Rmin)
print('COLD best fit:')
print(ij_best_all)


# this might have to be defined by hand >
#ij_best = [ij_best_all[0],ij_best_all[1]]
ij_best_COLD = ij_best_all[0][0],ij_best_all[1][0] # this makes a tuple:
# print(Res_lab_Q_mat[ij_best_all][0])
print(ij_best_COLD)
print(Rmin,Res_lab_Q_mat_COLD[ij_best_COLD])
Res_lab_etaLT_mat_COLD, ind_zLAB_LT_mat_COLD = flab.find_LAB_LT_Res(box,eta_LAB,zLAB_obs_km_COLD)

# HOT PLATE:
Res_lab_Q_mat_HOT, ind_zLAB_Q_mat_HOT = flab.find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km_HOT,i_fmin,i_fmax)
Rmin = Res_lab_Q_mat_HOT.min()
ij_best_all = np.where(Res_lab_Q_mat_HOT==Rmin)
print('HOT best fit:')
print(ij_best_all)
# this might have to be defined by hand >
#ij_best = [ij_best_all[0],ij_best_all[1]]
ij_best_HOT = ij_best_all[0][0],ij_best_all[1][0] # this makes a tuple:
# print(Res_lab_Q_mat[ij_best_all][0])
print(ij_best_HOT)
print('Rmin and best fitting val from Res matrix (should be equal):')
print(Rmin,Res_lab_Q_mat_HOT[ij_best_HOT])
Res_lab_etaLT_mat_HOT, ind_zLAB_LT_mat_HOT = flab.find_LAB_LT_Res(box,eta_LAB,zLAB_obs_km_HOT)

# or pick it off Residual Map, then export
# i_best = 5 ; % z_plate, plate thickness
# j_best = 5 ; % T_pot

# save("zLAB_mat.mat", "zLAB_mat")
LAB_bestfits = pd.DataFrame(
             {"COLDplate" : ij_best_COLD,
              "HOTplate" : ij_best_HOT})
LAB_bestfits.to_pickle('LAB_bestfits.pkl')

# ===================================================================
#
# ===================================================================
# ======================
# PLOT THE CURVES !
# ======================


plt.figure(1,figsize=(9,7))
# plasma colormap
cmap = plt.cm.get_cmap('autumn')

# =======================
#f, axs = plt.subplots(1,2,figsize=(7,7))
#plt.subplot(1, 2, 1)
L1 = 0.1
B1 = 0.2
W1 = 0.25
H1 = 0.7

ax = plt.axes([L1,B1,W1,H1])

bk = [0.0,0.0,0.0,1.0]

for i_var1 in range(n_var1):
    for i_var2 in range(n_var2):
        Z_km = box[i_var1,i_var2].run_info.Z_km
        frame = box[i_var1,i_var2].Frames[-1]
        Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
        Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
        Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
        ax.plot(np.log10(Qs_mnstd[:,0]),Z_km,color='black',lw=1.0, ls='-', alpha=0.2 )

# add standard deviations to these lines !
Z_km = box[ij_best_COLD].run_info.Z_km
frame = box[ij_best_COLD].Frames[-1]
Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
ax.plot(np.log10(Qs_mnstd[:,0]),Z_km,color='blue',lw=2.0, ls='-', alpha=0.9 )

Z_km = box[ij_best_HOT].run_info.Z_km
frame = box[ij_best_HOT].Frames[-1]
Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
ax.plot(np.log10(Qs_mnstd[:,0]),Z_km,color='red',lw=2.0, ls='-', alpha=0.9 )

#axes = plt.gca() # why this way ?
xmin = 1.0
xmax = 8.0
ax.set_xlim([xmin,xmax]) # used to be axes

ymin = 0.0
ymax = 350.0
ax.set_ylim([ymin,ymax]) # used to be axes'

# add COLD plate Q LAB line:
Z_LAB_LT_km_COLD = Z_km[int(ind_zLAB_Q_mat_COLD[ij_best_COLD])]
zLAB_LT_line_COLD = Z_LAB_LT_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_COLD,color='blue',lw=2.0, linestyle=':', alpha=0.5)
# add HOT plate Q LAB line:
Z_LAB_LT_km_HOT = Z_km[int(ind_zLAB_LT_mat_COLD[ij_best_HOT])]
zLAB_LT_line_HOT = Z_LAB_LT_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_HOT,color='red',lw=2.0, linestyle=':', alpha=0.5)

# add COLD plate LT viscosity LAB line:
Z_LAB_LT_km_COLD = Z_km[int(ind_zLAB_LT_mat_COLD[ij_best_COLD])]
zLAB_LT_line_COLD = Z_LAB_LT_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_COLD,color='blue',lw=2.0, linestyle='--', alpha=0.5)
# add HOT plate LT viscosity LAB line:
Z_LAB_LT_km_HOT = Z_km[int(ind_zLAB_LT_mat_COLD[ij_best_HOT])]
zLAB_LT_line_HOT = Z_LAB_LT_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_HOT,color='red',lw=2.0, linestyle='--', alpha=0.5)


# add (OBSERVED!) COLD plate thickness line :
zLAB_line_COLD = zLAB_obs_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_COLD,color='blue',lw=6.0, linestyle='-', alpha=0.2)
# add HOT plate thickness line:
zLAB_line_HOT = zLAB_obs_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_HOT,color='red',lw=6.0, linestyle='-', alpha=0.2)

plt.title('Fit the LAB with Qs')
plt.xlabel('log Qs')
plt.ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax.invert_yaxis()

# =================================
# plot the residuals for Q_LAB HOT
L2 = 0.45
B2 = B1+0.4
W2 = 0.3
H2 = W2

ax2 = plt.axes([L2,B2,W2,H2])
ax2.imshow(Res_lab_Q_mat_HOT, cmap=plt.cm.gray)
ax2.scatter(ij_best_HOT[1],ij_best_HOT[0], color='red', s=8)
ax2.set_title('Res: Z_LAB (Hot plate)')
ax2.set_xticklabels(zPlate_vec)
ax2.set_xlabel('z_plate')
ax2.set_yticklabels(Tpot_vec)
ax2.set_ylabel('T_pot')

# =================================
# plot the residuals for Q_LAB COLD
L3 = L2
B3 = B1
W3 = W2
H3 = H2

ax3 = plt.axes([L3,B3,W3,H3])
ax3.imshow(Res_lab_Q_mat_COLD, cmap=plt.cm.gray)
ax3.scatter(ij_best_COLD[1],ij_best_COLD[0], color='blue', s=8)

#ax3.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax3.set_title('Res: Z_LAB (Cold plate)')
ax3.set_xticklabels(zPlate_vec)
ax3.set_xlabel('z_plate')
ax3.set_yticklabels(Tpot_vec)
ax3.set_ylabel('T_pot')

# =================================


plt.show()

#import GIA_transient_plates_mat_v2

# =============
#sys.exit()
# =============

# ============
# sys.exit()
# ============
