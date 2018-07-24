# PLOT THE GIA predictions...

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from time import time
from imp import reload  # this lets us reload a module without restarting the kernel !
import pandas as pd
import sys

# our functions:
import fit_LAB_fcns as flab
reload(flab)

# ================================================
# READ PRIOR MODELS IN FROM THE MATLAB STRUCTURES:

#path = '/Users/ben/Dropbox/0_VBR_WORK/0_y17_Projects/Boxes/2017-07-20-SNA_forGIA/'
#matobj = 'Box_2017-07-20-SNA_forGIA_VBR_py.mat'
#path = '/Users/ben/0_vbr_git/VBRcloset/2018-07-02-BigBox/'
path = '/Users/ben/0_vbr_projects/VBRcloset/2018-07-02-BigBox/'
matobj = 'Box_2018-07-02-BigBox_VBR_py.mat'

zLAB_obs_km_COLD = 180 # plot this line as the seismically identified plate thickness..
zLAB_obs_km_HOT = 80 # plot this line as the seismically identified plate thickness..

Vs_adavg_obs_COLD = 4.7
Vs_adavg_obs_HOT = 4.6
# -----------------
# read in the Box

t0 = time()
# load the mat file:
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)
#hoozit = sio.whosmat(path+matobj) # get a little info without loading it..
#print(hoozit)
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
# frequency band:
f_band = box[0,1].Frames[-1].VBR.input.SV.f
print('full frequency band:')
print(f_band[0],f_band[-1])
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

t1 = time()
print(str(t1-t0)+' seconds to load')

# # ================================================
# FIND BEST MODELS FOR LAB DEPTH
# # ================================================

Q_LAB = 800 ; # if it breaks, this might have to be higher--
# add a try-catch
eta_LAB = 1e24
# or can try a gradient? dQ/dz_LAB ?

# COLD PLATE:
Res_lab_Q_mat_COLD, ind_zLAB_Q_mat_COLD, Z_LAB_Q_mat_COLD = flab.find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km_COLD,i_fmin,i_fmax)
Rmin = Res_lab_Q_mat_COLD.min()
print(Rmin)
ij_best_all = np.where(Res_lab_Q_mat_COLD==Rmin)
print('COLD best fit:')
#print(ij_best_all)
ij_best_COLD = ij_best_all[0][0],ij_best_all[1][0] # this makes a tuple:
print(ij_best_COLD)
#print(Rmin,Res_lab_Q_mat_COLD[ij_best_COLD])
#Res_lab_etaLT_mat_COLD, ind_zLAB_LT_mat_COLD = flab.find_LAB_LT_Res(box,eta_LAB,zLAB_obs_km_COLD)

# HOT PLATE:
Res_lab_Q_mat_HOT, ind_zLAB_Q_mat_HOT, Z_LAB_Q_mat_HOT = flab.find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km_HOT,i_fmin,i_fmax)
Rmin = Res_lab_Q_mat_HOT.min()
ij_best_all = np.where(Res_lab_Q_mat_HOT==Rmin)
print('HOT best fit:')
#print(ij_best_all)
ij_best_HOT = ij_best_all[0][0],ij_best_all[1][0] # this makes a tuple:
print(ij_best_HOT)
print('Rmin and best fitting val from Res matrix (should be equal):')
#print(Rmin,Res_lab_Q_mat_HOT[ij_best_HOT])
#Res_lab_etaLT_mat_HOT, ind_zLAB_LT_mat_HOT = flab.find_LAB_LT_Res(box,eta_LAB,zLAB_obs_km_HOT)

# do change i_fmin, i_fmax !
# find average Vs best fitting (COLD):
Res_Vs_adavg_mat_COLD, Vs_adavg_mat_COLD = flab.find_Vs_adavg_Res(box,Vs_adavg_obs_COLD,i_fmin,i_fmax)

# find average Vs best fitting (HOT):
Res_Vs_adavg_mat_HOT, Vs_adavg_mat_HOT = flab.find_Vs_adavg_Res(box,Vs_adavg_obs_HOT,i_fmin,i_fmax)

# FIND THE JOINT BEST FITTING MODEL !
# normalize and multiply pointwise the residual matrixes?
Res_Vs_N_HOT = Res_Vs_adavg_mat_HOT/np.max(Res_Vs_adavg_mat_HOT)
Res_Vs_N_COLD = Res_Vs_adavg_mat_COLD/np.max(Res_Vs_adavg_mat_COLD)

Res_zLAB_N_HOT = Res_lab_Q_mat_HOT/np.max(Res_lab_Q_mat_HOT)
Res_zLAB_N_COLD = Res_lab_Q_mat_COLD/np.max(Res_lab_Q_mat_COLD)

P_Vs_HOT = (1-Res_Vs_N_HOT)/np.max((1-Res_Vs_N_HOT))
P_zPlate_HOT = (1-Res_zLAB_N_HOT)/np.max((1-Res_zLAB_N_HOT))
P_Vs_COLD = (1-Res_Vs_N_COLD)/np.max((1-Res_Vs_N_COLD))
P_zPlate_COLD= (1-Res_zLAB_N_COLD)/np.max((1-Res_zLAB_N_COLD))

P_JOINT_HOT = (P_Vs_HOT*P_zPlate_HOT)**2
P_JOINT_COLD = (P_Vs_COLD*P_zPlate_COLD)**2


# Rmin = Res_JOINT_HOT.min()
# ij_best_all = np.where(Res_JOINT_HOT==Rmin)
# ij_best_HOT = ij_best_all[0][0],ij_best_all[1][0]
#
# Rmin = Res_JOINT_COLD.min()
# ij_best_all = np.where(Res_JOINT_COLD==Rmin)
# ij_best_COLD = ij_best_all[0][0],ij_best_all[1][0]

Pmax = P_JOINT_HOT.max()
ij_best_all = np.where(P_JOINT_HOT==Pmax)
ij_best_HOT = ij_best_all[0][0],ij_best_all[1][0]

Pmax = P_JOINT_COLD.max()
ij_best_all = np.where(P_JOINT_COLD==Pmax)
ij_best_COLD = ij_best_all[0][0],ij_best_all[1][0]


# ADD IN HERE THE JOINT FITTING WITH Vs !

LAB_bestfits = pd.DataFrame(
             {"COLDplate" : ij_best_COLD,
              "HOTplate" : ij_best_HOT})
LAB_bestfits.to_pickle('LAB_bestfits.pkl')

t2 = time()
print(str(t2-t1)+' seconds to do all the finding')

# NOW HERE RUN THE AVERAGE Vs at freq band of Vs surface waves (use kernels to WEIGHT freq to depth !
# ( a job for emily this !) !

# ===================================================================
#
# ===================================================================
# ======================
# PLOT THE CURVES !
# ======================
#
# import fit_LAB_PLOT as flabplot
# reload(flabplot)
# flabplot


plt.figure(1,figsize=(9,7))
cmap = plt.cm.get_cmap('autumn')

# LAYOUT:
# DEPTH PLOTS:
# Q(Z) plot:
L1 = 0.06
B1 = 0.2
W1 = 0.18
H1 = 0.7

# Vs(Z) plot:
L2 = L1+W1+0.01
B2 = B1
W2 = W1
H2 = H1


# =======================
#f, axs = plt.subplots(1,2,figsize=(7,7))

# Q(Z) PLOT
ax = plt.axes([L1,B1,W1,H1])

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
Z_LAB_Q_km_COLD = Z_LAB_Q_mat_COLD[ij_best_COLD]
zLAB_Q_line_COLD = Z_LAB_Q_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_Q_line_COLD,color='blue',lw=2.0, linestyle=':', alpha=0.5)
# add HOT plate Q LAB line:
Z_LAB_Q_km_HOT = Z_LAB_Q_mat_HOT[ij_best_HOT]
zLAB_Q_line_HOT = Z_LAB_Q_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_Q_line_HOT,color='red',lw=2.0, linestyle=':', alpha=0.5)

# add (OBSERVED!) COLD plate thickness line :
zLAB_line_COLD = zLAB_obs_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_COLD,color='blue',lw=7.0, linestyle='-', alpha=0.2)
# add HOT plate thickness line:
zLAB_line_HOT = zLAB_obs_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_HOT,color='red',lw=7.0, linestyle='-', alpha=0.2)

plt.title('Fit LAB depth with $Q_{LAB}$.')
plt.xlabel('log $Q$')
ax.set_yticklabels([])
plt.ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax.invert_yaxis()



# =======================
# PLOT THE Vs STRUCTURES

ax2 = plt.axes([L2,B2,W2,H2])

for i_var1 in range(n_var1):
    for i_var2 in range(n_var2):
        Z_km = box[i_var1,i_var2].run_info.Z_km
        frame = box[i_var1,i_var2].Frames[-1]
        Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va
        Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
        Vs_mnstd = flab.zip_meanstd_fband(Vs_f_band)
        ax2.plot((Vs_mnstd[:,0])/1e3,Z_km,color='black',lw=1.0, ls='-', alpha=0.2 )

# add standard deviations to these lines !
Z_km = box[ij_best_COLD].run_info.Z_km
frame = box[ij_best_COLD].Frames[-1]
Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va
Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
Vs_mnstd = flab.zip_meanstd_fband(Vs_f_band)
ax2.plot((Vs_mnstd[:,0])/1e3,Z_km,color='blue',lw=2.0, ls='-', alpha=0.9 )

Z_km = box[ij_best_HOT].run_info.Z_km
frame = box[ij_best_HOT].Frames[-1]
Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va
Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
Vs_mnstd = flab.zip_meanstd_fband(Vs_f_band)
ax2.plot((Vs_mnstd[:,0])/1e3,Z_km,color='red',lw=2.0, ls='-', alpha=0.9 )

#axes = plt.gca() # why this way ?
xmin = 4.2
xmax = 5.05
ax2.set_xlim([xmin,xmax]) # used to be axes

ymin = 0.0
ymax = 350.0
ax2.set_ylim([ymin,ymax]) # used to be axes'

# Vs lines MEASURED:

Vs_MEASRD_line_COLD = Vs_adavg_obs_COLD*np.ones(len(Z_km))
ax2.plot(Vs_MEASRD_line_COLD,np.linspace(ymin,ymax,len(Z_km)),color='blue',lw=5.0, linestyle='-', alpha=0.2)
# add HOT plate Vs line:
Vs_MEASRD_line_HOT = Vs_adavg_obs_HOT*np.ones(len(Z_km))
ax2.plot(Vs_MEASRD_line_HOT,np.linspace(ymin,ymax,len(Z_km)),color='red',lw=5.0, linestyle='-', alpha=0.2)
# Vs lines CALCULATED
# add COLD plate Vs line:
# Res_Vs_adavg_mat_COLD, Vs_adavg_mat_COLD
Vs_COLD = Vs_adavg_mat_COLD[ij_best_COLD]
Vs_line_COLD = Vs_COLD*np.ones(len(Z_km))
ax2.plot(Vs_line_COLD,np.linspace(ymin,ymax,len(Z_km)),color='blue',lw=2.0, linestyle=':', alpha=0.5)
# add HOT plate Vs line:
Vs_HOT = Vs_adavg_mat_HOT[ij_best_HOT]
Vs_line_HOT = Vs_HOT*np.ones(len(Z_km))
ax2.plot(Vs_line_HOT,np.linspace(ymin,ymax,len(Z_km)),color='red',lw=2.0, linestyle=':', alpha=0.5)



# # add (OBSERVED!) COLD plate thickness line :
# zLAB_line_COLD = zLAB_obs_km_COLD*np.ones(len(Z_km))
# ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_COLD,color='blue',lw=7.0, linestyle='-', alpha=0.2)
# # add HOT plate thickness line:
# zLAB_line_HOT = zLAB_obs_km_HOT*np.ones(len(Z_km))
# ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_HOT,color='red',lw=7.0, linestyle='-', alpha=0.2)

plt.title('Fit $V_s$ in asth.')
plt.xlabel('$V_s$ [km/s]')
ax2.set_yticklabels([])
#plt.autoscale(enable=True, axis='y', tight=True)
ax2.invert_yaxis()


# ====================================================
# RESIDUAL PLOTS:
hdel = 0.01
vdel = 0.02
hdX = 14*hdel

rhL1 = L2+W2+hdel
rhB1 = B1+0.4
rhW1 = 0.3
rhH1 = rhW1

rhL2 = rhL1+hdX+hdel
rhB2 = rhB1
rhW2 = rhW1
rhH2 = rhW2

rhL3 = rhL2+(hdX+hdel)
rhB3 = rhB1
rhW3 = rhW1
rhH3 = rhW2

# =============
rbL1 = L2+W2+hdel
rbB1 = B1
rbW1 = 0.3
rbH1 = rbW1

rbL2 = rbL1+hdX+hdel
rbB2 = rbB1
rbW2 = rbW1
rbH2 = rbW2

rbL3 = rbL2+(hdX+hdel)
rbB3 = rbB1
rbW3 = rbW1
rbH3 = rbW2

ddd = 0.0000000001
# =================================
# plot the residuals for Q_LAB HOT

ax2 = plt.axes([rhL1,rhB1,rhW1,rhH1])
xy_rangelist = [zPlate_vec[0],zPlate_vec[-1],Tpot_vec[-1],Tpot_vec[0]] # flip vertical directions?
ax2.imshow(np.log10(Res_lab_Q_mat_HOT), cmap=plt.cm.gray, extent=xy_rangelist) #, aspect="auto")
ax2.scatter(zPlate_vec[ij_best_HOT[1]], Tpot_vec[ij_best_HOT[0]], color='red', s=10)

ax2.set_title('Res: $Z_{LAB}$ (HOT)')
#ax2.set_xticklabels(zPlate_vec)
#ax2.xaxis.set_major_locator(plt.LinearLocator(5))
ax2.set_xlabel('$z_{plate}$')
#ax2.set_yticklabels(Tpot_vec)
#ax2.yaxis.set_major_locator(plt.LinearLocator(5))
ax2.set_ylabel('$T_{pot}$ [C]')


# =================================
# plot the residuals for Vs_adavg HOT

ax3 = plt.axes([rhL2,rhB2,rhW2,rhH2])
ax3.imshow(np.log10(Res_Vs_adavg_mat_HOT), cmap=plt.cm.gray, extent=xy_rangelist)
ax3.scatter(zPlate_vec[ij_best_HOT[1]], Tpot_vec[ij_best_HOT[0]], color='red', s=10)

ax3.set_title('Res: $V_s$')
#ax4.set_xticklabels(zPlate_vec)
ax3.set_xlabel('$z_{plate}$')
#ax4.set_yticklabels(Tpot_vec)
#ax4.set_ylabel('T_pot')
ax3.set_yticklabels([])

# =================================
# plot the residuals for Vs_adavg HOT


ax4 = plt.axes([rhL3,rhB3,rhW3,rhH3])
#ax4.imshow(np.log10(Res_JOINT_HOT), cmap=plt.cm.gray, extent=xy_rangelist)
ax4.imshow(P_JOINT_HOT, cmap=plt.cm.gray, extent=xy_rangelist)
ax4.scatter(zPlate_vec[ij_best_HOT[1]], Tpot_vec[ij_best_HOT[0]], color='red', s=10)

ax4.set_title('Joint Prob.')
#ax4.set_xticklabels(zPlate_vec)
ax4.set_xlabel('$z_{plate}$')
#ax4.set_yticklabels(Tpot_vec)
#ax4.set_ylabel('T_pot')
ax4.set_yticklabels([])


# BOTTOM ROW: COLD:
# =================================
# plot the residuals for Q_LAB COLD

ax5 = plt.axes([rbL1,rbB1,rbW1,rbH1])
ax5.imshow(np.log10(Res_lab_Q_mat_COLD), cmap=plt.cm.gray, extent=xy_rangelist)
ax5.scatter(zPlate_vec[ij_best_COLD[1]], Tpot_vec[ij_best_COLD[0]], color='blue', s=10)

#ax3.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax5.set_title('Res: $Z_{LAB}$ (COLD)')
#ax3.set_xticklabels(zPlate_vec)
ax5.set_xlabel('$z_{plate}$')
#ax3.set_yticklabels(Tpot_vec)
ax5.set_ylabel('$T_{pot}$ [C]')

# BOTTOM ROW: COLD
# =================================
# plot the residuals for Vs_adavg COLD

ax5 = plt.axes([rbL2,rbB2,rbW2,rbH2])
ax5.imshow(np.log10(Res_Vs_adavg_mat_COLD), cmap=plt.cm.gray, extent=xy_rangelist)
ax5.scatter(zPlate_vec[ij_best_COLD[1]], Tpot_vec[ij_best_COLD[0]], color='blue', s=10)

#ax5.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax5.set_title('Res: $V_s$')
#ax5.set_xticklabels(zPlate_vec)
ax5.set_xlabel('$z_{plate}$')
ax5.set_yticklabels([])
#ax5.set_ylabel('T_pot')

# =================================
# plot the residuals for JOINT PROB

ax5 = plt.axes([rbL3,rbB3,rbW3,rbH3])
#ax5.imshow(np.log10(Res_JOINT_COLD), cmap=plt.cm.gray, extent=xy_rangelist)
ax5.imshow(P_JOINT_COLD, cmap=plt.cm.gray, extent=xy_rangelist)
ax5.scatter(zPlate_vec[ij_best_COLD[1]], Tpot_vec[ij_best_COLD[0]], color='blue', s=10)

#ax5.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax5.set_title('Joint prob.')
#ax5.set_xticklabels(zPlate_vec)
ax5.set_xlabel('$z_{plate}$')
ax5.set_yticklabels([])
#ax5.set_ylabel('T_pot')

figname = 'fitLAB_YTmaxwell.png'
#figname = 'fitLAB_YTtempdep.png'
#figname = 'fitLAB_JFandrade.png'
#figname = 'fitLAB_JFeburgers.png'
plt.savefig(figname)
plt.show()

# ============
# sys.exit()
# ============
