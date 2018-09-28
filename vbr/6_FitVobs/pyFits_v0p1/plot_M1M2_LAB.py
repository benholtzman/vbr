# plot storage and loss moduli as a function of frequency,
# and different criteria for defining the LAB:


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

path = '/Users/ben/Dropbox/00_Test_Boxes/VBRcloset/2018-07-02-BigBox/'
matobj = 'Box_2018-07-02-BigBox_VBR_py.mat'

zLAB_obs_km_COLD = 180 # plot this line as the seismically identified plate thickness..
zLAB_obs_km_HOT = 80 # plot this line as the seismically identified plate thickness..

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



per_bw_max = 30. # 200
per_bw_min = 10. # 20
f_bw_min = 1/per_bw_max
f_bw_max = 1/per_bw_min
print('subset frequency band:')
print(f_bw_min,f_bw_max)
fmin = f_band[f_band>=f_bw_min][0]
fmax = f_band[f_band<=f_bw_max][-1]
i_fmax = int(np.where(f_band==fmax)[0])
i_fmin = int(np.where(f_band==fmin)[0])
print('BODY WAVE frequency band and indexes:')
print(fmin,fmax,i_fmin,i_fmax)

per_sw_max = 250
per_sw_min = 100
f_sw_min = 1/per_sw_max
f_sw_max = 1/per_sw_min
print('subset frequency band:')
print(f_sw_min,f_sw_max)
fmin = f_band[f_band>=f_sw_min][0]
fmax = f_band[f_band<=f_sw_max][-1]
i_fmax = int(np.where(f_band==fmax)[0])
i_fmin = int(np.where(f_band==fmin)[0])
print('SURFACE WAVE frequency band and indexes:')
print(fmin,fmax,i_fmin,i_fmax)
# ===========================================================
# Extract the properties from the box
# ===========================================================
LAB_bestfits = pd.read_pickle('LAB_bestfits.pkl')
#print(LAB_bestfits.HOTplate)
i_var1 = 5 #LAB_bestfits.HOTplate[0]
i_var2 = 5 #LAB_bestfits.COLDplate[1]

Z_km = box[i_var1,i_var2].run_info.Z_km
frame = box[i_var1,i_var2].Frames[-1]
z_plate = float(box[i_var1,i_var2].info.var2val)
# kill this stuff !

J1_f_AndPsP_mat = frame.VBR.out.anelastic.AndradePsP.J1
J2_f_AndPsP_mat = frame.VBR.out.anelastic.AndradePsP.J2

J1_f_YTmax_mat = frame.VBR.out.anelastic.YT_maxwell.J1
J2_f_YTmax_mat = frame.VBR.out.anelastic.YT_maxwell.J2

Mf_YTmax_mat = frame.VBR.out.anelastic.YT_maxwell.M
Mf_AndPsP_mat = frame.VBR.out.anelastic.AndradePsP.Ma

Qf_YTmax_mat = frame.VBR.out.anelastic.YT_maxwell.Q
Qf_AndPsP_mat = frame.VBR.out.anelastic.AndradePsP.Qa

Vf_YTmax_mat = frame.VBR.out.anelastic.YT_maxwell.V
Vf_AndPsP_mat = frame.VBR.out.anelastic.AndradePsP.Va

eta_ss_LH_mat = frame.VBR.out.viscous.LH2012.eta_total
eta_ss_HK_mat = frame.VBR.out.viscous.HK2003.eta_total

# print(frame.VBR.out.anelastic._fieldnames)
# print(frame.VBR.out.anelastic.AndradePsP._fieldnames)
# print(frame.VBR.out.anelastic.YT_maxwell._fieldnames)
# Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
# Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
# Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
# Eta_app_AndPsP = M2_f_AndPsP_mat[]/(2*np.pi*f)

# ===========================================================
# Find the LABs
# ===========================================================

#Res_lab_Q_mat_HOT, ind_zLAB_Q_mat_HOT, Z_LAB_Q_mat_HOT = find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km_HOT,i_fmin,i_fmax)

nn_pts = 1000
Z_km_interp = np.linspace(Z_km[0],Z_km[-1],nn_pts)
ind_LAB_YTmax_fmat = [] #(f_band.size)
ind_LAB_AndPsP_fmat = [] #np.zeros(f_band.size)
fQ_LAB = 20

for ifreq in range(f_band.size):
    Qs_mat = Qf_YTmax_mat
    Qz_interp = np.interp(Z_km_interp,Z_km,Qs_mat[:,ifreq])
    Z_LAB_Qs_km, ind_Qlab = flab.find_LAB_Q_one(Qz_interp,Z_km_interp,z_plate,fQ_LAB,zLAB_obs_km_HOT) #Qz,Z_km,z_plate,fQ_LAB,zLAB_obs_km
    ind_LAB_YTmax_fmat.append(ind_Qlab)

    Qs_mat = Qf_AndPsP_mat
    Qz_interp = np.interp(Z_km_interp,Z_km,Qs_mat[:,ifreq])
    Z_LAB_Qs_km, ind_Qlab = flab.find_LAB_Q_one(Qz_interp,Z_km_interp,z_plate,fQ_LAB,zLAB_obs_km_HOT)
    ind_LAB_AndPsP_fmat.append(ind_Qlab)


# ===========================================================
# PLOTTING
# ===========================================================
# PLOT M1 and M2, as a function of frequency,
# along with possible criteria for defining the LAB
# using Q(f) and M(f), and making sure they are equivalent.

# how different will LAB estimate be for body and surface waves?

plt.figure(1,figsize=(9,7))
cmap = plt.cm.get_cmap('autumn')

# LAYOUT:
# DEPTH PLOTS:
# M1(Z) plot:
L1 = 0.06
B1 = 0.2
W1 = 0.16
H1 = 0.7

dx = 0.02
# M2(Z) plot:
L2 = L1+W1+dx
B2 = B1
W2 = W1
H2 = H1

# M(f)(Z) plot:
L3 = L2+W1+dx
B3 = B1
W3 = W1
H3 = H1

# Q(f)(Z) plot:
L4 = L3+W1+dx
B4 = B1
W4 = W1
H4 = H1

# Vs(f)(Z) plot:
L5 = L4+W1+dx
B5 = B1
W5 = W1
H5 = H1

# eta_ss(Z) plot:
L6 = L5+W1+dx
B6 = B1
W6 = W1
H6 = H1

# ======================================
# color
from matplotlib import cm
cmap1 = cm.Greens
cm1_range = np.linspace(0.5,0.99,f_band.size)
cmap2 = cm.Blues
cm2_range = np.linspace(0.5,0.99,f_band.size)

# MAKE THE x_limits adaptive !!



# ======================================
# M1,M2(Z) PLOT
ax = plt.axes([L1,B1,W1,H1])

for f in range(f_band.size):
    #plotting!
    J1_interp_AndPsP = np.interp(Z_km_interp,Z_km,J1_f_AndPsP_mat[:,f])
    ax.plot(J1_interp_AndPsP,Z_km_interp,lw=1.0, ls='-', alpha=0.8, color=cmap1(cm1_range[f]))
    ax.plot(J1_interp_AndPsP[ind_LAB_AndPsP_fmat[f]],Z_km_interp[ind_LAB_AndPsP_fmat[f]],'r.')

    J1_interp_YTmax = np.interp(Z_km_interp,Z_km,J1_f_YTmax_mat[:,f])
    ax.plot(J1_interp_YTmax,Z_km_interp,lw=1.0, ls='-', alpha=0.8,color=cmap2(cm2_range[f]) )
    ax.plot(J1_interp_YTmax[ind_LAB_YTmax_fmat[f]],Z_km_interp[ind_LAB_YTmax_fmat[f]],'k.')

xmin = 1.3e-11
xmax = 1.45e-11
ax.set_xlim([xmin,xmax]) # used to be axes

ymin = 40.0
ymax = 100 # 350.0
ax.set_ylim([ymin,ymax]) # used to be axes'

ax.set_title('Storage compliance')
ax.set_xlabel('$J_1$ [Pa$^{-1}$]')
#ax.set_yticklabels([])
ax.set_ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax.invert_yaxis()

# =======================
# PLOT THE J2 STRUCTURES

ax2 = plt.axes([L2,B2,W2,H2])

for f in range(f_band.size):
    #plotting!
    J2_interp_AndPsP = np.interp(Z_km_interp,Z_km,J2_f_AndPsP_mat[:,f])
    ax2.plot(np.log10(J2_interp_AndPsP),Z_km_interp,lw=1.0, ls='-', alpha=0.8, color=cmap1(cm1_range[f]) )
    ax2.plot(np.log10(J2_interp_AndPsP[ind_LAB_AndPsP_fmat[f]]),Z_km_interp[ind_LAB_AndPsP_fmat[f]],'r.')

    J2_interp_YTmax = np.interp(Z_km_interp,Z_km,J2_f_YTmax_mat[:,f])
    ax2.plot(np.log10(J2_interp_YTmax),Z_km_interp,lw=1.0, ls='-', alpha=0.8, color=cmap2(cm2_range[f]) )
    ax2.plot(np.log10(J2_interp_YTmax[ind_LAB_YTmax_fmat[f]]),Z_km_interp[ind_LAB_YTmax_fmat[f]],'k.')

xmin = -20.0
xmax = -10.0
ax2.set_xlim([xmin,xmax]) # used to be axes

ymin = 40.0
ymax = 100.0
ax2.set_ylim([ymin,ymax]) # used to be axes'

ax2.set_title('Loss compliance')
ax2.set_xlabel('log $J_2$ [Pa$^{-1}$]')
ax2.set_yticklabels([])
#ax2.set_ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax2.invert_yaxis()


# =======================
# PLOT THE M(f) STRUCTURES

ax3 = plt.axes([L3,B3,W3,H3])

for f in range(f_band.size):
    #plotting!
    Mf_interp_AndPsP = np.interp(Z_km_interp,Z_km,Mf_AndPsP_mat[:,f])
    ax3.plot(Mf_interp_AndPsP/1e9,Z_km_interp,lw=1.0, ls='-',alpha=0.8,color=cmap1(cm1_range[f]) )
    ax3.plot(Mf_interp_AndPsP[ind_LAB_AndPsP_fmat[f]]/1e9,Z_km_interp[ind_LAB_AndPsP_fmat[f]],'r.')

    Mf_interp_YTmax = np.interp(Z_km_interp,Z_km,Mf_YTmax_mat[:,f])
    ax3.plot(Mf_interp_YTmax/1e9,Z_km_interp,lw=1.0, ls='-',alpha=0.8,color=cmap2(cm2_range[f]) )
    ax3.plot(Mf_interp_YTmax[ind_LAB_YTmax_fmat[f]]/1e9,Z_km_interp[ind_LAB_YTmax_fmat[f]],'k.')

xmin = 68
xmax = 76
ax3.set_xlim([xmin,xmax]) # used to be axes

ax3.set_ylim([ymin,ymax]) #

ax3.set_title('$M(f)|_z$')
ax3.set_xlabel('$M(f)$ [GPa]')
ax3.set_yticklabels([])
#ax2.set_ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax3.invert_yaxis()



# =======================
# PLOT THE Q STRUCTURES

ax4 = plt.axes([L4,B4,W4,H4])

for f in range(f_band.size):
    #plotting!
    Qf_interp_AndPsP = np.interp(Z_km_interp,Z_km,Qf_AndPsP_mat[:,f])
    ax4.plot(np.log10(Qf_interp_AndPsP),Z_km_interp,lw=1.0,ls='-',alpha=0.8,color=cmap1(cm1_range[f]))
    ax4.plot(np.log10(Qf_interp_AndPsP[ind_LAB_AndPsP_fmat[f]]),Z_km_interp[ind_LAB_AndPsP_fmat[f]],'r.')

    Qf_interp_YTmax = np.interp(Z_km_interp,Z_km,Qf_YTmax_mat[:,f])
    ax4.plot(np.log10(Qf_interp_YTmax),Z_km_interp,lw=1.0, ls='-',alpha=0.8,color=cmap2(cm2_range[f]))
    ax4.plot(np.log10(Qf_interp_YTmax[ind_LAB_YTmax_fmat[f]]),Z_km_interp[ind_LAB_YTmax_fmat[f]],'k.')


xmin = 1.0
xmax = 9.0
ax4.set_xlim([xmin,xmax]) # used to be axes


ax4.set_ylim([ymin,ymax]) # used to be axes'

ax4.set_title('Q [find LAB]')
ax4.set_xlabel('log $Q$')
ax4.set_yticklabels([])
#ax2.set_ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax4.invert_yaxis()


# =======================
# PLOT THE Q STRUCTURES

ax5 = plt.axes([L5,B5,W5,H5])

for f in range(f_band.size):
    #plotting!
    Vf_interp_AndPsP = np.interp(Z_km_interp,Z_km,Vf_AndPsP_mat[:,f])
    ax5.plot(Vf_interp_AndPsP/1000,Z_km_interp,lw=1.0, ls='-',alpha=0.8,color=cmap1(cm1_range[f]) )
    ax5.plot(Vf_interp_AndPsP[ind_LAB_AndPsP_fmat[f]]/1000,Z_km_interp[ind_LAB_AndPsP_fmat[f]],'r.')

    Vf_interp_YTmax = np.interp(Z_km_interp,Z_km,Vf_YTmax_mat[:,f])
    ax5.plot(Vf_interp_YTmax/1000,Z_km_interp,lw=1.0, ls='-', alpha=0.8,color=cmap2(cm2_range[f]))
    ax5.plot(Vf_interp_YTmax[ind_LAB_YTmax_fmat[f]]/1000,Z_km_interp[ind_LAB_YTmax_fmat[f]],'k.')


xmin = 4.5
xmax = 4.9
ax5.set_xlim([xmin,xmax]) # used to be axes


ax5.set_ylim([ymin,ymax]) # used to be axes'

ax5.set_title('$V_s$')
ax5.set_xlabel('$V_s$ [km/sec]')
ax5.set_yticklabels([])
#ax2.set_ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax5.invert_yaxis()


plt.show()

# ===========
#plt.show()
#sys.exit()
# ===========
