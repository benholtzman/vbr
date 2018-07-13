# PLOT THE GIA predictions...

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import sys
from time import time
from imp import reload  # this lets us reload a module without restarting the kernel ! very convenient now...
import pandas as pd

# our functions:
import fit_LAB_fcns as flab
reload(flab)

# ============================================================
path = '/Users/ben/0_vbr_git/VBRcloset/2018-07-02-BigBox/'
matobj = 'Box_2018-07-02-BigBox_VBR_GIA_py.mat'
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)


zLAB_obs_km_COLD = 180 # plot this line as the seismically identified plate thickness..
zLAB_obs_km_HOT = 80 # plot this line as the seismically identified plate thickness..

Z_find_M_COLD = 0.8*zLAB_obs_km_COLD # find the values of M at this depth--
Z_find_M_HOT = 0.8*zLAB_obs_km_HOT # find the values of M at this depth--
# this is the line at which to find the values of M
# BUT I SHOULD PROLLY KILL THIS !!

# -----------------
# read in the Box

t0 = time() # see how long it takes !
hoozit = sio.whosmat(path+matobj) # get a little info without loading it..
print(hoozit)
# load the mat file:
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)
print(type(b))
print(b.keys())
# The box itself is a structured object (a class operating on numpy array?)
# access field names by [name]._fieldnames
box = b['Box']
#print('box is a ' + str(type(box)) + ', with shape:')
#print(box.shape)
n_var1 = box.shape[0]
n_var2 = box.shape[1]
#print(box[0,2]._fieldnames)
#print(box[0,2].info._fieldnames)
#print(box[0,2].run_info._fieldnames)

zPlate_vec = box[0,0].info.var2range # this flipped from var1range
Tpot_vec = box[0,0].info.var1range
print('Tpot_vec = ', str(Tpot_vec) )
print('zPlate_vec = ', str(zPlate_vec) )

# testing the box:
# frame_test = box[0,0].Frames[-1]
# print(frame._fieldnames)
# print(frame.VBR._fieldnames)
# print('VBR input:')
# print(frame.VBR.input._fieldnames ) # shit !! .in does not work ! already a built in function !
# print('VBR out:')
# print(frame.VBR.out._fieldnames)
# print("That took " + str(time() - t0) + " s")
# f_band = frame_test.VBR.input.SV.f
# print('full frequency band:')
# print(f_band)


# ====================================
# Pull the best fits for HOT and COLD plates !
# (after restructuring, this may not be needed.. stays in memory in %run in ipython)

LAB_bestfits = pd.read_pickle('LAB_bestfits.pkl')
print(LAB_bestfits.HOTplate)
print(LAB_bestfits.COLDplate)

# later, make a class that does all this? hmmmm.... based on the animation codes!

# ================================================
# GET THE DATA OUT:
# ================================================

# COLD:
ij_best = [LAB_bestfits.COLDplate[0], LAB_bestfits.COLDplate[1]]# from the previous code !
frame = box[ij_best[0],ij_best[1]].Frames[-1]
Z_km_COLD = box[ij_best[0],ij_best[1]].run_info.Z_km
# frequency band:
freq_vec = f_band = frame.VBR.input.SV.f
time_vec = 1/freq_vec[::-1] # reverse direction, short time to long time.
num_freqs = len(freq_vec)
# print(f_band)

# # extract the Frequency-INDEPENDENT properties we need:
# eta = frame.VBR.out.viscous.LH2012.eta_total # composite viscosity
# eta_diff = frame.VBR.out.viscous.LH2012.diff.eta # viscosity for diff creep
# Gu = frame.VBR.out.elastic.anharmonic.Gu # unrelaxed shear modulus
# TauMxw = eta/Gu
# TauMxw_diff = eta_diff/Gu

# extract the Frequency-DEPENDENT properties we need over f band of interest:
Mod_f_mat_COLD = frame.VBR.out.anelastic.AndradePsP.Ma
dims = Mod_f_mat_COLD.shape
print('cold_M dims'+str(dims))


# HOT:
ij_best = [LAB_bestfits.HOTplate[0], LAB_bestfits.HOTplate[1]]# from the previous code !
frame = box[ij_best[0],ij_best[1]].Frames[-1]
Z_km_HOT = box[ij_best[0],ij_best[1]].run_info.Z_km
# frequency band:
#freq_vec = f_band = frame.VBR.input.SV.f
#time_vec = 1/freq_vec[::-1] # reverse direction, short time to long time.
#num_freqs = len(freq_vec)
# print(f_band)

# # extract the Frequency-INDEPENDENT properties we need:
# eta = frame.VBR.out.viscous.LH2012.eta_total # composite viscosity
# eta_diff = frame.VBR.out.viscous.LH2012.diff.eta # viscosity for diff creep
# Gu = frame.VBR.out.elastic.anharmonic.Gu # unrelaxed shear modulus
# TauMxw = eta/Gu
# TauMxw_diff = eta_diff/Gu

# extract the Frequency-DEPENDENT properties we need over f band of interest:
Mod_f_mat_HOT = frame.VBR.out.anelastic.AndradePsP.Ma
dims = Mod_f_mat_HOT.shape
print('cold_M dims'+str(dims))

# ================================================
#  PICK A THRESHOLD M value # TO DEFINE APPARENT PLATE THICKNESS !
# make sure it coincides at lowest freq with the eta_LAB = 10*eta_asth_avg.
Mod_LAB = 20.0  # THIS NEEDS TO BE SELF CONSISTENT !

#np.savetxt(pathtodata+'zplate_freq.txt', zplate_freq)
Z_find_M = Z_find_M_COLD
zLAB_obs_km = zLAB_obs_km_COLD
zplate_freq, zplate_time_COLD, i_zM, i_zLAB = flab.find_zplate_GIA_freq(Z_km_COLD,Z_find_M,zLAB_obs_km,Mod_f_mat_COLD,Mod_LAB,freq_vec)
#print(zplate_freq)

Z_find_M = Z_find_M_HOT
zLAB_obs_km = zLAB_obs_km_HOT
zplate_freq, zplate_time_HOT, i_zM, i_zLAB = flab.find_zplate_GIA_freq(Z_km_HOT,Z_find_M,zLAB_obs_km,Mod_f_mat_HOT,Mod_LAB,freq_vec)
#print(zplate_freq)


# ===================================================================
# ===================================================================
# >>>>>>>>>>>>>>>>>>
# PLOT THE CURVES !
# >>>>>>>>>>>>>>>>>>
plt.figure(1)
# plasma colormap
cmap = plt.cm.get_cmap('autumn')

# # =======================
# #f, axs = plt.subplots(1,2,figsize=(7,7))
# #plt.subplot(1, 2, 1)
L1 = 0.1
B1 = 0.2
W1 = 0.25
H1 = 0.6
#
# ax = plt.axes([L,B,W,H])
#
#
# bk = [0.0,0.0,0.0,1.0]
#
# for i_freq,freq in enumerate(freq_vec):
#     color_num = i_freq/float(num_freqs)
#     #color_num = freq / max(freq_vec)
#     col = cmap(color_num)
#     time_line = np.log10(time_vec[-i_freq]/(np.pi*1e7)*np.ones(len(Z_km)))
#     #plt.plot(time_line,Z_km,color=col,lw=2.0 )
#     ax.plot(time_line,Z_km,color=col,lw=2.0 )
#
# time_line = np.log10(time_vec[0]/(np.pi*1e7)*np.ones(len(Z_km)))
# ax.plot(time_line,Z_km,color=bk,lw=1.0 )
# time_line = np.log10(time_vec[-1]/(np.pi*1e7)*np.ones(len(Z_km)))
# ax.plot(time_line,Z_km,color=bk,lw=1.0 )
#
#
# ax.plot(np.log10(TauMxw_diff/(np.pi*1e7)),Z_km,color='green',lw=2.0, ls='--' )
# ax.plot(np.log10(TauMxw/(np.pi*1e7)),Z_km,color='blue',lw=3.0 )
#
# axes = plt.gca()
# xmin = 1.0
# xmax = 8.0
# axes.set_xlim([xmin,xmax])
#
# # add plate thickness line:
# zLAB_obs_km = zLAB_obs_km_COLD
# zLAB_line = zLAB_obs_km*np.ones(len(Z_km))
# ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line,color='black',lw=5.0, alpha=0.4)
#
# plt.gca().invert_yaxis()
# plt.title('Maxwell time')
# plt.xlabel('log10 time [yrs]')
# plt.ylabel('Depth [km]')
# plt.autoscale(enable=True, axis='y', tight=True)

# =======================
#plt.subplot(1, 2, 2)
#plt.figure(figsize=(5,10))
L2 = L1
B2 = B1
W2 = W1
H2 = H1
ax2 = plt.axes([L2,B2,W2,H2])

cmap_HOT = 'Oranges'
cmap_COLD = 'Blues'



#for i_freq in range(num_freqs):
# scrolling through BACKWARDS (short to long time ! )
#for i_freq in range(num_freqs-1,0,-1):
for i_freq,freq in enumerate(freq_vec):
    #print(i_freq)
    #color_num = freq / max(freq_vec)
    color_num = i_freq/float(num_freqs)
    cmap=plt.get_cmap(cmap_COLD)
    col = cmap(color_num)
    m_f = (Mod_f_mat_COLD[:,i_freq])/1e9
    ax2.plot(m_f,Z_km_COLD,color=col,lw=2.0 )

for i_freq,freq in enumerate(freq_vec):
    color_num = i_freq/float(num_freqs)
    cmap=plt.get_cmap(cmap_HOT)
    col = cmap(color_num)
    m_f = (Mod_f_mat_HOT[:,i_freq])/1e9
    ax2.plot(m_f,Z_km_HOT,color=col,lw=2.0 )

m_f = (Mod_f_mat_HOT[:,0])/1e9
ax2.plot(m_f,Z_km_HOT,color='black',lw=1.0 )
m_f = (Mod_f_mat_HOT[:,num_freqs-1])/1e9
ax2.plot(m_f,Z_km_HOT,color='black',lw=1.0 )

m_f = (Mod_f_mat_COLD[:,0])/1e9
ax2.plot(m_f,Z_km_COLD,color='black',lw=1.0 )
m_f = (Mod_f_mat_COLD[:,num_freqs-1])/1e9
ax2.plot(m_f,Z_km_COLD,color='black',lw=1.0 )

xmin = 0.0
xmax = 80.0

zLAB_line_HOT = zLAB_obs_km_HOT*np.ones(len(Z_km_HOT))
ax2.plot(np.linspace(xmin,xmax,len(Z_km_HOT)),zLAB_line_HOT,color='red',lw=6.0, alpha=0.4)
zLAB_line_COLD = zLAB_obs_km_COLD*np.ones(len(Z_km_COLD))
ax2.plot(np.linspace(xmin,xmax,len(Z_km_COLD)),zLAB_line_COLD,color='blue',lw=6.0, alpha=0.4)
#
# ModLAB_line = Mod_LAB*np.ones(len(Z_km[0:i_zLAB]))
# ax2.plot(ModLAB_line,Z_km[0:i_zLAB],color='purple',lw=10.0, alpha=0.2)
#
# # add M @ depth line:
# z_Mt_line = Z_find_M *np.ones(len(Z_km))
# ax2.plot(np.linspace(xmin,xmax,len(Z_km)),z_Mt_line,color='black',lw=1.0, alpha=0.7)

axes = plt.gca()
axes.set_xlim([xmin,xmax])
axes.invert_yaxis()
axes.set_yticklabels([])

plt.title('Modulus(f)')
plt.xlabel('Modulus [GPa]')
plt.ylabel('')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)




# ================================================
# PLOT THE EXTRACTED DATA !
# ================================================
#f0 = plt.subplots(2,1,figsize=(6,8))
#up = 0.36
L2 = L1 + 0.33
#B2 = B + up
W2 = W1
H2 = 0.4*H1
#ax3 = plt.axes([L,B,W,H])


# >>>>>>>>>>>>>>>>>>
# PLOT THE modulus change with frequency !
## >>>>>>>>>>>>>>>>>>
#plt.subplot(2,1,1)

# ax3.plot(freq_vec,(Mod_f_mat[i_zM,:])/1e9, c='black')
#
# for i_freq,freq in enumerate(freq_vec):
#     color_num = i_freq/float(num_freqs)
#     #color_num = freq / max(freq_vec)
#     col = cmap(color_num)
#     plt.plot(freq_vec[i_freq],(Mod_f_mat[i_zM,i_freq])/1e9,c=col,
#                        marker='o',
#                        markersize=10, markeredgecolor = 'black' )
#
# plt.title('M(f) at ' + str(Z_find_M) + ' km')
# plt.xlabel('frequency')
# plt.ylabel('M [GPa]')


# >>>>>>>>>>>>>>>>>>
# PLOT THE LAB EVOLUTION !
# >>>>>>>>>>>>>>>>>>
#plt.subplot(2,1,2)

L3 = L2
B3 = B1 + 0.25*H1
W3 = W1
H3 = H2
ax4 = plt.axes([L3,B3,W3,H3])

sec_yr = np.pi*1e7
time_vec = zplate_time_COLD[:,0]/sec_yr
zlab_ara_t = zplate_time_COLD[:,2]
n_time = len(time_vec)
ax4.plot(time_vec,zlab_ara_t, c='black')
ax4.plot([time_vec[0], time_vec[-1]],[zLAB_obs_km_COLD,zLAB_obs_km_COLD], color='blue',lw=6.0, alpha=0.4)

time_vec = zplate_time_HOT[:,0]/sec_yr
zlab_ara_t = zplate_time_HOT[:,2]
n_time = len(time_vec)
ax4.plot(time_vec,zlab_ara_t, c='black')
ax4.plot([time_vec[0], time_vec[-1]],[zLAB_obs_km_HOT,zLAB_obs_km_HOT], color='red',lw=6.0, alpha=0.4)

#for i_freq in range(num_freqs-1):
#    color_num = i_freq/float(num_freqs)
#for i_freq,freq in enumerate(freq_vec):
for i_time,time in enumerate(time_vec):
    color_num = i_time/float(n_time)
    #color_num = freq / max(freq_vec)
    cmap=plt.get_cmap(cmap_COLD)
    col = cmap(color_num)
    ax4.plot(time,zplate_time_COLD[i_time,2],c=col,
                       marker='o',
                       markersize=10, markeredgecolor = 'black' )

for i_time,time in enumerate(time_vec):
    color_num = i_time/float(n_time)
    #color_num = freq / max(freq_vec)
    cmap=plt.get_cmap(cmap_HOT)
    col = cmap(color_num)
    ax4.plot(time,zplate_time_HOT[i_time,2],c=col,
                       marker='o',
                       markersize=10, markeredgecolor = 'black' )

#axes = plt.gca()
#ax4.set_ylim([0.9*min(zlab_ara_t),1.1*max(zlab_ara_t)])
ax4.invert_yaxis()

ax4.set_title('apparent plate thickness')
ax4.set_xlabel('time [yrs]')
ax4.set_ylabel('thickness [km]')


#
plt.show()
# # ============
# plt.show()
# sys.exit()
# # ============



# ===================================================================
# ===================================================================
# PLOT and write out plate thickness vs frequency:
# for the laplace transform to time domain:

# plt.figure(0)
# plt.plot(freq_vec,zlab_ara)
# plt.plot([freq_vec[0], freq_vec[-1]],[zLAB_obs_km,zLAB_obs_km], lw=10, color='gray', alpha=0.5)
# plt.xlabel('Frequency (or strain rate) [Hz]')
# plt.ylabel('Apparent plate thickness [km]')
# plt.title('Apparent plate thickness with frequency')
# plt.autoscale(enable=True, axis='y', tight=True)
# axes = plt.gca()
# ymin = zlab_ara[0]
# ymax = zLAB_obs_km+0.05*zLAB_obs_km
# axes.set_ylim([ymin,ymax])
