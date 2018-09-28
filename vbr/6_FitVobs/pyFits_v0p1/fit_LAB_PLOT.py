# ===================================================================
#
# ===================================================================
# ======================
# PLOT THE CURVES !
# ======================
import numpy as np
import matplotlib.pyplot as plt

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
Z_LAB_Q_km_COLD = Z_LAB_Q_mat_COLD[ij_best_COLD]
zLAB_Q_line_COLD = Z_LAB_Q_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_Q_line_COLD,color='blue',lw=2.0, linestyle=':', alpha=0.5)
# add HOT plate Q LAB line:
Z_LAB_Q_km_HOT = Z_LAB_Q_mat_HOT[ij_best_HOT]
zLAB_Q_line_HOT = Z_LAB_Q_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_Q_line_HOT,color='red',lw=2.0, linestyle=':', alpha=0.5)

# # add COLD plate LT viscosity LAB line:
# Z_LAB_LT_km_COLD = Z_km[int(ind_zLAB_LT_mat_COLD[ij_best_COLD])]
# zLAB_LT_line_COLD = Z_LAB_LT_km_COLD*np.ones(len(Z_km))
# ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_COLD,color='blue',lw=2.0, linestyle='--', alpha=0.5)
# # add HOT plate LT viscosity LAB line:
# Z_LAB_LT_km_HOT = Z_km[int(ind_zLAB_LT_mat_COLD[ij_best_HOT])]
# zLAB_LT_line_HOT = Z_LAB_LT_km_HOT*np.ones(len(Z_km))
# ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_LT_line_HOT,color='red',lw=2.0, linestyle='--', alpha=0.5)


# add (OBSERVED!) COLD plate thickness line :
zLAB_line_COLD = zLAB_obs_km_COLD*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_COLD,color='blue',lw=7.0, linestyle='-', alpha=0.2)
# add HOT plate thickness line:
zLAB_line_HOT = zLAB_obs_km_HOT*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line_HOT,color='red',lw=7.0, linestyle='-', alpha=0.2)

plt.title('Fit the LAB with Qs')
plt.xlabel('log Qs')
plt.ylabel('Depth [km]')
#plt.autoscale(enable=True, axis='y', tight=True)
ax.invert_yaxis()


# =================================
# plot the residuals for Q_LAB HOT
L2 = 0.4
B2 = B1+0.4
W2 = 0.3
H2 = W2

ax2 = plt.axes([L2,B2,W2,H2])
xy_rangelist = [zPlate_vec[0],zPlate_vec[-1],Tpot_vec[-1],Tpot_vec[0]] # flip vertical directions?
ax2.imshow(np.log10(Res_lab_Q_mat_HOT), cmap=plt.cm.gray, extent=xy_rangelist) #, aspect="auto")
#ax2.set_aspect(1) #
ax2.scatter(zPlate_vec[ij_best_HOT[1]], Tpot_vec[ij_best_HOT[0]], color='red', s=10)

ax2.set_title('Res: $Z_{LAB}$ (Hot plate)')
#ax2.set_xticklabels(zPlate_vec)
#ax2.xaxis.set_major_locator(plt.LinearLocator(5))
ax2.set_xlabel('z_plate')
#ax2.set_yticklabels(Tpot_vec)
#ax2.yaxis.set_major_locator(plt.LinearLocator(5))
ax2.set_ylabel('$T_{pot}$')

# =================================
# plot the residuals for Q_LAB COLD
L3 = L2
B3 = B1
W3 = W2
H3 = H2

ax3 = plt.axes([L3,B3,W3,H3])
ax3.imshow(np.log10(Res_lab_Q_mat_COLD), cmap=plt.cm.gray, extent=xy_rangelist)
ax3.scatter(zPlate_vec[ij_best_COLD[1]], Tpot_vec[ij_best_COLD[0]], color='blue', s=10)

#ax3.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax3.set_title('Res: Z_LAB (Cold plate)')
#ax3.set_xticklabels(zPlate_vec)
ax3.set_xlabel('z_plate')
#ax3.set_yticklabels(Tpot_vec)
ax3.set_ylabel('T_pot')

# =================================

# =================================
# plot the residuals for Vs_adavg HOT
L2 = 0.7
B2 = B1+0.4
W2 = 0.3
H2 = W2

ax4 = plt.axes([L2,B2,W2,H2])
ax4.imshow(np.log10(Res_Vs_adavg_mat_HOT), cmap=plt.cm.gray, extent=xy_rangelist)
ax4.scatter(zPlate_vec[ij_best_HOT[1]], Tpot_vec[ij_best_HOT[0]], color='red', s=10)

ax4.set_title('Res: Vs (Hot plate)')
#ax4.set_xticklabels(zPlate_vec)
ax4.set_xlabel('z_plate')
#ax4.set_yticklabels(Tpot_vec)
ax4.set_ylabel('T_pot')

# =================================
# plot the residuals for Vs_adavg COLD
L3 = L2
B3 = B1
W3 = W2
H3 = H2

ax5 = plt.axes([L3,B3,W3,H3])
ax5.imshow(np.log10(Res_Vs_adavg_mat_COLD), cmap=plt.cm.gray, extent=xy_rangelist)
ax5.scatter(zPlate_vec[ij_best_COLD[1]], Tpot_vec[ij_best_COLD[0]], color='blue', s=10)

#ax5.invert_yaxis()
#orientation confused because of imshow and THEN scatter?
# this doesnt work: ax3.colorbar()
ax5.set_title('Res: Vs (Cold plate)')
#ax5.set_xticklabels(zPlate_vec)
ax5.set_xlabel('z_plate')
#ax3.set_yticklabels(Tpot_vec)
ax5.set_ylabel('T_pot')

plt.show()

#import GIA_transient_plates_mat_v2

# =============
#sys.exit()
# =============

# ============
# sys.exit()
# ============
