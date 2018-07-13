# functions for fitting LAB depth !
import numpy as np

# function to get mean and standard deviation at each depth:  (std = sqrt(mean(abs(x - x.mean())**2)))
# Pp_f_mnstd = zip_meanstd_fband(Pprop_fband)
def zip_meanstd_fband(Pprop_fband):
    dims = Pprop_fband.shape
    Pp_f_mnstd = np.zeros([dims[0],2])
    for i in range(dims[0]):
        Pp_f_mnstd[i,0] = np.mean(Pprop_fband[i,:])
        Pp_f_mnstd[i,1] = np.std(Pprop_fband[i,:])
    return Pp_f_mnstd

# ==================
# function to find best fitting models (calculate residuals):
# Res_lab_Q_mat, ind_zLAB_Q_mat = flab.find_LAB_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax)
def find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax):

    # get the size of the box
    dims = box.shape
    n_var1 = dims[0]
    n_var2 = dims[1]
# ALSO MAKE AN ARRAY FOR THE INDEX OF LAB IN EACH
    Res_lab_Q_mat = np.zeros([n_var1,n_var2])
    ind_zLAB_Q_mat = np.zeros([n_var1,n_var2])
    #Res_lab_etaLT_mat = np.zeros([n_var1,n_var2])

    for i_var1 in range(n_var1):
        for i_var2 in range(n_var2):

            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]

            Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
            Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
            Qs_mnstd = zip_meanstd_fband(Qs_f_band)
            # INTERPOLATE TO HIGHER RES ?

            Qlab_near = Qs_mnstd[:,0][Qs_mnstd[:,0]<=Q_LAB][0] # going from high values to low wi increasing Z
            ind_Qlab = int(np.where(Qs_mnstd[:,0]==Qlab_near)[0])

            Z_LAB_Qs_km = Z_km[ind_Qlab]

    #        Z_LAB_LT_km = box[i_var1,i_var2].run_info.zLAB[-1]/1e3

            # calculate the residual between the predicted and the measured:
            Res = ((Z_LAB_Qs_km - zLAB_obs_km)**2)/zLAB_obs_km
            Res_lab_Q_mat[i_var1,i_var2] = Res # np.log10(Res)
            ind_zLAB_Q_mat[i_var1,i_var2] = ind_Qlab
    #        Res = ((Z_LAB_LT_km- zLAB_obs_km)**2)/zLAB_obs_km
    #        Res_lab_etaLT_mat[i_var1,i_var2] = Res # np.log10(Res)

    return Res_lab_Q_mat, ind_zLAB_Q_mat


# ==================
# function to find the LAB defined by steady state viscosity (labLT):
# Res_lab_etaLT_mat, ind_zLAB_LT_mat = flab.find_LAB_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax)

def find_LAB_LT_Res(box,eta_LAB,zLAB_obs_km):

    # get the size of the box
    dims = box.shape
    n_var1 = dims[0]
    n_var2 = dims[1]
    Res_lab_etaLT_mat = np.zeros([n_var1,n_var2])
    ind_zLAB_LT_mat = np.zeros([n_var1,n_var2])

    for i_var1 in range(n_var1):
        for i_var2 in range(n_var2):

            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]
            eta = frame.VBR.out.viscous.LH2012.eta_total # composite viscosity

            #fmax = f_band[f_band<=f_bw_max][-1]
            #i_fmax = int(np.where(f_band==fmax)[0])

            # how did these switch ?
            # z_plate = float(box[i_var1,i_var2].info.var1val)
            z_plate = float(box[i_var1,i_var2].info.var2val)
            #print('z_plate = ', str(z_plate))
            ind_zplate = int(np.where(Z_km >= z_plate)[0][0])
            eta_adbt_vec = eta[ind_zplate:-1]
            eta_adbt_mean = np.mean(eta_adbt_vec)
            # -------------------------
            eta_LAB = 1000*eta_adbt_mean
            #eta_LAB = eta_LAB
            # -------------------------
            ind_zLAB_LT = int(np.where(eta <= eta_LAB)[0][0])
            zLAB_LT = Z_km[ind_zLAB_LT]

            # find the residual:
            Res = ((zLAB_LT- zLAB_obs_km)**2)/zLAB_obs_km
            Res_lab_etaLT_mat[i_var1,i_var2] = Res # np.log10(Res)
            ind_zLAB_LT_mat[i_var1,i_var2] = ind_zLAB_LT

    return Res_lab_etaLT_mat, ind_zLAB_LT_mat

# ================================================
# GIA TIMESCALE !
# ================================================

# get the nearest value in an array to val:
#(is this duplicated in purpose but different from somewhere else?)
# yes : int(np.where(eta <= eta_LAB)[0][0] -- compare results !
def getnearest(array,val):
    idx = (np.abs(array-val)).argmin()
    return idx

# =========================
# function to find the plate thickness through the GIA timescale sweep:
# zplate_freq, zplate_time, i_zM, i_zLAB = flab.find_zplate_GIA_freq(Z_km,Z_find_M,zLAB_obs_km,Mod_f_mat,Mod_LAB)

def find_zplate_GIA_freq(Z_km,Z_find_M,zLAB_obs_km,Mod_f_mat,Mod_LAB,freq_vec):
    n_freqs = Mod_f_mat.shape[1]
    ilab_ara = np.zeros(len(freq_vec))
    zlab_ara = np.zeros(len(freq_vec))

    # find the index at which Z=Z where you find the M(time) curve.
    i_zM = getnearest(Z_km,Z_find_M)
    # find the index at which Z=Z_LAB.
    i_zLAB = getnearest(Z_km,zLAB_obs_km)

    for i_freq in range(n_freqs):
        m_f = (Mod_f_mat[:,i_freq])/1e9  # remember, these go from low freq to high freq !
        i_lab = getnearest(m_f[0:i_zLAB+20],Mod_LAB)  # 3/26/18, extended the search below i_zLAB..
        # could be source of future bugs too!

        ilab_ara[i_freq] = i_lab
        zlab_ara[i_freq] = Z_km[i_lab]

    # make a little array of zPlate as a function of frequency
    zplate_freq = np.zeros((n_freqs,2))
    zplate_freq[:,0] = freq_vec
    zplate_freq[:,1] = zlab_ara
    # flip these now so they correspond to increasing time...
    zplate_time = np.zeros((n_freqs,3))
    zplate_time[:,0] = time_vec = 1/(freq_vec[::-1])
    zplate_time[:,1] = ilab_ara[::-1]
    zplate_time[:,2] = zlab_ara[::-1]


    return zplate_freq, zplate_time, i_zM, i_zLAB
