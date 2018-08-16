# functions for fitting LAB depth !
import numpy as np

# To Do:

# (1) the anelastic model is hard wired in here..
# better to let that choice be an option, in some way.

# (2) option for LAB finding method,
# between normalized and absolute values-- need both,
# at least in find_LAB_Q_one , for GIA band.
# to pass into find_zplate_GIA_freq

# function to get mean and standard deviation at each depth:
# (std = sqrt(mean(abs(x - x.mean())**2)))
# Pp_f_mnstd = zip_meanstd_fband(Pprop_fband)
def zip_meanstd_fband(Pprop_fband):
    dims = Pprop_fband.shape
    Pp_f_mnstd = np.zeros([dims[0],2])
    for i in range(dims[0]):
        Pp_f_mnstd[i,0] = np.mean(Pprop_fband[i,:])
        Pp_f_mnstd[i,1] = np.std(Pprop_fband[i,:])
    return Pp_f_mnstd

# get the nearest value in an array to val:
#(is this duplicated in purpose but different from somewhere else?)
# yes : int(np.where(eta <= eta_LAB)[0][0] -- compare results !
def getnearest(array,val):
    idx = (np.abs(array-val)).argmin()
    # OR
    #idx = int(np.where(array >= val)[0][0])
    return idx


def find_LAB_Q_one(Qz,Z_km,z_plate,fQ_LAB,zLAB_obs_km):
    # find the average Q in the adiabatic part:
    ind_zplate =  getnearest(Z_km,z_plate)
    Q_adbt_vec = Qz[ind_zplate:-1]
    Q_adbt_mean = np.mean(Q_adbt_vec)
    # calculate Q at the LAB:
    Q_LAB = fQ_LAB*Q_adbt_mean
    Qlab_near = Qz[Qz<=Q_LAB][0]
    ind_Qlab = int(getnearest(Qz,Qlab_near)) # int(np.where(Qs_mnstd[:,0]==Qlab_near)[0])

    Z_LAB_Qs_km = Z_km[ind_Qlab]

    return Z_LAB_Qs_km, ind_Qlab

# ==================
# function to find best fitting models (calculate residuals):
# Res_lab_Q_mat, ind_zLAB_Q_mat = flab.find_LAB_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax)
def find_LAB_Q_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax):

    # get the size of the box
    dims = box.shape
    n_var1 = dims[0]
    n_var2 = dims[1]
    # MAKE ARRAYS for the residual, the depth and the index of LAB depth in each model
    Res_lab_Q_mat = np.zeros([n_var1,n_var2])
    ind_zLAB_Q_mat = np.zeros([n_var1,n_var2])
    Z_LAB_Q_mat = np.zeros([n_var1,n_var2])

    for i_var1 in range(n_var1):
        for i_var2 in range(n_var2):

            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]

            # BH: TAKE THIS STEP OUT-- put in its own method..
            # only pass in to this method the 1D values...
            #Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
            Qs_f_mat = frame.VBR.out.anelastic.YT_maxwell.Q
            Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
            Qs_mnstd = zip_meanstd_fband(Qs_f_band)
            # INTERPOLATE Qz and Z_km TO HIGHER RES:
            nn_pts = 1000
            Z_km_interp = np.linspace(Z_km[0],Z_km[-1],nn_pts)
            Qz_interp = np.interp(Z_km_interp,Z_km,Qs_mnstd[:,0])

            # IDENTIFY THE LAB depth in each model:
            # method 0: an absolute value of Q (use the Q_LAB read in)
            # method 1: a factor above the adiabatic average
            method=1
            if method==1:
                fQlab = 20
                # find the average Q in the adiabatic part:
                z_plate = float(box[i_var1,i_var2].info.var2val)
                ind_zplate =  getnearest(Z_km_interp,z_plate)
                Q_adbt_vec = Qz_interp[ind_zplate:-1]
                Q_adbt_mean = np.mean(Q_adbt_vec)
                # calculate Q at the LAB:
                Q_LAB = fQlab*Q_adbt_mean
            # find the position of the nearest Q value to Q_LAB
            # Qlab_near = Qs_mnstd[:,0][Qs_mnstd[:,0]<=Q_LAB][0] # going from high values to low wi increasing Z
            Qlab_near = Qz_interp[Qz_interp<=Q_LAB][0]
            ind_Qlab = getnearest(Qz_interp,Qlab_near) # int(np.where(Qs_mnstd[:,0]==Qlab_near)[0])

            Z_LAB_Qs_km = Z_km_interp[ind_Qlab]

    #       Z_LAB_LT_km = box[i_var1,i_var2].run_info.zLAB[-1]/1e3

            # calculate the residual between the predicted and the measured:
            Res = ((Z_LAB_Qs_km - zLAB_obs_km)**2)/zLAB_obs_km
            Res_lab_Q_mat[i_var1,i_var2] = Res # np.log10(Res)
            ind_zLAB_Q_mat[i_var1,i_var2] = ind_Qlab
            Z_LAB_Q_mat[i_var1,i_var2] = Z_LAB_Qs_km

    #        Res = ((Z_LAB_LT_km- zLAB_obs_km)**2)/zLAB_obs_km
    #        Res_lab_etaLT_mat[i_var1,i_var2] = Res # np.log10(Res)

    return Res_lab_Q_mat, ind_zLAB_Q_mat, Z_LAB_Q_mat


# ==================
# function to find model that best fits AVERAGE Vs in the adiabatic part:
# Res_lab_Q_mat, ind_zLAB_Q_mat = flab.find_LAB_Res(box,Q_LAB,zLAB_obs_km,i_fmin,i_fmax)
def find_Vs_adavg_Res(box,Vs_adavg_obs,i_fmin,i_fmax):

    # get the size of the box
    dims = box.shape
    n_var1 = dims[0]
    n_var2 = dims[1]
    # MAKE ARRAYS for the residual, the depth and the index of LAB depth in each model
    Res_Vs_adavg_mat = np.zeros([n_var1,n_var2])
    Vs_adavg_mat = np.zeros([n_var1,n_var2])

    for i_var1 in range(n_var1):
        for i_var2 in range(n_var2):

            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]

            #Vs_f_mat = (frame.VBR.out.anelastic.AndradePsP.Va)/1e3
            Vs_f_mat = (frame.VBR.out.anelastic.YT_maxwell.V)/1e3
            Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
            Vs_mnstd = zip_meanstd_fband(Vs_f_band)
            # INTERPOLATE Qz and Z_km TO HIGHER RES:
            #nn_pts = 1000
            #Z_km_interp = np.linspace(Z_km[0],Z_km[-1],nn_pts)
            #Qz_interp = np.interp(Z_km_interp,Z_km,Qs_mnstd[:,0])

            # Calculate average Vs in the adiabatic region:
            z_plate = float(box[i_var1,i_var2].info.var2val)
            ind_zplate =  getnearest(Z_km,z_plate)
            Z_ad_int = 20
            ind_bot_asth = ind_zplate+Z_ad_int
            if (ind_bot_asth >= len(Z_km)):
                ind_bot_asth = -1
            Vs_adbt_vec = Vs_mnstd[ind_zplate:ind_bot_asth,0]
            Vs_adavg = np.mean(Vs_adbt_vec)

            # calculate the residual between the predicted and the measured:
            Res = ((Vs_adavg - Vs_adavg_obs)**2)/Vs_adavg_obs
            Res_Vs_adavg_mat[i_var1,i_var2] = Res # np.log10(Res)
            Vs_adavg_mat[i_var1,i_var2] = Vs_adavg

    return Res_Vs_adavg_mat, Vs_adavg_mat



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


# =========================
# function to find the plate thickness through the GIA timescale sweep:
# zplate_freq, zplate_time, i_zM, i_zLAB = flab.find_zplate_GIA_freq(Z_km,Z_find_M,zLAB_obs_km,Mod_f_mat,Mod_LAB)

def find_zplate_GIA_freq(Z_km,Z_find_M,zLAB_obs_km,Mod_f_mat,Mod_LAB,freq_vec):
    n_freqs = Mod_f_mat.shape[1]
    ilab_ara = np.zeros(len(freq_vec))
    zlab_ara = np.zeros(len(freq_vec))

    # INTERPOLATE Qz and Z_km TO HIGHER RES:
    nn_pts = 1000
    Z_km_interp = np.linspace(Z_km[0],Z_km[-1],nn_pts)

    # find the index at which Z=Z where you find the M(time) curve.
    i_zM = getnearest(Z_km_interp,Z_find_M)
    # find the index at which Z=Z_LAB.
    i_zLAB = getnearest(Z_km_interp,zLAB_obs_km)
    print(i_zLAB)
    for i_freq in range(n_freqs):
        m_f = (Mod_f_mat[:,i_freq])/1e9  # remember, these go from low freq to high freq !
        Mf_interp = np.interp(Z_km_interp,Z_km,m_f)
        del_ind = int(np.round(0.4*i_zLAB))
        i_lab = getnearest(Mf_interp[i_zLAB-del_ind:i_zLAB+del_ind],Mod_LAB)+(i_zLAB-del_ind)  # 3/26/18, extended the search below i_zLAB..
        # could be source of future bugs too!

        ilab_ara[i_freq] = i_lab
        Z_LAB_f = Z_km_interp[i_lab]
        zlab_ara[i_freq] = Z_LAB_f
        print(Z_LAB_f)

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
