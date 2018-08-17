import numpy as np
import matplotlib.pyplot as plt
import fit_LAB_fcns as flab

def plot_Q_profiles(layout,BoxObj,fits,obs,i_fmin,i_fmax):
    ''' plots all the Q profiles and the best fitting profiles '''
    axes_location=[layout['L'],layout['B'],layout['W'],layout['H']]
    ax = plt.axes(axes_location)
    box=BoxObj.box
    for i_var1 in range(BoxObj.n_var1):
        for i_var2 in range(BoxObj.n_var2):
            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]
            Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
            Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
            Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
            ax.plot(np.log10(Qs_mnstd[:,0]),Z_km,color='black',lw=1.0, ls='-', alpha=0.2 )

    # add standard deviations to these lines !
    temp_colors={'hot':'red','cold':'blue'}
    axis_lims={'x':[1.0,8.0],'y':[0,350]}
    axis_titles={'title':'Fit LAB depth with $Q_{LAB}$.','x':'log $Q$','y':'Depth [km]'}
    for temp in ['cold','hot']:
        # best fittin Q vs z profile
        clr=temp_colors[temp]
        ij_best=fits[temp]['ij_best']
        Z_km = box[ij_best].run_info.Z_km
        frame = box[ij_best].Frames[-1]
        Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa
        Qs_f_band = Qs_f_mat[:,i_fmin:i_fmax]
        Qs_mnstd = flab.zip_meanstd_fband(Qs_f_band)
        ax.plot(np.log10(Qs_mnstd[:,0]),Z_km,color=clr,lw=2.0, ls='-', alpha=0.9 )

        # best fititng Z LAB from best Q profile
        xmin=axis_lims['x'][0]; xmax=axis_lims['x'][1]
        Z_LAB_Q_km =  fits[temp]['LAB'][ij_best]
        zLAB_Q_line = Z_LAB_Q_km*np.ones(len(Z_km))
        print("Z_LAB_Q_km for " + temp + ' is ' + str(Z_LAB_Q_km))
        ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_Q_line,color=clr,lw=2.0, linestyle=':', alpha=0.5)

        # OBSERVED plate thickness line
        zLAB_line = obs[temp]['zLAB']*np.ones(len(Z_km))
        ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line,color=clr,lw=7.0, linestyle='-', alpha=0.2)

    ax.set_xlim(axis_lims['x'])
    ax.set_ylim(axis_lims['y'])
    plt.title(axis_titles['title'])
    plt.xlabel(axis_titles['x'])
    plt.ylabel(axis_titles['y'])
    ax.set_yticklabels([])
    ax.invert_yaxis()

    return ax


def plot_Vs_profiles(layout,BoxObj,fits,obs,i_fmin,i_fmax):

    box=BoxObj.box
    axes_location=[layout['L'],layout['B'],layout['W'],layout['H']]
    ax2 = plt.axes(axes_location)

    for i_var1 in range(BoxObj.n_var1):
        for i_var2 in range(BoxObj.n_var2):
            Z_km = box[i_var1,i_var2].run_info.Z_km
            frame = box[i_var1,i_var2].Frames[-1]
            Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va
            Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
            Vs_mnstd = flab.zip_meanstd_fband(Vs_f_band)
            ax2.plot((Vs_mnstd[:,0])/1e3,Z_km,color='black',lw=1.0, ls='-', alpha=0.2 )

    # add standard deviations to these lines !
    temp_colors={'hot':'red','cold':'blue'}
    axis_lims={'x':[4.2,5.05],'y':[0,350]}
    axis_titles={'title':'Fit $V_s$ in asth.','x':'$V_s$ [km/s]',}

    for temp in ['cold','hot']:
        # best Vs profile
        clr=temp_colors[temp]
        ymin=axis_lims['y'][0]
        ymax=axis_lims['y'][1]
        ij_best=fits[temp]['ij_best']
        Z_km = box[ij_best].run_info.Z_km
        frame = box[ij_best].Frames[-1]
        Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va
        Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
        Vs_mnstd = flab.zip_meanstd_fband(Vs_f_band)
        ax2.plot((Vs_mnstd[:,0])/1e3,Z_km,color=clr,lw=2.0, ls='-', alpha=0.9 )

        # OBSERVED Vs lines:
        Vs_MEASRD_line = obs[temp]['Vs_adavg']*np.ones(len(Z_km))
        ax2.plot(Vs_MEASRD_line,np.linspace(ymin,ymax,len(Z_km)),color=clr,lw=5.0, linestyle='-', alpha=0.2)

        # Vs lines CALCULATED
        Vs_line = fits[temp]['Vs_adavg_mat'][ij_best]*np.ones(len(Z_km))
        ax2.plot(Vs_line,np.linspace(ymin,ymax,len(Z_km)),color=clr,lw=2.0, linestyle=':', alpha=0.5)

    ax2.set_xlim(axis_lims['x'])
    ax2.set_ylim(axis_lims['y'])
    plt.title(axis_titles['title'])
    plt.xlabel(axis_titles['x'])
    ax2.set_yticklabels([])
    ax2.invert_yaxis()

    return ax2

def residPlots(layouts,layoutid,Tpot_vec,zPlate_vec,fits,temp,clr):
    ij_best=fits[temp]['ij_best']

    # =================================
    # plot the residuals for ZLAB
    layout=layouts[layoutid+'1']
    axes_locations=[layout['L'],layout['B'],layout['W'],layout['H']]
    ax3 = plt.axes(axes_locations)
    xy_rangelist = [zPlate_vec[0],zPlate_vec[-1],Tpot_vec[-1],Tpot_vec[0]] # flip vertical directions?
    ax3.imshow(np.log(fits[temp]['Res_LAB']), cmap=plt.cm.gray, extent=xy_rangelist) #, aspect="auto")
    ax3.scatter(zPlate_vec[ij_best[1]], Tpot_vec[ij_best[0]], color=clr, s=10)

    ax3.set_title('Res: $Z_{LAB}$ ('+temp.upper()+')')
    #ax2.set_xticklabels(zPlate_vec)
    #ax2.xaxis.set_major_locator(plt.LinearLocator(5))
    ax3.set_xlabel('$z_{plate}$')
    #ax2.set_yticklabels(Tpot_vec)
    #ax2.yaxis.set_major_locator(plt.LinearLocator(5))
    ax3.set_ylabel('$T_{pot}$ [C]')


    # =================================
    # plot the residuals for Vs_adavg
    layout=layouts[layoutid+'2']
    axes_locations=[layout['L'],layout['B'],layout['W'],layout['H']]
    ax4 = plt.axes(axes_locations)
    ax4.imshow(np.log(fits[temp]['Res_Vs_adavg_mat']), cmap=plt.cm.gray, extent=xy_rangelist)
    ax4.scatter(zPlate_vec[ij_best[1]], Tpot_vec[ij_best[0]], color=clr, s=10)

    ax4.set_title('Res: $V_s$')
    #ax4.set_xticklabels(zPlate_vec)
    ax4.set_xlabel('$z_{plate}$')
    #ax4.set_yticklabels(Tpot_vec)
    #ax4.set_ylabel('T_pot')
    ax4.set_yticklabels([])

    # =================================
    # plot the residuals for Vs_adavg HOT
    layout=layouts[layoutid+'3']
    axes_locations=[layout['L'],layout['B'],layout['W'],layout['H']]
    ax5 = plt.axes(axes_locations)
    #ax5.imshow(np.log10(Res_JOINT_HOT), cmap=plt.cm.gray, extent=xy_rangelist)
    ax5.imshow(fits[temp]['P_JOINT'], cmap=plt.cm.gray, extent=xy_rangelist)
    ax5.scatter(zPlate_vec[ij_best[1]], Tpot_vec[ij_best[0]], color=clr, s=10)

    ax5.set_title('Joint Prob.')
    #ax5.set_xticklabels(zPlate_vec)
    ax5.set_xlabel('$z_{plate}$')
    #ax5.set_yticklabels(Tpot_vec)
    #ax5.set_ylabel('T_pot')
    ax5.set_yticklabels([])

    return ax3,ax4,ax5

def buildLayout():

    # DEPTH PLOTS:
    L1 = 0.06
    B1 = 0.2
    W1 = 0.18
    H1 = 0.7
    Layout={'Qz':{'L':L1,'B':B1,'W':W1,'H':H1},
            'Vsz':{'L':L1+W1+0.01,'B':B1,'W':W1,'H':H1}}
            
    # RESIDUAL PLOTS:
    hdel = 0.01
    vdel = 0.02
    hdX = 14*hdel

    # top row
    Layout['rh1']={'L': Layout['Vsz']['L']+Layout['Vsz']['W']+hdel,
                   'B': B1+0.4,'W':0.3,'H':0.3}

    Layout['rh2']={'L': Layout['rh1']['L']+hdX+hdel,
                   'B': Layout['rh1']['B'],'W':Layout['rh1']['W'],
                   'H': Layout['rh1']['H']}

    Layout['rh3']={'L': Layout['rh2']['L']+hdX+hdel,
                   'B': Layout['rh1']['B'],'W':Layout['rh1']['W'],
                   'H': Layout['rh1']['H']}

    # bottom row
    Layout['rb1']={'L': Layout['Vsz']['L']+Layout['Vsz']['W']+hdel,
                   'B': B1,'W':0.3,'H':0.3}

    Layout['rb2']={'L': Layout['rb1']['L']+hdX+hdel,
                   'B': Layout['rb1']['B'],'W':Layout['rb1']['W'],
                   'H': Layout['rb1']['H']}

    Layout['rb3']={'L': Layout['rb2']['L']+hdX+hdel,
                   'B': Layout['rb1']['B'],'W':Layout['rb1']['W'],
                   'H': Layout['rb1']['H']}

    return Layout
