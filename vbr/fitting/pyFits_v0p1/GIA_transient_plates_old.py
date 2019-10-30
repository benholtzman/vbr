# PLOT THE GIA predictions...
# not using notebooks !

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import sys

# ================================================
# READ THESE IN DIRECTLY FROM  THE MATLAB STRUCTURES !
#pathtodata = 'VBR_Tp1375_Z200_An_lowETA/'
#pathtodata = './SNA_fits/VBR_Tp1375_Z200_An/'
#pathtodata = './TNA_fits/VBR_Tp1450_Z80_An/'
SNAfit = True
if SNAfit==True:
    pathtodata = './SNA_fits_Qfit/VBR_Tp1325_Z200_An/'
    zLAB = 180 # plot this line as the seismically identified plate thickness..
elif SNAfit==False:
    pathtodata = './TNA_fits_Qfit/VBR_Tp1425_Z100_An/'
    zLAB = 80 # plot this line as the seismically identified plate thickness..

Z_find_M = 0.8*zLAB # find the values of M at this depth.


Mod_f_mat = np.loadtxt(pathtodata + 'M_f.txt')
dims = Mod_f_mat.shape
print(dims)
eta = np.loadtxt(pathtodata + 'eta.txt') # composite viscosity
eta_diff = np.loadtxt(pathtodata + 'eta_diff.txt') # viscosity for diff creep
G_u = np.loadtxt(pathtodata + 'G_u.txt') # unrelaxed shear modulus
TauMxw = eta/G_u
TauMxw_diff = eta_diff/G_u

Z_km = np.loadtxt(pathtodata + 'Z_km.txt') # composite viscosity
freq_vec = np.loadtxt(pathtodata + 'freqs.txt')

time_vec = 1/freq_vec[::-1] # reverse direction, short time to long time.
num_freqs = len(freq_vec)


# ================================================
#  PICK A THRESHOLD M value
# make sure it coincides at lowest freq with the eta_LAB = 10*eta_asth_avg. DEFINE
# TO DEFINE APPARENT PLATE THICKNESS !

Mod_LAB = 20.0
ilab_ara = np.zeros(len(freq_vec))
zlab_ara = np.zeros(len(freq_vec))

def getnearest(array,val):
    idx = (np.abs(array-val)).argmin()
    return idx

i_zM = getnearest(Z_km,Z_find_M)
i_zLAB = getnearest(Z_km,zLAB)

for i_freq in range(len(freq_vec)):
    m_f = (Mod_f_mat[:,i_freq])/1e9  # remember, these go from low freq to high freq !
    i_lab = getnearest(m_f[0:i_zLAB+20],Mod_LAB)  # 3/26/18, extended the search below i_zLAB..
    # could be source of future bugs too!

    ilab_ara[i_freq] = i_lab
    zlab_ara[i_freq] = Z_km[i_lab]

# flip these now so they correspond to increasing time...
ilab_ara_t = ilab_ara[::-1]
zlab_ara_t = zlab_ara[::-1]

# find the index at which Z=Z where you find the M(time) curve.
i_zM = getnearest(Z_km,Z_find_M)

# make a little array of zPlate as a function of frequency
zplate_freq = np.zeros((len(freq_vec),2))
zplate_freq[:,0] = freq_vec
zplate_freq[:,1] = zlab_ara
np.savetxt(pathtodata+'zplate_freq.txt', zplate_freq)

# ===================================================================
# ===================================================================
# PLOT and write out plate thickness vs frequency:
# for the laplace transform to time domain:
plt.figure(0)
plt.plot(freq_vec,zlab_ara)
plt.plot([freq_vec[0], freq_vec[-1]],[zLAB,zLAB], lw=10, color='gray', alpha=0.5)
plt.xlabel('Frequency (or strain rate) [Hz]')
plt.ylabel('Apparent plate thickness [km]')
plt.title('Apparent plate thickness with frequency')
plt.autoscale(enable=True, axis='y', tight=True)
axes = plt.gca()
ymin = zlab_ara[0]
ymax = zLAB+0.05*zLAB
axes.set_ylim([ymin,ymax])


# ===================================================================
# ===================================================================
# >>>>>>>>>>>>>>>>>>
# PLOT THE CURVES !
# >>>>>>>>>>>>>>>>>>
plt.figure(1)
# plasma colormap
cmap = plt.cm.get_cmap('autumn')

# =======================
#f, axs = plt.subplots(1,2,figsize=(7,7))
#plt.subplot(1, 2, 1)
L = 0.1
B = 0.2
W = 0.25
H = 0.6

ax = plt.axes([L,B,W,H])


bk = [0.0,0.0,0.0,1.0]

for i_freq,freq in enumerate(freq_vec):
    color_num = i_freq/float(num_freqs)
    #color_num = freq / max(freq_vec)
    col = cmap(color_num)
    time_line = np.log10(time_vec[-i_freq]/(np.pi*1e7)*np.ones(len(Z_km)))
    #plt.plot(time_line,Z_km,color=col,lw=2.0 )
    ax.plot(time_line,Z_km,color=col,lw=2.0 )

time_line = np.log10(time_vec[0]/(np.pi*1e7)*np.ones(len(Z_km)))
ax.plot(time_line,Z_km,color=bk,lw=1.0 )
time_line = np.log10(time_vec[-1]/(np.pi*1e7)*np.ones(len(Z_km)))
ax.plot(time_line,Z_km,color=bk,lw=1.0 )


ax.plot(np.log10(TauMxw_diff/(np.pi*1e7)),Z_km,color='green',lw=2.0, ls='--' )
ax.plot(np.log10(TauMxw/(np.pi*1e7)),Z_km,color='blue',lw=3.0 )

axes = plt.gca()
xmin = 1.0
xmax = 8.0
axes.set_xlim([xmin,xmax])

# add plate thickness line:
zLAB_line = zLAB*np.ones(len(Z_km))
ax.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line,color='black',lw=5.0, alpha=0.4)

plt.gca().invert_yaxis()
plt.title('Maxwell time')
plt.xlabel('log10 time [yrs]')
plt.ylabel('Depth [km]')
plt.autoscale(enable=True, axis='y', tight=True)

# =======================
#plt.subplot(1, 2, 2)
#plt.figure(figsize=(5,10))
L = L + 0.26
B = B
ax2 = plt.axes([L,B,W,H])
#for i_freq in range(num_freqs):
# scrolling through BACKWARDS (short to long time ! )
#for i_freq in range(num_freqs-1,0,-1):
for i_freq,freq in enumerate(freq_vec):
    #print(i_freq)
    color_num = i_freq/float(num_freqs)
    #color_num = freq / max(freq_vec)
    col = cmap(color_num)
    m_f = (Mod_f_mat[:,i_freq])/1e9
    ax2.plot(m_f,Z_km,color=col,lw=2.0 )

m_f = (Mod_f_mat[:,0])/1e9
ax2.plot(m_f,Z_km,color=bk,lw=1.0 )

m_f = (Mod_f_mat[:,num_freqs-1])/1e9
ax2.plot(m_f,Z_km,color=bk,lw=1.0 )

axes = plt.gca()
xmin = 0.0
xmax = 80.0
axes.set_xlim([xmin,xmax])

ax2.plot(np.linspace(xmin,xmax,len(Z_km)),zLAB_line,color='black',lw=5.0, alpha=0.4)

ModLAB_line = Mod_LAB*np.ones(len(Z_km[0:i_zLAB]))
ax2.plot(ModLAB_line,Z_km[0:i_zLAB],color='purple',lw=10.0, alpha=0.2)

# add M @ depth line:
z_Mt_line = Z_find_M *np.ones(len(Z_km))
ax2.plot(np.linspace(xmin,xmax,len(Z_km)),z_Mt_line,color='black',lw=1.0, alpha=0.7)


plt.gca().invert_yaxis()
plt.title('Modulus(f)')
plt.xlabel('Modulus [GPa]')
plt.ylabel('')
ax2.set_yticklabels([])
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
# ================================================
# PLOT THE EXTRACTED DATA !
# ================================================
#f0 = plt.subplots(2,1,figsize=(6,8))
up = 0.36
L = L + 0.33
B = B + up
H = 0.4*H
ax3 = plt.axes([L,B,W,H])


# >>>>>>>>>>>>>>>>>>
# PLOT THE modulus change with frequency !
## >>>>>>>>>>>>>>>>>>
#plt.subplot(2,1,1)

ax3.plot(freq_vec,(Mod_f_mat[i_zM,:])/1e9, c='black')

for i_freq,freq in enumerate(freq_vec):
    color_num = i_freq/float(num_freqs)
    #color_num = freq / max(freq_vec)
    col = cmap(color_num)
    plt.plot(freq_vec[i_freq],(Mod_f_mat[i_zM,i_freq])/1e9,c=col,
                       marker='o',
                       markersize=10, markeredgecolor = 'black' )

plt.title('M(f) at ' + str(Z_find_M) + ' km')
plt.xlabel('frequency')
plt.ylabel('M [GPa]')


# >>>>>>>>>>>>>>>>>>
# PLOT THE LAB EVOLUTION !
# >>>>>>>>>>>>>>>>>>
#plt.subplot(2,1,2)

L = L
B = B - up
H = H
ax4 = plt.axes([L,B,W,H])

sec_yr = np.pi*1e7

ax4.plot(time_vec/sec_yr/100,zlab_ara_t, c='black')
ax4.plot([time_vec[0]/sec_yr/100, time_vec[-1]/sec_yr/100],[zLAB,zLAB], color='black',lw=5.0, alpha=0.4)
axes = plt.gca()
axes.set_ylim([0.9*min(zlab_ara_t),1.1*max(zlab_ara_t)])
plt.gca().invert_yaxis()

#for i_freq in range(num_freqs-1):
for i_freq,freq in enumerate(freq_vec):
    color_num = i_freq/float(num_freqs)
    #color_num = freq / max(freq_vec)
    col = cmap(color_num)
    ax4.plot(time_vec[-i_freq]/sec_yr/100,zlab_ara_t[-i_freq],c=col,
                       marker='o',
                       markersize=10, markeredgecolor = 'black' )

plt.title('apparent plate thickness')
plt.xlabel('time [yrs x 100]')
plt.ylabel('thickness [km]')


#
plt.show()
# ============
sys.exit()
# ============
