# reading VBR BOXES .mat data structures into python:



# -----------------
# read in the Box

t0 = time()
# when you do stuff that takes time, see how long it takes !
hoozit = sio.whosmat(path+matobj) # get a little info without loading it..
print(hoozit)
# load the mat file:
b=sio.loadmat(path+matobj,squeeze_me=True,struct_as_record=False)
print(type(b))
print(b.keys())
# The box itself is a structured object (a class operating on numpy array?)
# access field names by [name]._fieldnames
box = b['Box']
print('box is a ' + str(type(box)) + ', with shape:')
print(box.shape)
n_var1 = box.shape[0]
n_var2 = box.shape[1]
print(box[0,2]._fieldnames)
print(box[0,2].info._fieldnames)
print(box[0,2].run_info._fieldnames)
# CH's has ['info', 'Movie'] -- Movie is now gone (y18-06)
# BH's has ['info', 'run_info', 'Frames']

# Test frame:
frame_test = box[0,1].Frames[-1]
print(frame_test._fieldnames)
print(frame_test.VBR._fieldnames)
print('VBR input:')
print(frame_test.VBR.input._fieldnames ) # shit !! .in does not work ! already a built in function !
print('VBR out:')
print(frame_test.VBR.out._fieldnames)



===============================================================
EXTRACT PHYSICAL PROPERTIES FROM EACH BOX
testing:

# extract the Frequency-INDEPENDENT properties we need:
eta = frame_test.VBR.out.viscous.LH2012.eta_total # composite viscosity
eta_diff = frame_test.VBR.out.viscous.LH2012.diff.eta # viscosity for diff creep
Gu = frame_test.VBR.out.elastic.anharmonic.Gu # unrelaxed shear modulus
TauMxw = eta/Gu
TauMxw_diff = eta_diff/Gu

# extract the Frequency-DEPENDENT properties we need over f band of interest:
M_f_mat = frame_test.VBR.out.anelastic.AndradePsP.Ma
Vs_f_mat = frame_test.VBR.out.anelastic.AndradePsP.Va
# dims = M_f_mat.shape  # print(dims)
M_f_band = M_f_mat[:,i_fmin:i_fmax] # slice: [start:stop:step], # dimsliced = M_f_band.shape ; print(dimsliced)
Vs_f_band = Vs_f_mat[:,i_fmin:i_fmax]
