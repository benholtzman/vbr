
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from time import time
from imp import reload  # this lets us reload a module without restarting the kernel !
import pandas as pd
import sys

fi ='./VBR_GIA_LUT.mat'

lutX=sio.loadmat(fi,squeeze_me=True,struct_as_record=False)
vbr = lutX['VBR']
print(vbr._fieldnames)
print(vbr.out.elastic._fieldnames)
print(vbr.out.anelastic._fieldnames)
print(vbr.out.viscous._fieldnames)
print(vbr.out.viscous.LH2012._fieldnames)
eta_total = vbr.out.viscous.LH2012.eta_total
print(eta_total)
