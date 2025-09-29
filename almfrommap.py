import numpy as np
import healpy as hp
from pathlib import Path

filename='COM_CMB_IQU-smica_2048_R3.00_full.fits' 
file_no_ext=Path(filename).stem


polmap,header = hp.read_map(filename,field=(0,1,2),h=True)

cls,alms = hp.anafast(polmap, alm=True, use_pixel_weights=True, lmax=2048) #contains cl's, a_lm's for T,E,B

f = open("alm_E_"+file_no_ext+".dat", "w")
for l in np.arange(2049):
 for m in np.arange(-l,0):
  f.write(format(np.imag(alms[1][hp.Alm.getidx(2048,l,abs(m))]),'.19g')+'\n')
 for m  in np.arange(0,l+1):
  f.write(format(np.real(alms[1][hp.Alm.getidx(2048,l,abs(m))]),'.19g')+'\n')
f.close()

f = open("alm_B_"+file_no_ext+".dat", "w")
for l in np.arange(2049):
 for m in np.arange(-l,0):
  f.write(format(np.imag(alms[2][hp.Alm.getidx(2048,l,abs(m))]),'.19g')+'\n')
 for m  in np.arange(0,l+1):
  f.write(format(np.real(alms[2][hp.Alm.getidx(2048,l,abs(m))]),'.19g')+'\n')
f.close()

