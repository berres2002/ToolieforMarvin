from astropy.table import Table
import pandas as pd
import numpy as np
import os
print('Reading fits file...')
f = Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
r = f[f['mngtarg1'] > 0]

print('Appending columns....')
ef=[]
en=[]
u=[]
g=[]
ar=[]
eye=[]
z=[]
pl=[]
s=r.columns['plateifu'].size
for i in range(s):
    ef.append(r.columns['nsa_petro_flux'][i][0])
    en.append(r.columns['nsa_petro_flux'][i][1])
    u.append(r.columns['nsa_petro_flux'][i][2])
    g.append(r.columns['nsa_petro_flux'][i][3])
    ar.append(r.columns['nsa_petro_flux'][i][4])
    eye.append(r.columns['nsa_petro_flux'][i][5])
    z.append(r.columns['nsa_petro_flux'][i][6])
print('Creating and Saving DataFrame...')
df=pd.DataFrame({'Plate-IFU':r.columns['plateifu'],'F':ef,'N':en,'u':u,'g':g,'r':ar,'i':eye,'z':z})
df.to_csv('/uufs/chpc.utah.edu/common/home/u6030555/data_folder/csv_data/drp_bands1.csv')
print('Done!')
