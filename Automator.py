from unwise import unwise as uw
from astropy.table import Table
import numpy as np
import os
from toolie import MFTOOLIE as mft
import pandas as pd
from datetime import datetime as dt

class auto:

    def run(self,plates,file):
        dist = []
        dist = np.array(dist)
        lum = []
        lum = np.array(lum)
        sfr = []
        sfr = np.array(sfr)
        imfl = []
        imfl = np.array(imfl)
        scales = [1.5, 2, 2.5, 3]  # Radius, 1.5,2,2.5,3
        s1 = []
        s1 = np.array(s1)
        b = []
        b = np.array(b)
        pl = []
        pl = np.array(pl)
        print('Starting For loop for an array of '+str(len(plates))+' elements')
        for i in range(len(plates)):
            gal = mft(plates[i],locality='local') #Make sure to change locality to local when run in Utah
            d = gal.findDist()
            dist = np.append(dist, np.full(8, d))
            pl = np.append(pl, np.full(8, plates[i]))
            l = gal.findLum()
            lum = np.append(lum, np.full(8, l))
            strf = gal.findSFR()
            sfr = np.append(sfr, np.full(8, strf))
            axy = uw(plates[i], 3)
            imfl = np.append(imfl, axy.getApFlux(scales))
            b = np.append(b, [axy.band, axy.band, axy.band, axy.band])
            s1 = np.append(s1, scales)
            axy = uw(plates[i], 4)
            imfl = np.append(imfl, axy.getApFlux(scales))
            b = np.append(b, [axy.band, axy.band, axy.band, axy.band])
            s1 = np.append(s1, scales)
        data = pd.DataFrame(
            {'Plate-IFU': pl, 'Distance': dist, 'Luminosity': lum, 'SFR': sfr, 'WISE Band': b, 'Scale of radius': s1,
             'Flux': imfl},columns=['Plate-IFU','Distance','Luminosity','SFR','WISE Band','Scale of radius','Flux'])
        data.to_csv(file,index=False)
        print('Data table saved as = '+file)
    def __init__(self,temp):
        if temp=='0':
             f = Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
             r = f[f['mngtarg1'] > 0]
             # ra=np.random.randint(0,len(r))
             k = 0
             a = []
             while k < 5:
                 a.append(r[np.random.randint(0, len(r))]['plateifu'])
                 k = k + 1
             a = np.array(a)
             b = a.astype(str)
             self.run(b,'/uufs/chpc.utah.edu/common/home/u6030555/data_folder/csv_data/test1.csv')
             print('DONE!')
        else:
            umm=input('****WARNING: You are about to run the automator for the WHOLE drpall file****\n'
                      'Do you want to continue (y/[n])?')
            if umm=='y':
                f = Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
                r = f[f['mngtarg1'] > 0]
                k = 0
                a = []
                for i in range(len(r)):
                    a.append(r[i]['plateifu'])
                a = np.array(a)
                b = a.astype(str)
                self.run(b,'/uufs/chpc.utah.edu/common/home/u6030555/data_folder/csv_data/full_MPL9.csv')
            if umm=='n' or '':
                print('Exiting Program')
                return
class auto2:
    def run(self,plates,file):
        t1=dt.now()
        imfl = []
        imfl = np.array(imfl)
        scales = [1.5, 2, 2.5, 3]
        pl = []
        pl = np.array(pl)
        s1 = np.array([])
        backg=np.array([])
        efl=np.array([])
        print('Starting For loop for an array of '+str(len(plates))+' elements')
        for i in range(len(plates)):
            wis=uw(plates[i],4)
            fl,bg = wis.getApFlux(scales,bg=True)
            imfl=np.append(imfl,fl)
            backg=np.append(backg,bg)
            efl=np.append(efl,fl-bg)
            s1=np.append(s1,scales)
            pl = np.append(pl, np.full(4, plates[i]))
        data = pd.DataFrame(
            {'Plate-IFU': pl, 'Scale_of_radius': s1,
             'flux_dens': imfl,'bg_flux_dens':backg,'fd_diff':efl})
        data.to_csv(file,index=False)
        t2=dt.now()
        print('Data table saved as = '+file)
        print('This took',t2-t1)
    def __init__(self):
        f = Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
        r = f[f['mngtarg1'] > 0]
        # ra=np.random.randint(0,len(r))
        #k = 0
        a = []
        f = Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
        r = f[f['mngtarg1'] > 0]
        #k = 0
        a = []
        for i in range(len(r)):
            a.append(r[i]['plateifu'])
        a = np.array(a)
        b = a.astype(str)
        run(b,'/uufs/chpc.utah.edu/common/home/u6030555/data_folder/csv_data/wise_test1.csv')
if __name__=='__main__':
    inp=input('\nType 0 for Demo, "exit" to exit, and 1 for the full thing and 2 for something special -> ')
    if inp =='exit':
        print('Exiting Program')
    if inp=='1':
        auto(inp)
    if inp=='2':
        auto2()
    else:
        print('Exiting Program')
