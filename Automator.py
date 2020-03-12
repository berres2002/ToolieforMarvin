from unwise import unwise as uw
from astropy.table import Table
import numpy as np
import os
from toolie import MFTOOLIE as mft
import pandas as pd

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
        if temp==0:
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
                self.run(b,'full_test.csv')
            if umm=='n' or '':
                print('Exiting Program')
                return

if __name__=='__main__':
    inp=input('\n0 for Demo, anything else for full thing -> ')
    if inp =='exit':
        print('Exiting Program')
    else:
        auto(int(inp))
