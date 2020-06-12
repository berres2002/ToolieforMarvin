from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from photutils import aperture as ap
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
#from datetime import datetime as dt
import logging as log
import numpy as np
import os
#filePath='./fits_files/'
#errtxt='ERR_MSG = {}, OCCURRED AT TIME = {}'
log.basicConfig(filename='Errors.log',level=log.ERROR,format='%(message)s at %(asctime)s', datefmt='%m/%d/%Y %I:%M:%S %p')

class unwise:
    def __init__(self,plate,band):
        try:
            warnings.simplefilter('ignore', AstropyWarning)
            warnings.simplefilter('ignore', ResourceWarning)
            self.plate=plate
            self.band=band
            self.path=os.getenv('WISE_DATA')+'/' #Make sure this is the correct folder
            self.fN = self.plate +'_'+str(self.band)+ '.fits'
            self.fP = self.path + self.fN
            fT= Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits')
            try:
                self.r = fT[fT['plateifu'] == self.plate]
                self.ra, self.dec = float(self.r['ifura']), float(self.r['ifudec'])
            except:
                e2='Something went wrong accessing the drpall file for plate-ifu: ' + self.plate
                log.error(e2)
                print(e2)
                return
        except:
            e1='The unwise initializer failed for plate-ifu: '+plate+' band: '+str(band)
            log.error(e1)
            print(e1)
            return
    def getFits(self):
        try:
            f=fits.open(self.fP)
            bruh=self.fN+' already exists at '+self.path
            print(bruh)
            log.error(bruh)
            f.close()
            return
        except:
            url1='https://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd?POS={},{}'

            try:
                re=Table.read(url1.format(self.ra,self.dec),format='ipac')
            except:
                e3='Something went wrong searching for plate-ifu: '+self.plate
                log.error(e3)
                print(e3)
                return
            ci=re[0]['coadd_id']
            params={'coadd_id':ci,'band':self.band,}
            params['coaddgrp'] = params['coadd_id'][:2]
            params['coadd_ra'] = params['coadd_id'][:4]
            path = str.format(
                '{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits?',
                **params)
            path=path+str.format('center={},{}&size=200pix',self.ra,self.dec)

            url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + path

            try:
                f2=fits.open(url)
            except:
                e4='Something went wrong calling the file for plate-ifu: '+self.plate+' in band '+str(self.band)
                log.error(e4)
                print(e4)
                return
            try:

                f2.writeto(self.fP)
            except:
                e5='Something went wrong saving the file to the folder "'+self.path+'" for plate-ifu: '+self.plate+' in band '+str(self.band)
                log.error(e5)
                print(e5)
                return
            try:
               g= fits.open(self.fP);
               g.close()
            except:
                e6=('Something went wrong trying to open the file for plate-ifu: '+self.plate+' in band '+str(self.band))
                log.error(e6)
                print(e6)
                return
            #print('Saved successfully as: '+self.fN)

    def viewImage(self):
        bd={1:r'$3.4 \mu m$',2:r'$4.6 \mu m$',3:r'$12 \mu m$',4:r'$22 \mu m$'}
        try:
            image_file = get_pkg_data_filename(self.fP)
        except:
            raise ValueError('Something went wrong accessing the file,\n '
                             'make sure the .fits file is downloaded and is in the correct band, you can use .getFits() to get the image')
        image_data = fits.getdata(image_file, ext=0)
        plt.imshow(image_data, cmap='gray')
        plt.grid(False)
        plt.title(self.plate+' Band-'+str(self.band)+" ("+bd[self.band]+')')
        plt.colorbar()
    def getApFlux(self,n,**kwargs:{'plot':False,'bg':False}):
         if kwargs.get('bg'):
            posbg=SkyCoord(ra=self.ra,dec=self.dec-0.01,unit='deg')
            abg=float(self.r['nsa_elpetro_th50_r'])*2.0
            apbg=ap.SkyEllipticalAperture(posbg,abg*u.arcsec,abg*u.arcsec)
            f=fits.open(self.fP)
            self.viewImage()
            fw=WCS(f[0].header)
            apbg.to_pixel(fw).plot(color='red')
            flux=ap.aperture_photometry(f[0],apbg)
            f.close()
            return flux['aperture_sum'][0]
        aper=[]
        self.getFits()
        #try:
        for i in range(len(n)):
            a=float(self.r['nsa_elpetro_th50_r'])*n[i]
            b=a*float(self.r['nsa_elpetro_ba'])
            phi=float(self.r['nsa_elpetro_phi'])*u.deg
                #print(a,b,phi,self.ra,self.dec)
            pos=SkyCoord(ra=self.ra,dec=self.dec,unit='deg')
            aper.append(ap.SkyEllipticalAperture(pos,a*u.arcsec,b*u.arcsec,theta=phi))
            #return aper
        try:
            f=fits.open(self.fP)
            flux=ap.aperture_photometry(f[0],aper)
        except:
            er='Something went wrong in getAPFlux() for the file "'+self.fN+'" for scale radius='+str(n)
            log.error(er)
            print(er)
            bad=[999.0,999.0,999.0,999.0]
            return bad
            #return flux[0]['aperture_sum_0']
        if kwargs.get('plot'):
            self.viewImage()
            fw=WCS(f[0].header)
            aper[0].to_pixel(fw).plot(color='blue')
            aper[1].to_pixel(fw).plot(color='red')
            aper[2].to_pixel(fw).plot(color='green')
            aper[3].to_pixel(fw).plot(color='yellow')
            plt.savefig('/uufs/chpc.utah.edu/common/home/sdss09/mangawork/users/u6030555/'+'wise_pdf/'+self.plate +'_'+str(self.band)+'.pdf')
        if self.band==3:
            arr=np.array([flux[0]['aperture_sum_0'],flux[0]['aperture_sum_1'],flux[0]['aperture_sum_2'],flux[0]['aperture_sum_3']])
            b3=arr*1.8326e-06
            f.close()
            return b3
        if self.band==4:
            arr=np.array([flux[0]['aperture_sum_0'],flux[0]['aperture_sum_1'],flux[0]['aperture_sum_2'],flux[0]['aperture_sum_3']])
            b4=arr*5.2269E-05
            f.close()
            return b4
     #   except:
     #       er='Something went wrong in getAPFlux() for the file "'+self.fN+'" for scale radius='+str(n)
     #       log.error(er)
     #       print(er)
     #       u=[0.0,0.0,0.0,0.0]
     #       return u

#/uufs/chpc.utah.edu/common/home/sdss09/mangawork/users/u6030555/
