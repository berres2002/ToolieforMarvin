#Made by Aidan Berres
import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from marvin.tools import Maps
from astropy import units as u
import astropy.constants as const
import scipy.integrate as integrate
from numpy import sinh , sqrt, pi
from astropy.table import Table
import logging as log
import os
import marvin
marvin.config.forceDbOff()
marvin.config.access = 'collab'
marvin.config.login()
marvin.config.setRelease('MPL-9')
#Change level to ERROR
log.basicConfig(filename='Errors.log',level=log.ERROR,format='%(message)s at %(asctime)s', datefmt='%m/%d/%Y %I:%M:%S %p')
class dist:
    def __init__(self,z):
        self.z=float(z)
    def findDist(self):
        dh=299792.458/cosmo.H(0)
        om=0.27
        ok=1-om
        oa=0
        #E=lambda x: sqrt(om*(1+x)**3+ok*(1+x)**2+oa)
        dc=dh*integrate.quad(lambda x: 1/(sqrt(om*(1+x)**3+ok*(1+x)**2+oa)),0,self.z)
        dm=dh*(1/sqrt(ok)*sinh(sqrt(ok)*(dc/dh)))
        dl=dm*(1+self.z)
        a=cosmo.luminosity_distance(self.z).value
        return a




        #k=(cosmo.kpc_proper_per_arcmin(self.z).value)/60
        #g=k*(1+self.z)
        #a=cosmo.luminosity_distance(self.z)

class flux:
    def __init__(self,hA,hB):
        self.hA=hA
        self.hB=hB
    def extinct(self):
        e=0.934*np.log((self.hA/self.hB)/2.86)
        return e
    def fluxFind(self):
        #$print('THIS STATEMENT IS TRUE')
        e= flux.extinct()
        a= self.hA*10**(0.4*2.468*e)
        return a
 #       else:
  #         a = self.hA * 10 ** (0.4 * 2.468 * e)
   #         return a

#Use this class

class   MFTOOLIE:

    #Initializes object like a marvin Maps file, applies the star forming mask and then averages the arrays
    def __init__(self,s,**kwargs:{'locality':'remote'}):
        try:
            self.efail=False
            self.s=s
            m=Maps(s)
            #print(m.release)
            self.m=m
            masks=m.get_bpt(return_figure=False,show_plot=False)
            hal=m.getMap('emline_gflux',channel='ha_6564')
            hbe=m.getMap('emline_gflux',channel='Hb_4862')
            man = masks['sf']['global'] #* hal.pixmask.labels_to_value('DONOTUSE')
            self.m1=hal[man]
            ham=self.m1.value *10**(-17)
            #ham = hal.mask | man
            #h1=np.array(hal)
            self.ha=np.ma.sum(ham) * u.dimensionless_unscaled
            masb = masks['sf']['global'] #* hbe.pixmask.labels_to_value('DONOTUSE')
            m2=hbe[masb]
            hbm=m2.value*10**(-17)
            #hbm = hbe.mask | masb
           # h2=np.array(hbe)
            self.hb=np.ma.sum(hbm)* u.dimensionless_unscaled
           #ratio=m.getMapRatio('emline_gflux','ha_6564','Hb_4862')
           #Mas = ~masks['sf']['global'] * ratio.pixmask.labels_to_value('DONOTUSE')
           #Rat = ratio.mask | rMas
           #hr=np.ma.array(Rat)
           #self.fRat=np.ma.sum(hr)
            self.z=float(m.dapall['nsa_zdist'])
            #self.Mpc = u.parsec *1_000_000
            f=Table.read(os.getenv('MANGA_SPECTRO_REDUX')+'/MPL-9/drpall-v2_7_1.fits',format='fits') # Change this
            r=f[f['plateifu']==m.plateifu]
            self.sm=r['nsa_elpetro_mass']
        except:
            yo='The MFTOOLIE initializer failed for Plate-Ifu '+s
            log.error(yo)
            print(yo)
            self.efail = True
            return
    #Finds Dust extinction, can use marvins ratio or calculated ratio
    def extinct(self,**kwargs:{'use_mRatio':False}):
        if self.efail:
            return 0.0
        try:
            if kwargs.get('use_mRatio') == True:
                e = 0.934 * np.log((self.fRat) / 2.86)
                return e
            else:
                e = 0.934 * np.log((float(self.ha.value) / float(self.hb.value)) / 2.86)
                return e
        except:
            s1='The extinct() function failed for plate-ifu: '+self.s
            log.error(s1)
            print(s1)
            self.efail=True
            return 0.0
    #Finds Flux using extinct, can use marvins ratio or calculated ratio
    def fluxFind(self,**kwargs:{'use_mRatio':False}):
        # $print('THIS STATEMENT IS TRUE')
        try:
            if kwargs.get('use_mRatio') == True:
                if self.efail:
                    return 0.0
                e = self.extinct(use_mRatio=True)

                if e==0.0 :
                    return 0.0
                if e=='Fail':
                    return 'Fail'
                a = self.ha * 10 ** (0.4 * 2.468 * e)
                return a
            else:
                if self.efail:
                    return 0.0
                e = self.extinct()
                if e==0.0:
                    return 0.0
                if e=='Fail':
                    return 'Fail'
                a = self.ha * 10 ** (0.4 * 2.468 * e) * u.erg / (u.cm ** 2 * u.second)
                return a
        except:
            s2='The fluxFind() function failed for plate-ifu: '+self.s
            log.error(s2)
            print(s2)
            return 0.0
    #Finds distance using Redshift, can use astropy's function or calculate with our own parmeters for
    #omega m, omega k, omega A (dark matter), and Hubble constant
    def findDist(self,**kwargs:{'calc':False}):
        try:
            if kwargs.get('calc') ==True:
                #c=299792.458* u.km / u.second
                dh=const.c/cosmo.H(0)
                om=0.27
                ok=1-om
                oa=0
                #E=lambda x: sqrt(om*(1+x)**3+ok*(1+x)**2+oa)
                dc=dh*integrate.quad(lambda x: 1/(sqrt(om*(1+x)**3+ok*(1+x)**2+oa)),0,self.z)
                dm=dh*(1/sqrt(ok))*((np.exp(sqrt(ok)*(dc/dh))-np.exp(-sqrt(ok)*(dc/dh)))/2)
                dl=dm*(1+self.z)
                return dl
            else:
                a=cosmo.luminosity_distance(self.z)#.value * u.Mpc
                return a
        except:
            s3='The findDist() function failed for plate-ifu: '+self.s
            log.error(s3)
            print(s3)
            return 0.0
    #Finds luminosity using previous methods, have not implemented 'calc' parameter for distance calculations
    def findLum(self, **kwargs:{'use_mRatio':False,'nod':False}):
        #Might have to convert distance to something else
        try:
            if kwargs.get('nod')==True:
                f = cosmo.luminosity_distance(self.z).value * u.Mpc
                x = f.to(u.cm)
                if self.efail:
                    return 0.0
                fl = self.fluxFind()
                if fl  == 0.0:
                    return 0.0
                if fl=='Fail':
                    return 'Fail'# * u.erg / (u.cm ** 2 * u.second)
                L = fl * 4 * pi * x ** 2
                return L.value
            if kwargs.get('use_mRatio')==True:
                f = cosmo.luminosity_distance(self.z)#.value * u.Mpc
                x = f.to(u.cm)
                if self.efail:
                    return 0.0
                fl = self.fluxFind(use_mRatio=True) * u.erg / (u.cm ** 2 * u.second)
                if fl== 0.0:
                    return 0.0
                if fl=='Fail':
                    return 'Fail'
                L = fl * 4 * pi * x ** 2
                return L
            else:
                f = cosmo.luminosity_distance(self.z)#.value * u.Mpc
                x = f.to(u.cm)
                if self.efail:
                    return 0.0
                fl = self.fluxFind() #* u.erg / (u.cm ** 2 * u.second)
                if fl == 0.0:
                    return 0.0
                if fl=='Fail':
                    return 'Fail'
                L = fl * 4 * pi * x ** 2
                return L
        except:
            es='The findLum() function failed for plate-ifu: '+self.s
            print(es)
            log.error(es)
            return 0.0
    def findSFR(self):
        try:
            if self.efail:
                return 0.0
            L=self.findLum(nod=True)
            if L is None:
                return 0.0
            sf=np.log10(L)-41.27
            sfr=10**sf
            #M_sun/year
            return sfr
        except:
            es='The findSFR() function failed for plate-ifu: '+self.s
            print(es)
            log.error(es)
            return 0.0
    def findSSFR(self):
        try:
            if self.efail:
                return 0.0
            s=self.findSFR()/float(self.sm)
            #1/year
            if s==0.0:
                return 0.0
            return s
        except:
            yuh='The findSSFR() function failed for plate-ifu: '+self.s
            log.error(yuh)
            print(yuh)
            return 0.0
#MFTOOLIE PARTE DEUX using 'masks' function instead of masking non SF spaxels
class MFT2:
    def __init__(self,s):
        try:
            self.efail=False
            self.s=s
            m=Maps(s)
            #print(m.release)
            self.m=m
            hal=m.getMap('emline_gflux',channel='ha_6564')
            hbe=m.getMap('emline_gflux',channel='Hb_4862')
            self.ha=hal.masked.sum()
            self.hb=hbe.masked.sum()
            # add more if needed
        except:
            yo='The MFTOOLIE2 initializer failed for Plate-Ifu '+s
            log.error(yo)
            print(yo)
            self.efail = True

    def extinct(self,**kwargs:{'use_mRatio':False}):
        if self.efail:
            return 999.0
        try:
            if kwargs.get('use_mRatio') == True:
                e = 0.934 * np.log((self.fRat) / 2.86)
                return e
            else:
                e = 0.934 * np.log((float(self.ha) / float(self.hb)) / 2.86)
                return e
        except:
            s1='The extinct() function failed for plate-ifu: '+self.s
            log.error(s1)
            print(s1)
            self.efail=True
            return 999.0
