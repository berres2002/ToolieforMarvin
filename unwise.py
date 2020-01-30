from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
#from astropy.coordinates import SkyCoord
#import requests
#from PIL import Image
#from io import BytesIO

filePath='./fits_files/'

def getFits(plate):
    f=Table.read('drpall-v2_4_3.fits')
    try:
        r=f[f['plateifu']==plate]
        ra, dec = float(r['ifura']), float(r['ifudec'])
    except:
        raise ValueError('Something went wrong accessing the drpall file for plate-ifu: '+plate)


    url1='https://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd?POS={},{}'

    try:
        re=Table.read(url1.format(ra,dec),format='ipac')
    except:
        raise ValueError('Something went wrong searching for plate-ifu: '+plate)
    ci=re[0]['coadd_id']
    params={'coadd_id':ci,'band':4,}
    params['coaddgrp'] = params['coadd_id'][:2]
    params['coadd_ra'] = params['coadd_id'][:4]
    path = str.format(
        '{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits?',
        **params)
    path=path+str.format('center={},{}&size=200pix',ra,dec)

    url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + path

    try:
        f2=fits.open(url)
    except:
        raise ValueError('Something went wrong calling the file for plate-ifu: '+plate)
    try:
        f2.writeto(filePath+plate+'.fits')
    except:
        raise ValueError('Something went wrong saving the file to the folder "'+filePath+'" for plate-ifu: '+plate)

    try:
        fits.open(filePath+plate+'.fits')
    except:
        raise ValueError('Something went wrong trying to open the file for plate-ifu: '+plate)
    print('Saved successfully as: '+plate+'.fits')

def viewImage(fN):
    try:
        image_file = get_pkg_data_filename(filePath+fN)
    except:
        raise ValueError('Something went wrong accessing the file,\n '
                         'make sure the .fits file is downloaded, you can use findFits() to get the image')
    image_data = fits.getdata(image_file, ext=0)
    plt.imshow(image_data, cmap='gray')
    plt.title(fN)
    plt.colorbar()