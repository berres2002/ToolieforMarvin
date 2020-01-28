from astropy.table import Table
from astropy.io import fits
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

    #e=requests.get(url1.format(ra,dec))
    #file=open('f.txt','w+')
    #file.write(e.text)
    #file.close()
    try:
        re=Table.read(url1.format(ra,dec),format='ipac')
    except:
        raise ValueError('Something went wrong searching for plate-ifu: '+plate)
    ci=re[0]['coadd_id']
    params={'coadd_id':ci,'band':4,}
    params['coaddgrp'] = params['coadd_id'][:2]
    params['coadd_ra'] = params['coadd_id'][:4]
    path = str.format(
        '{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',
        **params)
    url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + path

    #y=requests.get(url)
    #with open('7443-12701.fits', 'w+') as fd:
    #    for chunk in y.iter_content(chunk_size=128):
    #        fd.write(str(chunk))
    try:
        f2=fits.open(url)
    except:
        raise ValueError('Something went wrong calling the file for plate-ifu: '+plate)
    try:
        f2.writeto(filePath+plate+'.fits')
    except:
        raise ValueError('Something went wrong saving the file to the folder "'+filePath+'" for plate-ifu: '+plate)
    #print(f[0].header)
    #file=open('f.txt','w+')
    #file.write('')
    #file.close()
    try:
        fits.open(filePath+plate+'.fits')
    except:
        raise ValueError('Something went wrong trying to open the file for plate-ifu: '+plate)
    print('Saved successfully as: '+plate+'.fits')