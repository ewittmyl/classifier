from subprocess import Popen
import pandas as pd
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from pathlib import Path
import glob
from astropy.io import fits


def prep(path,date,image):
    """
    This is used to copy the image to the present directory for 
    processing. Including funpack the .fz file.
    ------------------------------------------------------------
    path: string
        input path for the pipeline
    date: string
        input date for processing
    image: string
        input image in fits format
    """
    
    # run bash code with 'Popen'
    P = Popen('cp '+path+date+'/final/'+image+' ./', shell=True)
    P.wait()
    P = Popen('mv '+image+' '+image+'.fz', shell=True)
    P.wait()
    P = Popen('funpack *.fz', shell=True)
    P.wait()
    P = Popen('rm -rf *.fz', shell=True)
    P.wait()

def check_image(image, psf_criteria=6):
    """
    Check the FWHM of all images at a given date to see which images are good images. 
    (Criteria of being a good image: median of FWHM < psf_median)
    -----------------------------------------------------------------------------------
    image: string
        input image
    psf_criteria: float
        PSF criteria for classifying a good image (in pixel)
    """
    
    print('checking image: {}'.format(image))
    tab = pd.DataFrame(np.array(Table.read(image, format='fits')))
    # check the mean of FWHM < psf_criteria for good image
    if tab['FWHM_IMAGE'].mean() < psf_criteria:
        return 1
    else:
        return 0

def crop(image, x_low=0.3, x_up=0.7, y_low=0.3, y_up=0.7):
    """
    Cropping image with boundaries arguments in terms of the 
    percentage of the image coordinates.
    -------------------------------------------------------------
    table: variable
        input pandas table 
    x_low: float
        percentage of lower boundary of x-coordinate (image)
    x_up: float
        percentage of upper boundary of x-coordinate (image) 
    y_low: float
        percentage of lower boundary of y-coordinate (image)
    y_up: float
        percentage of upper boundary of y-coordinate (image) 
    """

    x_l, x_h = image['x'].max() * x_low, image['x'].max() * x_up
    y_l, y_h = image['y'].max() * y_low, image['y'].max() * y_up
    image = image[(image.x > x_l) & (image.x < x_h)]
    image = image[(image.y > y_l) & (image.y < y_h)]
    return image

def sex(image):
    """
    Re-run sextractor to generate necessary parameters for creating feature table.
    --------------------------------------------------------------------------------
    image: string
      input image
    """
   
    # run sextractor with different default parameters
    print('running SExtractor to {}...'.format(image))
    P = Popen('sextractor -c goto.sex '+image+' -CATALOG_NAME '+image[:-5]+'.cat', shell=True)
    P.wait()

def gen_tab(cat):
    """
    Generate the feature tables by combining all .cat made by SExtractor.
    --------------------------------------------------------------------------
    file: string
      input file
    """

    col = ['FLUX_APER2','FLUX_APER4','FLUX_APER5','FLUX_APER8','FLUX_APER10','FLUX_APER14',
            'MAG_APER2','MAG_APER4','MAG_APER5','MAG_APER8','MAG_APER10','MAG_APER14',
            'MAG_AUTO','MAG_PETRO','KRON_RADIUS',
            'PETRO_RADIUS','FLUX_MAX','ISOAREAF_IMAGE','x',
            'y','ra','dec','X2_IMAGE','Y2_IMAGE','XY_IMAGE',
            'THETA_IMAGE','X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE','AWIN_IMAGE','BWIN_IMAGE',
            'THETAWIN_IMAGE','AWIN_WORLD','BWIN_WORLD','THETAWIN_WORLD',
            'MU_MAX','FLAGS','FWHM_IMAGE','ELONGATION','SEX_CLASS','FLUX_RADIUS25',
            'FLUX_RADIUS50','FLUX_RADIUS85','FLUX_RADIUS95','FLUX_RADIUS99']
    print('generating features table: {}'.format(cat))
    tab = pd.read_table(cat,skiprows=41,sep=r'\s+',header=None, names=col)

    # crop the image for just using the central part of the image
    tab = crop(tab)

    # add concentration column by subtracting mag10 by mag5, rejecting the detections with negative concentration
    tab['CONCENT'] = tab.MAG_APER5 - tab.MAG_APER10
    tab = tab[tab.CONCENT > 0]

    # normalizing the columns
    print('normalizing features...')
    seesq_norm = ['X2_IMAGE','Y2_IMAGE','X2WIN_IMAGE',
                  'Y2WIN_IMAGE','XY_IMAGE','XYWIN_IMAGE',
                  'ISOAREAF_IMAGE']
    see_norm = ['AWIN_WORLD','AWIN_WORLD','FWHM_IMAGE',
                'KRON_RADIUS','PETRO_RADIUS','FLUX_RADIUS25',
                'FLUX_RADIUS50','FLUX_RADIUS85',
                'FLUX_RADIUS95','FLUX_RADIUS99']
    mag_norm = ['MAG_APER4','MAG_APER5','MAG_APER8',
                'MAG_APER10','MAG_APER14','MAG_AUTO',
                'MAG_PETRO','MU_MAX','CONCENT']
    flux_norm = ['FLUX_APER2','FLUX_APER4','FLUX_APER5',
                 'FLUX_APER8','FLUX_APER10','FLUX_APER14']
    fwhm_mean = tab.FWHM_IMAGE.mean()
    for seesq_col in seesq_norm:
        tab[seesq_col] = tab[seesq_col] / (fwhm_mean**2)
    for see_col in see_norm:
        tab[see_col] = tab[see_col] / fwhm_mean
    for mag_col in mag_norm:
        tab[mag_col] = tab[mag_col] * tab['MAG_APER2']
    for flux_col in flux_norm:
        tab[flux_col] = tab[flux_col] * tab['FLUX_MAX']
    tab['CONCENT'] = -1 * tab['CONCENT']

    # add column for galactic latitude
    print('calculating galactic latitude...')
    ra = np.array(tab['ra'].values)
    dec = np.array(tab['dec'].values)
    pos = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    tab['b'] = list(pos.galactic.b.deg)

    tab.drop(['MAG_APER2','FLUX_MAX','x','y'], axis=1, inplace=True)
    tab.to_csv(cat[:-4]+'.csv', index=False, header=True)


def merge_cat(UT="UT4"):
    """
    This is the code for merging all individual catalog of each image together.
    """
    csv_path = Path("./catalog"+UT+".csv")
    if csv_path.exists() != 1:
        Popen('rm -rf merged'+UT+'.log', shell=True)
        Popen('touch merged'+UT+'.log', shell=True)
        all_files = glob.glob("./results/20*/"+UT+"/*")
        print('merging table: {} (1/{})'.format(all_files[0],len(all_files)))
        tab = pd.read_csv(all_files[0])
        cat = tab.copy()
        merged = open('merged'+UT+'.log','a+')
        merged.write(all_files[0]+'\n')
        try:
            for i, file in enumerate(all_files[1:]):
                print('merging table: {} ({}/{})'.format(file,i+2,len(all_files)))
                tab = pd.read_csv(file)
                cat = pd.merge(cat, tab, how='outer')
                merged.write(file+'\n')
            cat.to_csv('catalog'+UT+'.csv', index=False, header=True)
            merged.close()
        except:
            cat.to_csv('catalog'+UT+'.csv', index=False, header=True)
            merged.close()
    else:
        cat = pd.read_csv('catalog'+UT+'.csv')
        all_files = glob.glob("./results/20*/"+UT+"/*")
        merged = list(pd.read_table('merged'+UT+'.log', header=None).values)
        merged = [i[0] for i in merged]
        if set(all_files) == set(merged):
            print('GOOD NEWS: No new table is needed to be merged.')
        else:
            non_processed = list(set(all_files) - set(merged))
            merged = open('merged'+UT+'.log','a+')
            try:
                for i, new_img in enumerate(non_processed):
                    print('merging table: {} ({}/{})'.format(new_img,i+1,len(non_processed)))
                    tab = pd.read_csv(new_img)
                    cat = pd.merge(cat, tab, how='outer')
                    merged.write(new_img+'\n')
                cat.to_csv('catalog'+UT+'.csv', index=False, header=True)
                merged.close()
            except:
                cat.to_csv('catalog'+UT+'.csv', index=False, header=True)
                merged.close()
    cat = pd.read_csv('catalog'+UT+'.csv')
    m = Table(cat.values, names=cat.columns)
    hdu = fits.table_to_hdu(m)
    hdulist = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdulist.writeto('catalog'+UT+'.fits', overwrite=True)



        





