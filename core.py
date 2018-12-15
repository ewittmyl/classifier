import glob
import os
from functions import *
from pathlib import Path

def main(path):
    """
    This is the main script for operating on the 
    pipeline data and performing machine learning to 
    classify between stars and galaxies from the detections.
    ------------------------------------------------------------
    path: string
        path of the pipelines
    """

    # define variable for all the ovbservatng dates
    all_dates = glob.glob(path+"20*")
    all_dates = [date[-10:] for date in all_dates]

    # only run if the dates have not been processed
    try:
        processed_dates = os.listdir("./results")
    except:
        processed_dates = []
    processing = list(set(all_dates)-set(processed_dates))
    for date in processing:
        date_dir = Path("./"+date)
        if date_dir.exists() != 1:
            Popen('mkdir '+date, shell=True)
            for i in ['1','2','3','4']:
                Popen('mkdir '+date+'/UT'+i, shell=True)
            images = glob.glob(path+date+"/final/*-median.fits")
            images = [i.split("/")[-1] for i in images]
            for i, image in enumerate(images):
                print("processing: {}/{}".format(i+1,len(images)))
                prep(path,date,image)
                # only process the image if check_image returns 1 as a good image
                try:
                    if check_image(image) == 1:
                        sex(image)
                        gen_tab(image[:-5]+'.cat')
                        P = Popen("rm -rf "+image, shell=True)
                        P = Popen("rm -rf "+image[:-5]+".cat", shell=True)
                        P.wait()
                        for i in ['1','2','3','4']:
                            try:
                                if image.find("UT"+i) != -1:
                                    P = Popen("mv "+image[:-5]+".csv "+date+"/UT"+i, shell=True)
                                    P.wait()
                            except:
                                pass
                    else:
                        P = Popen('rm -rf '+image, shell=True)
                        P.wait()
                except:
                    print("Error occurs when processing: {}".format(image))
                    print("Skip processing: {}".format(image))
                    P = Popen('rm -rf '+image, shell=True)
                    P.wait()
                    P = Popen('rm -rf '+image[:-5]+'.cat', shell=True)
                    P.wait()
        elif date_dir.exists() == 1:
            all_date_files = glob.glob(path+date+"/final/*-median.fits")
            all_date_files = [i.split("/")[-1].split(".")[0] for i in all_date_files]
            processed_date_files = glob.glob(date+"/*/*")
            processed_date_files = [i.split("/")[-1].split(".")[0] for i in processed_date_files]
            to_be_processed = list(set(all_date_files)-set(processed_date_files))
            for i, not_yet_processed in enumerate(to_be_processed):
                print("processing: {}/{}".format(i+1,len(to_be_processed)))
                image = not_yet_processed + ".fits"
                prep(path,date,image)
                # only process the image if check_image returns 1 as a good image
                try:
                    if check_image(image) == 1:
                        sex(image)
                        gen_tab(image[:-5]+".cat")
                        P = Popen("rm -rf "+image, shell=True)
                        P = Popen("rm -rf "+image[:-5]+".cat", shell=True)
                        P.wait()
                        for i in ['1','2','3','4']:
                            try:
                                if image.find("UT"+i) != -1:
                                    P = Popen("mv "+image[:-5]+".csv "+date+"/UT"+i, shell=True)
                                    P.wait()
                            except:
                                pass
                    else:
                        P = Popen('rm -rf '+image, shell=True)
                        P.wait()
                except:
                    print("Error occurs when processing: {}".format(image))
                    print("Skip processing: {}".format(image))
                    P = Popen('rm -rf '+image, shell=True)
                    P.wait()
                    P = Popen('rm -rf '+not_yet_processed+'.cat', shell=True)
                    P.wait()
        P = Popen("mv "+date+" ./results", shell=True)
        P.wait()
    merge_cat("UT4")



