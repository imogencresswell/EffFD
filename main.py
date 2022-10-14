"""
Main script for FFD generation from TESS lightcurve data


Usage:
    main.py [options]
    main.py <config> [options]

Options:
    --star_names=<star>            Names of stars  [default: [AU Mic]]

"""


import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import csv

def save_raw_lc(object):
    #function that saves lightkurve image
    plt.figure()
    pixelfile = lk.search_targetpixelfile(object)[2].download()
    lc = pixelfile.to_lightcurve(aperture_mask='all').flatten()
    lc.plot()
    plt.savefig(str(object)+'.png')

#reads in arguments from config file or command line
#also allows x = True commands
args   = docopt(__doc__)
if args['<config>'] is not None:
    config_file = Path(args['<config>'])
    config = ConfigParser()
    config.read(str(config_file))
    for n, v in config.items('parameters'):
        for k in args.keys():
            if k.split('--')[-1].lower() == n:
                if v == 'true': v = True
                args[k] = v
#PROBEMS/QUERIES:
#the search_target pixel won't find AU Mic if star name is a string... but if you put a string into the function it finds it... dont understand


#TODO:
#need to figure out what the index means in search target pixel file.. could use download all but takes forever, should be a user input but idk what it is
#need to make sure we add failsafes if there is a typo/star not found in one name i.e program doesnt error out but continues without that name


if args['--star_names'] is not None:
    star_names = str(args['--star_names'])
    star_names_list = star_names.split(',')

for star in star_names_list:
    save_raw_lc(star)


