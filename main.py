"""
Main script for FFD generation from TESS lightcurve data

Usage:
    main.py [options]
    main.py <config> [options]

Options:
    --star_names=<star>            Names of stars  [default: AU Mic]
    --root_dir=<dir>               Root directory for outputs [default: ./]

"""
import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import csv
import os


def save_raw_lc(object, save_path):
    # function that saves lightkurve image and csv file
    search_result = lk.search_targetpixelfile(object)
    if not search_result:
        raise FileNotFoundError('No results for {}.'.format(object))

    pixelfile = search_result[2].download()
    lc = pixelfile.to_lightcurve(aperture_mask='all').flatten(niters=9,
                                                              sigma=2)
    # Added more strenuous filter for flatten. It is a 9-iteration savgov
    # which throws away data points above 2sigma, rather than the default
    # 3-iter, 3-sigma filter.
    # We might want to add these as options in the config.

    plt.figure()
    lc.plot()

    save_string = save_path + '/' + object.replace(' ', '_')
    lc.to_csv(save_string+'.csv', overwrite=True)
    plt.savefig(save_string+'.png')


def main():
    # Reads in arguments from config file or command line.
    # Also allows x = True commands.
    args = docopt(__doc__)
    if args['<config>'] is not None:
        config_file = Path(args['<config>'])
        config = ConfigParser()
        config.read(str(config_file))
        for n, v in config.items('parameters'):
            for k in args.keys():
                if k.split('--')[-1].lower() == n:
                    if v == 'true':
                        v = True
                    args[k] = v

    # PROBEMS/QUERIES:

    # TODO:
    # Need to figure out what the index means in search target pixel file..
    # could use download all but takes forever,
    # should be a user input but idk what it is.

    if args['--star_names'] is not None:
        star_names = str(args['--star_names'])
        star_names_list = list(map(str.strip, star_names.split(',')))

    print('Starting individual star search.\n')
    for star in star_names_list:
        star_path = os.path.join(str(args['--root_dir']),
                                 star.replace(' ', '_'))

        try:  # Keeps program running if folders already exist
            os.mkdir(star_path)
        except FileExistsError:
            print('A folder for {} already exists.'.format(star))

        try:  # Keeps program running if a star name is not searchable
            save_raw_lc(star, star_path)
        except FileNotFoundError:
            print('No search results found for {}.'.format(star))

        print('Operations for {} finished.\n'.format(star))


if __name__ == '__main__':
    main()
