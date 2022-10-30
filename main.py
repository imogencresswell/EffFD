"""
Main script for FFD generation from TESS lightcurve data

Usage:
    main.py [options]
    main.py <config> [options]

Options:
    --star_names=<star>            Names of stars  [default: AU Mic]
    --root_dir=<dir>               Root directory for outputs [default: ./]
    --savgov_iterations=int        Iterations for lc.flatten [default: 9]
    --savgov_sigma_cutoff=int      Sigma cutoff for lc.flatten [default: 3]

"""
import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import csv
import os


def save_raw_lc(object, save_path, filter_iter, filter_sig):
    # function that saves lightkurve image and csv file
    search_result = lk.search_targetpixelfile(object)
    if not search_result:
        raise FileNotFoundError('No results for {}.'.format(object))

    pixelfile = search_result[2].download()
    lc = pixelfile.to_lightcurve(aperture_mask='all')
    lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

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

        print('\n###############################')
        print('Starting individual star search')
        print('###############################\n')

        for star in star_names_list:
            star_path = os.path.join(str(args['--root_dir']),
                                     star.replace(' ', '_'))

            try:  # Keeps program running if folders already exist
                os.mkdir(star_path)
            except FileExistsError:
                print('A folder for {} already exists.'.format(star))

            try:  # Keeps program running if a star name is not searchable
                save_raw_lc(star,
                            star_path,
                            int(args['--savgov_iterations']),
                            int(args['--savgov_sigma_cutoff']))
            except FileNotFoundError:
                print('No search results found for {}.'.format(star))
                # If we add other functions under this and they all fail
                # we could have an os.rmdir here.

            print('Operations for {} finished.\n'.format(star))


if __name__ == '__main__':
    main()
