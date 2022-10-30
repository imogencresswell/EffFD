"""
Main script for FFD generation from TESS lightcurve data

Usage:
    main.py [options]
    main.py <config> [options]

Options:
    --root_dir=<dir>               Root directory for outputs  [default: ./]

    --star_names=<star>            Names of stars
    --spectral_type=<type>         Spectral type to search for stars
    --teff_low=<temp>              Lower limit of Teff to search [default: 2000]
    --teff_high=<temp>             Lower limit of Teff to search [default: 3500]

    --savgov_iterations=<iter>     Iterations for lc.flatten  [default: 9]
    --savgov_sigma_cutoff=<sig>    Sigma cutoff for lc.flatten  [default: 3]

"""
import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import csv
import os


def get_spectral_temp(classification):
    if classification=='M':
        return 2000, 3500  # low temp limit, high temp limit, in K
    elif classification=='K':
        return 3500, 5000
    elif classification=='G':
        return 5000, 6000
    elif classification=='F':
        return 6000, 7500
    elif classification=='A':
        return 7500, 10000
    elif classification=='B':
        return 10000, 30000
    elif classification=='O':
        return 30000, 60000


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
    #print(args)

    # PROBEMS/QUERIES:

    # TODO:
    # Need to figure out what the index means in search target pixel file..
    # could use download all but takes forever,
    # should be a user input but idk what it is.

    if args['--star_names'] is not None:
        star_names = str(args['--star_names'])
        star_names_list = list(map(str.strip, star_names.split(',')))
    else:
        # Uses temperature-range star search if no names/spectral_type given
        teff_low = int(args['--teff_low'])
        teff_high = int(args['--teff_high'])

        if args['--spectral_type'] is not None:
            spectral_type = str(args['--spectral_type'])

            # Error catching for unsupported types. Stop run if unsupported.
            if spectral_type in 'L,T,Y'.split(','):
                raise ValueError('Brown dwarfs are not yet supported.')
            elif spectral_type not in 'O,B,A,F,G,K,M'.split(','):
                raise ValueError('Please use spectral type in OBAFGKM.')

            else:
                teff_low, teff_high = get_spectral_temp(spectral_type)

        # -- Astroquery search for list of stars by temperature
        star_names_list = []


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
            os.rmdir(star_path)

        print('Operations for {} finished.\n'.format(star))


if __name__ == '__main__':
    main()
