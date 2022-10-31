"""
Main script for FFD generation from TESS lightcurve data

Usage:
    main.py [options]
    main.py <config> [options]

Options:
    --root_dir=<dir>               Root directory for outputs  [default: ./]

    --star_names=<star>            Names of stars
    --spectral_type=<type>         Spectral type to search for stars
    --teff_low=<temp>              Low Teff limit to search [default: 2000]
    --teff_high=<temp>             High Teff limit to search [default: 3500]

    --savgov_iterations=<iter>     Iterations for lc.flatten  [default: 9]
    --savgov_sigma_cutoff=<sig>    Sigma cutoff for lc.flatten  [default: 3]

"""
import os
import glob
import utils as ut
from pathlib import Path
from docopt import docopt
from configparser import ConfigParser


def main():
    # Reads in arguments from config file or command line.
    # Also allows x = True commands.
    print('Loading parameters...')
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

    ###
    ### Create list of stars to search for
    ###
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
                teff_low, teff_high = ut.get_spectral_temp(spectral_type)

        # ASTROQUERY HERE
        star_names_list = []

    print('\n###############################')
    print('Starting individual star search')
    print('###############################\n')
    for placement, star in enumerate(star_names_list):
        print('Starting on {} ({}/{}).'.format(star,
                                               placement+1,
                                               len(star_names_list)))

        star_path = os.path.join(str(args['--root_dir']),
                                 star.replace(' ', '_'))

        try:  # Keeps program running if folders already exist
            os.mkdir(star_path)
        except FileExistsError:
            print('A folder for {} already exists.'.format(star))

        try:  # Keeps program running if a star name is not searchable
            ut.save_raw_lc(star,
                           star_path,
                           int(args['--savgov_iterations']),
                           int(args['--savgov_sigma_cutoff']))
        except FileNotFoundError:
            print('No search results found for {}.'.format(star))
            print('Continuing on...\n')
            os.rmdir(star_path)
            continue

        lc_path_list = glob.glob(os.path.join(star_path, '*.csv'))
        for lc_path in lc_path_list:
            ut.analyze_lc(star, lc_path)

        print('Operations for {} finished.\n'.format(star))


if __name__ == '__main__':
    main()
