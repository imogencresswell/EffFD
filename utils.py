import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table


def get_spectral_temp(classification):
    if classification == 'M':
        return 2000, 3500  # low temp limit, high temp limit, in K
    elif classification == 'K':
        return 3500, 5000
    elif classification == 'G':
        return 5000, 6000
    elif classification == 'F':
        return 6000, 7500
    elif classification == 'A':
        return 7500, 10000
    elif classification == 'B':
        return 10000, 30000
    elif classification == 'O':
        return 30000, 60000
    else:
        raise ValueError('Improper spectral type given.')


def save_raw_lc(object, save_path, filter_iter, filter_sig):
    # function that saves lightkurve image and csv file

    # SPOC == TESS data pipeline
    # Getting only the 120 second exposure light curves for consistency
    search_result = lk.search_lightcurve(object, author='SPOC', exptime=120)
    if not search_result:
        raise FileNotFoundError('No results for {}.'.format(object))

    for result in search_result:
        lc = result.download()
        lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

        ### Leaving the older way of getting the LC from the pixel file
        ### for a bit.
        # pixelfile = search_result[window].download()
        # lc = pixelfile.to_lightcurve(aperture_mask='all')
        # lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

        plt.figure()
        lc.plot()

        save_string = '{}/{}_{}'.format(save_path,
                                        object.replace(' ', '_'),
                                        result[0].mission[0][-2:])  # sector
        lc.to_csv(save_string+'.csv', overwrite=True)
        plt.savefig(save_string+'.png')
        plt.close()


def analyze_lc(csv_path):
    """Note: currently this creates toy data for each light curve
    inputted. The main flare-finding routine is based on a not-yet-published
    paper, and we are not able to share it. We have plans to implement other
    methods soon, but the toy data allows us to move forward with code
    development. Thank you for your understanding.

    """
    lc = ascii.read(csv_path, guess=False, format='csv')

    #
    # FLARE FINDING METHOD GOES HERE
    #

    # Toy data input
    flare_tbl = Table()
    flare_tbl['energy'] = np.random.randint(1, 300, size=30) * 1e29 * u.erg
    flare_tbl['total_lc_time'] = len(lc['time']) * 120.0 * u.second

    save_path = csv_path.replace('.csv', '_flares.ecsv')
    flare_tbl.write(save_path, overwrite=True)


def generate_ffd(object, save_path, list_of_paths):
    monitoring_time = 0.0
    flare_energy = np.array([])

    for file_path in list_of_paths:
        tbl = ascii.read(file_path, guess=False, format='ecsv')
        monitoring_time += tbl['total_lc_time'][0]
        flare_energy = np.append(flare_energy, tbl['energy'].value)
    flare_energy.sort()

    cumulative_number = np.arange(len(flare_energy)) + 1
    flare_frequency = cumulative_number / monitoring_time

    fig, ax = plt.subplots()
    ax.plot(np.log10(flare_energy),
            flare_frequency,
            marker='o',
            color='darkcyan')
    ax.set(xlabel=r'Log$_{10}$ $E_{TESS}$',
           ylabel=r'Cumulative Number of Flares $>E_{TESS}$ Per Day',
           title='EFFD for {}'.format(object))
    fig.savefig('{}/{}_FFD.png'.format{save_path, object.replace(' ', '_')})
    plt.close(fig)

    # will need to do a linear regression to find the slope of the curve, BUT
    # we have to use the middle regime of the FFD without smaller/super
    # flares. So, we will have to make a function to take care of that.
