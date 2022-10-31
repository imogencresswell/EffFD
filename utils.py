import os
import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from scipy.optimize import curve_fit


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
        sector = result[0].mission[0][-2:]
        save_string = '{}/{}_{}'.format(save_path,
                                        object.replace(' ', '_'),
                                        sector)

        if os.path.isfile(save_string+'.csv'):
            print('Sector {} CSV exists for {}'.format(sector, object))
            continue

        # Leaving the older way of getting the LC from the pixel file
        # for a bit.
        #
        # pixelfile = search_result[window].download()
        # lc = pixelfile.to_lightcurve(aperture_mask='all')
        # lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

        lc = result.download()
        lc = lc.flatten(niters=filter_iter, sigma=filter_sig)
        lc.to_csv(save_string+'.csv', overwrite=True)

        plt.figure()
        lc.plot()
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

    # Toy data input for now
    flare_tbl = Table()
    flare_tbl['energy'] = np.random.choice(range(1, 5000),
                                           size=20,
                                           replace=False) * 1e29 * u.erg
    flare_tbl['total_lc_time'] = len(lc['time']) * 120.0 * u.second

    save_path = csv_path.replace('.csv', '_flares.ecsv')
    flare_tbl.write(save_path, overwrite=True)


def get_middle_ffd_regime(x, y, min_slope=-2.0, max_slope=-0.40):
    dx = np.diff(x, 1)
    dy = np.diff(y, 1)
    yfirst = dy / dx
    condition = np.where((yfirst > min_slope) & (yfirst < max_slope))[0]

    for idx, place in enumerate(condition):
        if place - condition[idx+1] <= 2:
            starting_idx = place + 1
            break
    ending_idx = condition[-1] - 1

    new_x = x[starting_idx:ending_idx]
    new_y = y[starting_idx:ending_idx]
    return new_x, new_y


def func_powerlaw(x, a, b):
    return a + b*x


def calculate_slope_powerlaw(x, y):
    solution = curve_fit(func_powerlaw, x, y, maxfev=2000)
    a = solution[0][0]
    b = solution[0][1]
    b_err = b/np.sqrt(len(x))
    return a, b, b_err


def generate_ffd(object, save_path, list_of_paths):
    monitoring_time = 0.0
    flare_energy = np.array([])

    for file_path in list_of_paths:
        tbl = ascii.read(file_path, guess=False, format='ecsv')
        monitoring_time += tbl['total_lc_time'][0]
        flare_energy = np.append(flare_energy, tbl['energy'].value)

    flare_energy.sort()
    log_energy = np.log10(flare_energy)

    monitoring_time *= tbl['total_lc_time'].unit
    monitoring_time = monitoring_time.to(u.day)

    cumulative_number = np.arange(len(flare_energy)) + 1
    flare_frequency = cumulative_number[::-1] / monitoring_time.value
    log_frequency = np.log10(flare_frequency)

    # linear regression to get slope
    m_ene, m_fre = get_middle_ffd_regime(log_energy, log_frequency)
    intercept, slope, slope_err = calculate_slope_powerlaw(m_ene, m_fre)
    # alpha = np.abs(slope - 1)

    fig, ax = plt.subplots()
    ax.plot(log_energy,
            flare_frequency,
            marker='o',
            color='skyblue')
    ax.plot(m_ene,
            10**func_powerlaw(m_ene, intercept, slope),
            color='black',
            label=r'Slope: $%.2f\pm%.2f$' % (slope, slope_err))
    ax.set(xlabel=r'Log$_{10}$ $E_{TESS}$ [%s]' % tbl['energy'].unit,
           ylabel=r'Cumulative Number of Flares $>E_{TESS}$ Per Day',
           title='EFFD for {}'.format(object),
           yscale='log')
    ax.legend()
    fig.savefig('{}/{}_FFD.png'.format(save_path, object.replace(' ', '_')))
    plt.close(fig)
