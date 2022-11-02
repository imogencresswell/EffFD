import os
import urllib.request
import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from scipy.optimize import curve_fit
from astroquery.vizier import Vizier


def get_spectral_temp(classification):
    """Returns upper and lower temperature limits in
    Kelvin of required spectral class
    Parameters
    ----------
    classification : str
        Star classification

    Returns
    -------
    Tuple containing upper and lower temperature of given
    spectral class

    """
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


def save_sector_list(sector, search_path):
    """Function that retrieves data for each sector from
    TESS website
    Parameters
    ----------
    sector : str
        TESS sector to search for

    search_path: str
        path to search directory where
        previous queries are saved

    Returns
    -------
    Saves the sector data in the search folder

    """
    save_string = search_path+'sector{}.csv'.format(sector)
    if not os.path.isfile(save_string):
        print('Downloading Sector {} observation list.'.format(sector))
        url = "https://tess.mit.edu/wp-content/uploads/all_targets_S{}_v1.csv".format(sector.zfill(3))  # nopep8
        urllib.request.urlretrieve(url, save_string)


def build_names_from_sectors(sector_list, search_path):
    """Returns names of stars that TESS observed in given
    sector
    Parameters
    ----------
    sector_list : list of strings
        List TESS sectors to search through

    search_path: str
        Path to search directory where
        previous queries are saved

    Returns
    -------
    A list of star names found in the specified sectors

    """
    names_array = np.array([])
    for sector in sector_list:
        curr_csv = np.genfromtxt(search_path+'sector{}.csv'.format(sector),
                                 delimiter=',',
                                 skip_header=6)
        names_array = np.append(names_array, curr_csv[:, 0])
    names_list = names_array.astype(int).astype(str).tolist()
    tess_names_list = ['TIC '+name for name in names_list]
    return tess_names_list


def call_tess_catalog_names(search_path, Tmin, Tmax):
    """Function that uses Vizier query to seach through
    TESS data for stars within a given temperature range
    Note: Currently takes far too long, alternative options
    being explored
    Parameters
    ----------
    Tmin : int
        Minimum temperature

    Tmax : int
        Maximum temperature

    search_path: string
        path to search directory where
        previous queries are saved

    Returns
    -------
    Saves the data for stars within the given
    temperature range

    """
    T_str = '{}..{}'.format(Tmin, Tmax)
    save_path = search_path + T_str + '.csv'

    try:
        ascii.read(save_path, guess=False, format='csv')
    except FileNotFoundError:
        v = Vizier(columns=['TIC', '_RAJ2000', '_DEJ2000'],
                   catalog="IV/39/tic82",
                   row_limit=-1,
                   timeout=300)
        cat = v.query_constraints(Teff=T_str)[0]

        cat.write(save_path)


def save_raw_lc(object, save_path, filter_iter, filter_sig):
    """Uses lightkurve to retrieve the lightcurve from TESS
    object and saves the lightcurves plus the raw data
    Parameters
    ----------
    object : str
        TESS object

    save_path : str
        path to the save directory

    filter_iter: int
        number of iterations that lightkurve
        smooths data over

    filter_sigma: float
        statistical sigma at which lightkurve
        cuts off data

    Returns
    -------
    Saves the lightcurve in a directory named after the object
    and its sector inside the save directory

    """

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

        lc = result.download()
        lc = lc.flatten(niters=filter_iter, sigma=filter_sig)
        lc.to_csv(save_string+'.csv', overwrite=True)

        plt.figure()
        lc.plot()
        plt.savefig(save_string+'.png')
        plt.close()


def analyze_lc(csv_path):
    """Takes the light curve data, finds flares in the
    light curve and saves this data to an ecsv file
    Note: currently this creates toy data for each light curve
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
    #any u.unit converts the data into the correct
    #astrophysical unit
    flare_tbl['energy'] = np.random.choice(range(1, 5000),
                                           size=20,
                                           replace=False) * 1e29 * u.erg
    flare_tbl['total_lc_time'] = len(lc['time']) * 120.0 * u.second

    save_path = csv_path.replace('.csv', '_flares.ecsv')
    flare_tbl.write(save_path, overwrite=True)


def get_middle_ffd_regime(x, y, min_slope=-2.0, max_slope=-0.40):
    """Finds the location of the middle regime of the FFD using
    the min and max slope of where most middle regimes lie
    Parameters
    ----------
    x : numpy array
        Energy array

    y : numpy array
        Cumulative frequency array

    min_slope: float
        minimum value of the slope for the condition

    max_slope: float
        maximum value of the slope for the condition

    Returns
    -------
    new_x : Energy array within middle regime
    new_y : Cumulative frequency array within middle regime

    """
    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise TypeError('x and y must be numpy array')

    dx = np.diff(x, 1)
    dy = np.diff(y, 1)
    yfirst = dy / dx
    # finds points that satisfy the slope condition for the middle regime
    cond = np.where((yfirst > min_slope) & (yfirst < max_slope))[0]
    # in the first regime some points will satisfy the above condition
    # this loop is checking that the next two points also satisfy the
    # condition so that we only get data points in the middle regime and
    # not rogue points from first regime
    for idx, pla in enumerate(cond):
        if pla-cond[idx+1] <= 2 and cond[idx+3]-cond[idx+2] <= 2:
            starting_idx = pla + 1
            break
    ending_idx = cond[-1] - 1

    new_x = x[starting_idx:ending_idx]
    new_y = y[starting_idx:ending_idx]
    return new_x, new_y


def func_powerlaw(x, a, b):
    return a + b*x


def calculate_slope_powerlaw(x, y):
    """Find the slope of powerlaw
    Parameters
    ----------
    x : numpy array
        x array

    y : numpy array
        y array

    Returns
    -------
    a : intercept
    b : slope
    b_err : error on slope

    """
    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise TypeError('x and y must be numpy array')

    solution = curve_fit(func_powerlaw, x, y, maxfev=2000)
    a = solution[0][0]
    b = solution[0][1]
    b_err = b/np.sqrt(len(x))
    return a, b, np.abs(b_err)


def get_time_and_energy(paths):
    """Takes a list of data directory paths and finds the
    total time the object was observed and the energy of
    each flare
    Note: Some TESS sectors overlap therefore the same object
    might be in multiple sectors
    Parameters
    ----------
    paths : list of strings
        Path to object data, if the object
        was observed in multiple TESS sectors
        there will be multiple paths, otherwise
        it will be a single filepath

    Returns
    -------
    time : total time TESS observed the object
    flare_eng : list of energies of the flares

    """
    if not paths:
        raise ValueError('List of paths is empty')

    time = 0.0 * u.day
    flare_eng = np.array([])
    for file_path in paths:
        if type(file_path) is not str:
            raise TypeError('object path must be a string')
        try:
            tbl = ascii.read(file_path, guess=False, format='ecsv')
            time += tbl['total_lc_time'][0] * \
                (1.0 * tbl['total_lc_time'].unit).to(u.day)
            flare_eng = np.append(flare_eng, tbl['energy'].value)
        except FileNotFoundError:
            print('Flare filepath '+file_path+' not found')
            continue
    flare_eng.sort()
    return time, flare_eng, tbl['energy'].unit


def get_log_freq(flare_eng, tot_time):
    """Takes the flare energy array and the time it was observed
    over and returns the cumulative frequency for each energy
    Parameters
    ----------
    flare_eng : array
       array of flare energy

    tot_time: total time TESS observed the object

    Returns
    -------
    time : list of times that a flare occurred
    flare_eng : list of energies of the flares

    """
    energy = np.log10(flare_eng)
    cumulative_count = np.arange(len(energy)) + 1
    flare_frequency = cumulative_count[::-1] / tot_time
    frequency = np.log10(flare_frequency)
    return energy, frequency


def generate_ffd(object, save_path, list_of_paths):
    """This function generates and saves the FFD
    Parameters
    ----------
    object : str
       Name of object

    save_path: str
        Path to save directory

    list_of_paths: list of str
        Path to object data, if the object
        was observed in multiple TESS sectors
        there will be multiple paths, otherwise
        it will be a single filepath

    Returns
    -------
    Figure of FFD for the object

    """
    monitoring_time, flare_energy, e_unit = get_time_and_energy(list_of_paths)

    # THIS IS FOR THE TOY DATA ONLY - TO BE REMOVED
    flare_energy = np.unique(flare_energy)
    # END OF TOY DATA PART

    log_energy, log_frequency = get_log_freq(flare_energy,
                                             monitoring_time.value)

    # linear regression to get slope
    m_ene, m_fre = get_middle_ffd_regime(log_energy, log_frequency)
    intercept, slope, slope_err = calculate_slope_powerlaw(m_ene, m_fre)
    # alpha = np.abs(slope - 1)

    fig, ax = plt.subplots()
    ax.plot(log_energy,
            10**log_frequency,
            marker='o',
            color='skyblue')
    ax.plot(m_ene,
            10**func_powerlaw(m_ene, intercept, slope),
            color='black',
            label=r'Slope: $%.2f\pm%.2f$' % (slope, slope_err))
    ax.set(xlabel=r'Log$_{10}$ $E_{TESS}$ [%s]' % e_unit,
           ylabel=r'Cumulative Number of Flares $>E_{TESS}$ Per Day',
           title='EFFD for {}'.format(object),
           yscale='log')
    ax.legend()
    fig.savefig('{}/{}_FFD.png'.format(save_path, object.replace(' ', '_')))
    plt.close(fig)
