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


def type_error_catch(var, vartype, inner_vartype=None, err_msg=None):
    if not isinstance(var, vartype):
        if err_msg is None:
            err_msg = '{} is not a {}'.format(var, vartype.__name__)
        raise TypeError(err_msg)

    elif vartype is list:
        if not var:
            raise ValueError('No passed lists should be empty.')

        for inner in var:
            if not isinstance(inner, inner_vartype):
                if err_msg is None:
                    err_msg = '{} is not a {}'.format(inner,
                                                      inner_vartype.__name__)
                raise TypeError(err_msg)


def get_spectral_temp(classification):
    """Returns upper and lower temperature limits in
    Kelvin of required spectral class

    Parameters
    ----------
    classification : str
        Stellar spectral type to get temperatures for

    Returns
    -------
    temperature_limits : tuple of two int
        Upper and lower temperature of given spectral class

    """
    type_error_catch(classification, str)

    if classification in 'L,T,Y'.split(','):
        raise ValueError('Brown dwarfs are not yet supported.')
    elif classification not in 'O,B,A,F,G,K,M'.split(','):
        raise ValueError('Please use spectral type in OBAFGKM.')
    else:
        # low temp limit, high temp limit, in K
        temp_dict = {'M': (2000, 3500),
                     'K': (3500, 5000),
                     'G': (5000, 6000),
                     'F': (6000, 7500),
                     'A': (7500, 10000),
                     'B': (10000, 30000),
                     'O': (30000, 60000)}
        return temp_dict[classification]


def save_sector(sector, search_path):
    """Function that retrieves data for each sector from the TESS website

    Parameters
    ----------
    sector : str
        TESS sector to search for

    search_path : str
        Path to search directory where previous queries are saved

    Returns
    -------
    None
        Saves the sector data as a CSV in the search folder

    """
    type_error_catch(sector, str)
    type_error_catch(search_path, str)

    save_string = search_path + 'sector{}.csv'.format(sector)
    if not os.path.isfile(save_string):
        try:
            print('Downloading Sector {} observation list.'.format(sector))
            url = "https://tess.mit.edu/wp-content/uploads/all_targets_S{}_v1.csv".format(sector.zfill(3))  # nopep8
            urllib.request.urlretrieve(url, save_string)
        except urllib.error.HTTPError:
            raise ValueError('Inputted URL could not be found.')


def get_sector_tics(sector_list, search_path):
    """For given TESS Sectors, gets the unique TESS identifiers (TIC) of all
    observed objects

    Parameters
    ----------
    sector_list : list of str
        List of TESS Sectors

    search_path : str
        Path to search directory where previous queries are saved

    Returns
    -------
    names_list : list of str
        TESS identifiers (TIC) for unique objects in given Sectors

    """
    type_error_catch(sector_list, list, str)
    type_error_catch(search_path, str)

    tic_array = np.array([])

    for sector in sector_list:
        curr_csv = np.genfromtxt(search_path + 'sector{}.csv'.format(sector),
                                 delimiter=',',
                                 skip_header=6)
        tic_array = np.append(tic_array, curr_csv[:, 0])

    tic_array = np.unique(tic_array)
    tic_list = tic_array.astype(int).astype(str).tolist()
    return tic_list


def build_names_from_sectors(sector_list, search_path):
    """Returns TESS names (e.g. 'TIC 1234') of objects observed in the
    listed Sectors

    Parameters
    ----------
    sector_list : list of str
        List TESS sectors to search through

    search_path: str
        Path to search directory where previous queries are saved

    Returns
    -------
    tess_names_list : list of str
        TESS names for all observed objects in given Sectors

    """
    type_error_catch(sector_list, list, str)
    type_error_catch(search_path, str)

    tics = get_sector_tics(sector_list, search_path)
    tess_names_list = ['TIC ' + name for name in tics]
    return tess_names_list


def build_all_stars_table(tic_list, search_path):
    type_error_catch(tic_list, list, str)
    type_error_catch(search_path, str)

    save_path = search_path + 'all_stars_table.csv'
    star_tbl = None
    v = Vizier(catalog="IV/39/tic82", row_limit=-1, timeout=300)

    for count, tic in enumerate(tic_list):
        print('\nTIC {} {}/{}.'.format(tic, count + 1, len(tic_list)))
        cat = v.query_constraints(TIC=tic)[0]

        if cat['S_G'] == 'STAR' and cat['Teff'] > 0.0:
            if star_tbl is None:
                star_tbl = cat
            else:
                star_tbl.add_row(cat[0])

    star_tbl.write(save_path)


def save_raw_lc(obj, save_path, filter_iter, filter_sig):
    """Uses lightkurve to retrieve the lightcurve from TESS
    object and saves the light curve image plus the raw data

    Parameters
    ----------
    obj : str
        TESS object (mostly stars)

    save_path : str
        Path to the save directory

    filter_iter: int
        Number of iterations that lightkurve smooths data over

    filter_sig : float
        Statistical sigma at which lightkurve cuts off data

    """
    type_error_catch(obj, str)
    type_error_catch(save_path, str)
    type_error_catch(filter_iter, int)
    type_error_catch(filter_sig, float)

    # SPOC == TESS data pipeline
    # Getting only the 120 second exposure light curves for consistency
    search_result = lk.search_lightcurve(obj, author='SPOC', exptime=120)
    if not search_result:
        raise FileNotFoundError('No results for {}.'.format(obj))

    for result in search_result:
        sector = result[0].mission[0][-2:]

        # save files named after star+sector, in the star's output directory
        save_string = '{}/{}_{}'.format(save_path,
                                        obj.replace(' ', '_'),
                                        sector)

        if os.path.isfile(save_string+'.csv'):
            print('Sector {} CSV exists for {}'.format(sector, obj))
            continue

        lc = result.download()
        lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

        # Save light curve CSV file
        lc.to_csv(save_string + '.csv', overwrite=True)

        # Saves light curve PNG plot
        plt.figure()
        lc.plot()
        plt.savefig(save_string + '.png')
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
    type_error_catch(csv_path, str)

    lc = ascii.read(csv_path, guess=False, format='csv')

    # # # # # # # # # # # # # # # # # #
    # FLARE FINDING METHOD GOES HERE. #
    # # # # # # # # # # # # # # # # # #

    # # # # # TOY DATA INPUT FOR NOW # # # # #
    yy_v = np.array([29.62, 30.32, 30.22, 30.77, 31.54, 30.10, 29.05, 30.01,
                     30.65, 30.64, 30.28, 30.13, 30.68, 30.56, 30.03, 29.66,
                     29.89, 30.74, 30.51, 29.99, 30.03, 31.02, 31.35, 31.57,
                     29.88, 30.26, 29.57, 29.76, 30.19, 29.79, 30.16, 30.49,
                     30.67, 30.67, 31.75, 29.90, 29.83, 30.27, 29.81])
    yy_b = np.array([30.16, 28.37, 30.45, 31.03, 31.70, 30.75, 30.10, 30.25,
                     31.01, 31.08, 30.63, 30.65, 30.98, 30.68, 30.16, 29.93,
                     30.08, 30.83, 30.66, 29.93, 30.06, 31.25, 31.57, 31.87,
                     30.09, 30.60, 30.12, 29.94, 30.10, 29.26, 29.51, 30.64,
                     30.25, 29.70, 31.87, 29.59, 30.21, 30.15, 31.18, 30.17,
                     30.42, 29.75, 29.25, 29.87, 29.83, 29.61, 30.23, 30.26,
                     29.97, 30.18, 29.82, 29.42, ])
    yy_u = np.array([30.32, 29.96, 30.49, 31.00, 31.67, 30.93, 30.23, 30.33,
                     30.96, 31.06, 30.82, 30.44, 30.66, 30.76, 30.30, 30.02,
                     30.27, 30.96, 30.84, 30.30, 30.27, 31.34, 31.66, 32.01,
                     30.20, 30.66, 30.25, 30.02, 30.31, 29.61, 30.18, 30.88,
                     30.48, 30.48, 31.83, 29.96, 30.28, 29.94, 29.32, 29.77,
                     30.45, 29.91, 29.76, 30.34, 29.95, 30.33, 30.50, 30.49,
                     29.63, 28.69, 29.82])
    # yy_v = np.array([33.77, 33.5 , 32.98, 32.47, 32.32, 32.79])
    # yy_b = np.array([33.96, 33.71, 32.81, 32.35, 32.48, 32.78])
    # yy_u = np.array([34.09, 33.78, 32.44, 32.01, 32.44, 32.82])
    yy_b2u = np.log10(1.2 * 10**yy_b)
    yy_v2u = np.log10(1.79 * 10**yy_v)
    concatarray = np.concatenate([yy_u, yy_b2u, yy_v2u])

    flare_tbl = Table()
    # any u.unit converts the data into the correct astrophysical unit
    flare_tbl['energy'] = 10**concatarray * u.erg
    # flare_tbl['energy'] = np.random.choice(range(1, 5000),
    #                                        size=20,
    #                                        replace=False) * 1e29 * u.erg
    flare_tbl['total_lc_time'] = len(lc['time']) * 120.0 * u.second
    # # # # # END TOY DATA # # # # #

    save_path = csv_path.replace('.csv', '_flares.ecsv')
    flare_tbl.write(save_path, overwrite=True)


def get_middle_ffd_regime(x, y, min_slope=-5.0, max_slope=-1.0):
    """Finds the location of the middle regime of the flare frequency diagram
    (FFD) using the min and max slope of where most middle regimes lie

    Parameters
    ----------
    x : numpy array
        Energy array

    y : numpy array
        Cumulative frequency array

    min_slope : float
        minimum value of the slope for the condition

    max_slope : float
        maximum value of the slope for the condition

    Returns
    -------
    new_x : numpy array
        Energy array within middle regime
    new_y : numpy array
        Cumulative frequency array within middle regime

    """
    type_error_catch(x, np.ndarray)
    type_error_catch(y, np.ndarray)
    type_error_catch(min_slope, float)
    type_error_catch(max_slope, float)

    dx = np.diff(x, 1)
    dy = np.diff(y, 1)
    yfirst = dy / dx
    # finds points that satisfy the slope condition for the middle regime
    cond = np.where((yfirst > min_slope) & (yfirst < max_slope))[0]
    # in the first regime some points will satisfy the above condition
    # this loop is checking that the next two points also satisfy the
    # condition so that we only get data points in the middle regime and
    # not rogue points from first regime

    try:
        for count, idx in enumerate(cond):  # note: cond is an array of idxs
            if idx-cond[count+1] <= 2 and cond[count+3]-cond[count+2] <= 2:
                starting_idx = idx
                break
        ending_idx = cond[-1]
    except IndexError:
        starting_idx = cond[0]
        ending_idx = cond[-1]

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
    a : float
        intercept
    b : float
        slope
    b_err : float
        error on slope

    """
    type_error_catch(x, np.ndarray)
    type_error_catch(y, np.ndarray)

    # get the fit to func_powerlaw() using dataset (x, y)
    solution = curve_fit(func_powerlaw, x, y, maxfev=2000)

    a = solution[0][0]  # intercept
    b = solution[0][1]  # slope
    b_err = np.abs(b / np.sqrt(len(x)))  # slope_err

    return a, b, b_err


def get_time_and_energy(paths):
    """Takes a list of data directory paths and finds the
    total time the object was observed and the energy of
    each flare
    Note: Some TESS sectors overlap therefore the same object
    might be in multiple sectors

    Parameters
    ----------
    paths: list of str
        Path to object data.
        If the object was observed in multiple TESS sectors, there will be
        multiple paths, otherwise it will be a single filepath.

    Returns
    -------
    time : float
        Total time TESS observed the object, in units of days

    flare_eng : numpy array
        Flare energies, sorted by total size

    e_unit : astropy unit
        Unit used for flare energies, to be used in plotting

    """
    type_error_catch(paths, list, str)

    time = 0.0 * u.day
    flare_eng = np.array([])

    for file_path in paths:
        try:
            tbl = ascii.read(file_path, guess=False, format='ecsv')

            time += tbl['total_lc_time'][0] * \
                (1.0 * tbl['total_lc_time'].unit).to(u.day)

            flare_eng = np.append(flare_eng, tbl['energy'].value)

        except FileNotFoundError:
            print('Flare filepath ' + file_path + ' not found.')
            continue

    flare_eng.sort()  # Sort flares by size
    return time.value, flare_eng, tbl['energy'].unit


def get_log_freq(flare_eng, tot_time):
    """Takes the flare energy array and the time it was observed
    over and returns the cumulative frequency for each energy

    Parameters
    ----------
    flare_eng : numpy array
       array of flare energy

    tot_time : float
        total time TESS observed the object

    Returns
    -------
    energy : numpy array
        log10 flare energies

    frequency : numpy array
        log10 of the cumulative frequency

    """
    type_error_catch(flare_eng, np.ndarray)
    type_error_catch(tot_time, float)

    energy = np.log10(flare_eng)

    cumulative_count = np.arange(len(energy)) + 1
    flare_frequency = cumulative_count[::-1] / tot_time

    frequency = np.log10(flare_frequency)

    return energy, frequency


def generate_ffd(obj, save_path, list_of_paths):
    """This function generates and saves the flare freqeucny diagram (FFD)

    Parameters
    ----------
    obj : str
       Name of object (mainly stars)

    save_path: str
        Path to save directory

    list_of_paths: list of str
        Path to object data.
        If the object was observed in multiple TESS sectors, there will be
        multiple paths, otherwise it will be a single filepath.

    """
    type_error_catch(obj, str)
    type_error_catch(save_path, str)
    type_error_catch(list_of_paths, list, str)

    monitoring_time, flare_energy, e_unit = get_time_and_energy(list_of_paths)

    # # # # # THIS IS FOR THE TOY DATA ONLY - TBREMOVED # # # # #
    flare_energy = np.unique(flare_energy)
    # # # # # END OF TOY DATA PART # # # # #

    log_energy, log_frequency = get_log_freq(flare_energy, monitoring_time)

    # Linear regression to get slope
    m_ene, m_fre = get_middle_ffd_regime(log_energy, log_frequency)
    intercept, slope, slope_err = calculate_slope_powerlaw(m_ene, m_fre)
    # alpha = np.abs(slope - 1)

    # Saves FFD figure
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
           title='EFFD for {}'.format(obj),
           yscale='log')
    ax.legend()
    fig.savefig('{}/{}_FFD.png'.format(save_path, obj.replace(' ', '_')))
    plt.close(fig)
