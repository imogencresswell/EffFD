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


def group_by_missing(seq):
    """Takes input array and groups
    consecutive numbers
    Parameters
    ----------
    seq : list
        Data to be grouped

    Returns
    -------
    grouped : list
        List of list where each list
        is a group of consequative
        numbers
    """
    if not seq:
        return seq
    grouped = [[seq[0]]]
    for x in seq[1:]:
        if x == grouped[-1][-1] + 1:
            grouped[-1].append(x)
        else:
            grouped.append([x])
    return grouped


def analyze_lc(csv_path):
    """Takes the light curve data and finds flares
    using a 3 sigma method. If flares are found
    in the data the start time, end time, duration
    max flux, and time of the max flux are added
    to a table and saved.
    Parameters
    ----------
    csv_path : str
        Path to csv file for the sector data

    Returns
    -------
    _flares.ecsv file containing flare data
    """

    type_error_catch(csv_path, str)

    lc = ascii.read(csv_path, guess=False, format='csv')
    # Flare finding method
    criteria = 1 + 3*np.std(lc['flux'])
    criteria_index = np.where(lc['flux'] > criteria)[0]
    grouped_criteria = group_by_missing(criteria_index.tolist())

    flare_index = []
    for group in grouped_criteria:
        if len(group) >= 3:
            flare_index.append(group)

    if len(flare_index) == 0:
        print('No flares found in this sector data')

    else:

        flare_table = Table(names=['start_time', 'end_time', 'duration',
                                   'max_flux', 'max_flux_time', 'fluence'])

        for counts, flare in enumerate(flare_index):
            flare_flux = lc['flux'][flare] - 1  # quiscent=1; flare-only flux
            flare_time = lc['time'][flare]

            # Flare properties of interest
            t_start = lc['time'][flare[0]]
            t_end = lc['time'][flare[-1]]
            duration = t_start - t_end
            flux_max = np.max(flare_flux)
            t_flux_max = flare_time[(flare_flux == flux_max)]
            fluence = np.sum(flare_flux)

            flare_table.add_row([t_start, t_end, duration,
                                 flux_max, t_flux_max, fluence])

        # Save total light curve monitoring time for FFD statistics
        flare_table['total_lc_time'] = len(lc['time']) * 120.0 * u.second

        save_path = csv_path.replace('.csv', '_flares.ecsv')
        flare_table.write(save_path, overwrite=True)

        print(str(len(flare_table['fluence'])) +
              ' flares were found in this sector data!')


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

            flare_eng = np.append(flare_eng, tbl['fluence'].value)

        except FileNotFoundError:
            print('Flare filepath ' + file_path + ' not found.')
            continue

    flare_eng.sort()
    return time.value, flare_eng, tbl['fluence'].unit


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

    log_energy, log_frequency = get_log_freq(flare_energy, monitoring_time)

    # Linear regression to get slope
    try:
        m_ene, m_fre = get_middle_ffd_regime(log_energy, log_frequency)
        intercept, slope, slope_err = calculate_slope_powerlaw(m_ene, m_fre)
    except Exception:
        m_ene = []

    # alpha is used in some papers, but we don't need it for now
    # alpha = np.abs(slope - 1)

    # Saves FFD figure
    fig, ax = plt.subplots()
    ax.plot(log_energy,
            10**log_frequency,
            marker='o',
            color='skyblue')
    if len(m_ene) != 0:
        ax.plot(m_ene,
                10**func_powerlaw(m_ene, intercept, slope),
                color='black',
                label=r'Slope: $%.2f\pm%.2f$' % (slope, slope_err))
        ax.legend()
    ax.set(  # xlabel=r'Log$_{10}$ $E_{TESS}$ [%s]' % e_unit,
           xlabel=r'Log$_{10}$ TESS Fluence',
           ylabel=r'Cumulative Number of Flares $>E_{TESS}$ Per Day',
           title='EFFD for {}'.format(obj),
           yscale='log')

    fig.savefig('{}/{}_FFD.png'.format(save_path, obj.replace(' ', '_')))
    plt.close(fig)
