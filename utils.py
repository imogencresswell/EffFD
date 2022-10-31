import lightkurve as lk
import matplotlib.pyplot as plt


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

    # Note SPOC == TESS data pipeline
    search_result = lk.search_targetpixelfile(object, author='SPOC')
    if not search_result:
        raise FileNotFoundError('No results for {}.'.format(object))

    for window in range(len(search_result)):
        pixelfile = search_result[window].download()
        lc = pixelfile.to_lightcurve(aperture_mask='all')
        lc = lc.flatten(niters=filter_iter, sigma=filter_sig)

        plt.figure()
        lc.plot()

        save_string = '{}/{}_{}'.format(save_path,
                                        object.replace(' ', '_'),
                                        window)
        lc.to_csv(save_string+'.csv', overwrite=True)
        plt.savefig(save_string+'.png')
