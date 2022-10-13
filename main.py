import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt

pixelfile = lk.search_targetpixelfile("AU Mic")[2].download()
lc = pixelfile.to_lightcurve(aperture_mask='all').flatten()
lc.plot()
