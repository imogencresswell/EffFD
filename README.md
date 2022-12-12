# EffFD: Effective Flare Frequency Distribution Software

## Concept
“EffFD” will be a software package that works with light curves (observations of stellar light at regular time intervals [e.g. 20 sec]) from the Transiting Exoplanet Survey Satellite (TESS; a mission that observed many stars across the entire sky over the last few years). In most cases, these light curves look like sin waves with scattered temporary rises, which indicate that a stellar flare happened on the star (similar to the solar flares you might hear about in the news). Our plan is to have “EffFD” automatically find those flares using simple de-trending methods, return a table of the flares (with properties like duration, energy, etc.), and create flare-frequency distributions (FFDs), which are graphs that determine how often a flare of a certain energy or larger would appear. The FFDs will then be compared between different types of stars (mostly based on size) to see if there is an underlying trend.

## Reason
Flare frequency distributions are useful for various groups in astronomy. On the stellar physics side, it helps us understand the rate at which stars expel energy. This can be expanded to studying the statistical differences between different types of stars. On the exoplanet studies side, flaring rates help determine if a stellar system could be habitable or help theorize the current state of the planets there.

On software, many astronomers write their own codes from scratch for every project, which leads to work being unreproducible and unstandardized. Occasionally, flare-finding software is released, but they increasingly tend to rely on complicated Bayesian statistics and post-data modeling (e.g. Gaussian fitting), which may not be representative of the actual data. Many groups refuse to use these programs and end up relying on by-eye estimations. While the ‘preferred ways’ are still heavily debated, we believe a program that uses understandable de-trending/flare-finding methods (or at least some of the simpler ones in the field) would work well for many purposes.

## Usage 
To run EffFD on your local machine, first build the conda environment by typing 
`conda env create -f envrionment.yml`

Then type the command,
`python3 main.py ./config_file`,

replacing `./` with the path if your config file is stored in a different directory.

This will download TESS data, create light curves and FFDs for stars outlined in your config file 

### Dependencies
Some users have had trouble building the Conda environment using the provided `environment.yml`. In that case, this code should be usable by first installing `astropy`, `scipy`, `lightkurve`, `pycodestyle`, `wget`, and `docopt`. The packages `numpy`, `matplotlib`, and `astroquery` are also required, but these should be installed as dependencies of the others.

### Config file 
The options for the config file are as follows:
- `star_names`: This takes a list or a singular star name, or TIC number 
- `sectors`: TESS sector 
- `spectral_class`: Takes main sequence spectral classes O,B,A,F,G,K,M
- `teff_low`: Minumum temperature
- `teff_high`: Maximum temperature 
- `star_max`: The maximum number of stars to generate FFD's for

You can also change the lightkurve paramerters

Where the data is saved can be changed with the `out_dir`
The code will search for previous searches specified by the `search_dir` option
If these are left as default they will be created where the code is run  

When using the spectral class or temperature options, a .csv file containing TIC numbers and temperatures is downloaded. The TICs within the correct temperature range will then be analysed. This could be ~10,000 stars or more so the user should specify star_max which will choose a random selection from the given temperature range.



