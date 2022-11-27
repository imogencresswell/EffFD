# EffFD: Effective Flare Frequency Distribution Software

## Concept
“EffFD” will be a software package that works with light curves (observations of stellar light at regular time intervals [e.g. 20 sec]) from the Transiting Exoplanet Survey Satellite (TESS; a mission that observed many stars across the entire sky over the last few years). In most cases, these light curves look like sin waves with scattered temporary rises, which indicate that a stellar flare happened on the star (similar to the solar flares you might hear about in the news). Our plan is to have “EffFD” automatically find those flares using simple de-trending methods, return a table of the flares (with properties like duration, energy, etc.), and create flare-frequency distributions (FFDs), which are graphs that determine how often a flare of a certain energy or larger would appear. The FFDs will then be compared between different types of stars (mostly based on size) to see if there is an underlying trend.

## Reason
Flare frequency distributions are useful for various groups in astronomy. On the stellar physics side, it helps us understand the rate at which stars expel energy. This can be expanded to studying the statistical differences between different types of stars. On the exoplanet studies side, flaring rates help determine if a stellar system could be habitable or help theorize the current state of the planets there.

On software, many astronomers write their own codes from scratch for every project, which leads to work being unreproducible and unstandardized. Occasionally, flare-finding software is released, but they increasingly tend to rely on complicated Bayesian statistics and post-data modeling (e.g. Gaussian fitting), which may not be representative of the actual data. Many groups refuse to use these programs and end up relying on by-eye estimations. While the ‘preferred ways’ are still heavily debated, we believe a program that uses understandable de-trending/flare-finding methods (or at least some of the simpler ones in the field) would work well for many purposes.

## Current To-do
- What kind of unit tests could we do for the functions that are mainly parsing through folders?

## Potential To-do
- We have a general way of selecting the middle regime of the FFD, but it would be better to test this against real data to see how it holds up and possibly figure out a more rigorous algorithm. This one is decently simple to explain, though, so if tests against real data come out fine, it could be something to stick with.
- ISSUE: search_lightcurve does not seem to work for all stars. AF Psc has a light curve in the MAST archive, but it doesn't return results, even if you remove the exptime and author/mission qualifiers. Not sure if there's anything we can do about this.
- Some people might want to use just the easy FFD creation aspect with more complex flare finding routines. We could add a functionality to feed in just the property tables (e.g. star_01_flares.ecsv) and have figures created.

## Notes
- Searching by Sector does not care if stars have temperatures or not, but it would be cool if it could auto-reject non-star observations. Not sure how though, using astroquery would increase run times a lot. Lightkurve might not actually have LCs for non-stars, so it wouldn't matter.

## Dependencies
Some users have had trouble building the Conda environment using the provided `environment.yml`. In that case, this code should be usable by first installing `astropy`, `scipy`, `lightkurve`, `pycodestyle`, `wget`, and `docopt`. The packages `numpy`, `matplotlib`, and `astroquery` are also required, but these should be installed as dependencies of the others.
