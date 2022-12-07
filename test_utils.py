"""This tests the functions in utils.py

General goals:
- Check that functions require proper input types or raise errors for
  inputs that do not make sense
- Check that expected trends are returned given random data

"""
import sys
import os
import unittest
import numpy as np
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(src_path)
import utils as ut  # nopep8


class TestDataProcessor(unittest.TestCase):

    def setUp(self):
        # Generate fake FFD data
        self.y = np.unique(np.random.choice(range(1, 5000), size=100))
        self.y.sort()
        self.lo_y, self.lo_x = ut.get_log_freq(self.y, 80.0)
        self.re_x, self.re_y = ut.get_middle_ffd_regime(self.lo_x, self.lo_y)

        self.tics = np.array([['123', '1234'], [20, 30]])
        np.savetxt('./sector9999.csv', self.tics.T,
                   delimiter=',', fmt='%s', header='\n\n\n\n\n\n')

    def test_get_spectral_temp(self):
        # positive: check that a valid spectral type gives the right values
        self.assertTrue(ut.get_spectral_temp('M') == (2000, 3500))

        # negative: valid spectral type does not return random values
        self.assertFalse(ut.get_spectral_temp('M') != (2000, 3500))

        # asserts: Non-strings give TypeError.
        #          Non-recognized spectral types give ValueError.
        self.assertRaises(ValueError, ut.get_spectral_temp, 'banana')
        self.assertRaises(ValueError, ut.get_spectral_temp, 'Y')
        self.assertRaises(TypeError, ut.get_spectral_temp, 10)

    def test_save_sector(self):
        # positive/negative unncessary as function does file retrival

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.save_sector, [10, 10], './')
        self.assertRaises(TypeError, ut.save_sector, 10, './')
        self.assertRaises(TypeError, ut.save_sector, ['a', 'a'], 10)

    def test_get_sector_tics(self):
        tics_list = ut.get_sector_tics(['9999'], './')

        # positive:
        self.assertTrue(tics_list == ['123', '1234'])

        # negative:
        self.assertFalse(tics_list != ['123', '1234'])

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.get_sector_tics, [10, 10], './')
        self.assertRaises(TypeError, ut.get_sector_tics, 10, './')
        self.assertRaises(TypeError, ut.get_sector_tics, ['a', 'a'], 10)

    def test_build_names_from_sectors(self):
        tics_names = ut.build_names_from_sectors(['9999'], './')

        # positive:
        self.assertTrue(tics_names == ['TIC 123', 'TIC 1234'])

        # negative:
        self.assertFalse(tics_names != ['TIC 123', 'TIC 1234'])

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.build_names_from_sectors, [1], './')
        self.assertRaises(TypeError, ut.build_names_from_sectors, 10, './')
        self.assertRaises(TypeError, ut.build_names_from_sectors, ['a'], 10)

    def build_all_stars_table(self):
        # positive/negative unnecessary, as this is file generation
        # NOTE: This function will likely be removed or heavily modified,
        #       since this method takes too long and we have found an alt
        #       method to provide a list of TESS stars with temperatures

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.build_all_stars_table, [10, 10], './')
        self.assertRaises(TypeError, ut.build_all_stars_table, 10, './')
        self.assertRaises(TypeError, ut.build_all_stars_table, ['a', 'a'], 10)

    def test_save_raw_lc(self):
        # positive/negative unncessary - file generation mainly from lightkurve

        # assert: Check that incorrect types raise errors.
        #         Check that files that do not exist raise errors.
        self.assertRaises(TypeError, ut.save_raw_lc, 1, './', 2, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 'star', 5, 2, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 'star', './', 2.0, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 1, './', 2, 2)
        # The following gives a "ResourceWarning: unclosed <ssl.SSLSocket"
        # Need to look into this before leaving it.
        # self.assertRaises(FileNotFoundError, ut.save_raw_lc,
        #                   'ZZ top', './', 2, 2.0)

    def test_analyze_lc(self):
        # positive/negative unncessary - file generation

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.analyze_lc, 10)

    def test_get_middle_ffd_regime(self):
        # positive: Returned array is shorter than the input.
        self.assertTrue(len(self.re_x) <= len(self.lo_x) and
                        len(self.re_y) <= len(self.lo_y))

        # negative: returned array is not longer than the input one
        self.assertFalse(len(self.re_x) > len(self.lo_x) and
                         len(self.re_y) > len(self.lo_y))

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, [5, 4], [2, 2])
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 'banana', 'Y')
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 10.0, 20.0)

    def test_calculate_slope_powerlaw(self):
        # The first part of this test shows that the using this code with
        # published data returns results within error.
        def convert_energies(u_data, b_data, v_data, mon_time_days):
            """This function is only used to convert flares energies to the
            U band from B and V, as used in Lacy et al 1976. Standard procedure
            is done the slope using our code for these published results."""
            b2u = np.log10(1.2 * 10**b_data)
            v2u = np.log10(1.79 * 10**v_data)
            tot_flares = np.concatenate([u_data, b2u, v2u])
            tot_flares.sort()
            t_mon = mon_time_days/24.

            eng, freq = ut.get_log_freq(10**tot_flares, t_mon)
            m_eng, m_freq = ut.get_middle_ffd_regime(eng, freq)
            inte, slo, slo_e = ut.calculate_slope_powerlaw(m_eng, m_freq)
            return inte, slo, slo_e

        # YZ CMi and YY Gem data from Lacy et al 1976
        yz_v = np.array([29.62, 30.32, 30.22, 30.77,
                         31.54, 30.10, 29.05, 30.01,
                         30.65, 30.64, 30.28, 30.13,
                         30.68, 30.56, 30.03, 29.66,
                         29.89, 30.74, 30.51, 29.99,
                         30.03, 31.02, 31.35, 31.57,
                         29.88, 30.26, 29.57, 29.76,
                         30.19, 29.79, 30.16, 30.49,
                         30.67, 30.67, 31.75, 29.90,
                         29.83, 30.27, 29.81])
        yz_b = np.array([30.16, 28.37, 30.45, 31.03,
                         31.70, 30.75, 30.10, 30.25,
                         31.01, 31.08, 30.63, 30.65,
                         30.98, 30.68, 30.16, 29.93,
                         30.08, 30.83, 30.66, 29.93,
                         30.06, 31.25, 31.57, 31.87,
                         30.09, 30.60, 30.12, 29.94,
                         30.10, 29.26, 29.51, 30.64,
                         30.25, 29.70, 31.87, 29.59,
                         30.21, 30.15, 31.18, 30.17,
                         30.42, 29.75, 29.25, 29.87,
                         29.83, 29.61, 30.23, 30.26,
                         29.97, 30.18, 29.82, 29.42])
        yz_u = np.array([30.32, 29.96, 30.49, 31.00,
                         31.67, 30.93, 30.23, 30.33,
                         30.96, 31.06, 30.82, 30.44,
                         30.66, 30.76, 30.30, 30.02,
                         30.27, 30.96, 30.84, 30.30,
                         30.27, 31.34, 31.66, 32.01,
                         30.20, 30.66, 30.25, 30.02,
                         30.31, 29.61, 30.18, 30.88,
                         30.48, 30.48, 31.83, 29.96,
                         30.28, 29.94, 29.32, 29.77,
                         30.45, 29.91, 29.76, 30.34,
                         29.95, 30.33, 30.50, 30.49,
                         29.63, 28.69, 29.82])
        yz_int, yz_slo, yz_slo_e = convert_energies(yz_u, yz_b, yz_v, 76.)
        yy_v = np.array([33.77, 33.5,  32.98, 32.47, 32.32, 32.79])
        yy_b = np.array([33.96, 33.71, 32.81, 32.35, 32.48, 32.78])
        yy_u = np.array([34.09, 33.78, 32.44, 32.01, 32.44, 32.82])
        yy_int, yy_slo, yy_slo_e = convert_energies(yy_u, yy_b, yy_v, 120.)

        # Random FFD values
        intercept, slope, slope_err = ut.calculate_slope_powerlaw(self.re_x,
                                                                  self.re_y)

        # Positive Unit Tests
        # -------------------
        # Intercept exists and is a positive value.
        # The slope error is less than slope itself.
        # Note: Slope is always negative, so abs() is used to compare.
        # FFD slope should be between -0.4 and -5.0
        # The last two test that published values above are within error
        self.assertTrue(intercept)
        self.assertGreater(intercept, 0.0)
        self.assertGreater(abs(slope), slope_err)
        self.assertLessEqual(slope, -0.4)
        self.assertGreaterEqual(slope, -5.0)
        self.assertLess(np.abs(yz_slo + 0.71), yz_slo_e)
        self.assertLess(np.abs(yy_slo + 0.43), yy_slo_e)

        # Negative Unit Tests
        # -------------------
        # FFD slope is not beyond the bounds
        # intercept is not below 0
        self.assertFalse(slope > -0.4 and slope < -5.0)
        self.assertFalse(intercept < 0.0)

        # Error Assertions
        # ----------------
        # Check that incorrect types raise errors.
        # Check that empty lists do not pass.
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [5, 4], [2])
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 'ban', 'Y')
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 10.0, 20.0)
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [], [])
        self.assertRaises(ValueError, ut.calculate_slope_powerlaw,
                          np.array([]), np.array([]))

    def test_get_time_and_energy(self):
        # positive/negative unnecessary for now - reading from many files
        # will update once the code is switched to object-oriented format

        # assert:
        # Check that incorrect types raise errors.
        # Check that empty lists do not pass.
        self.assertRaises(TypeError, ut.get_time_and_energy, [-1, 1.0])
        self.assertRaises(ValueError, ut.get_time_and_energy, [])

    def test_get_log_freq(self):
        # positive: returned array is = to log10(original)
        self.assertTrue((self.lo_y == np.log10(self.y)).all())

        # negative: returned array is not different than log10(original)
        self.assertFalse((self.lo_y != np.log10(self.y)).all())

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.get_log_freq, np.array([0]), 'b')
        self.assertRaises(TypeError, ut.get_log_freq, 'b', 20.0)

    def generate_ffd(self):
        # positive/negative unnecessary as this outputs figures, which are
        # covered in the functional tests.

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.generate_ffd, 10, './', ['./'])
        self.assertRaises(TypeError, ut.generate_ffd, 'star', 10, ['./'])
        self.assertRaises(TypeError, ut.generate_ffd, 'star', './', 10)

    def tearDown(self):
        self.y = 0
        self.re_x = 0
        self.re_y = 0
        self.lo_x = 0
        self.lo_y = 0

        self.tics = 0
        os.remove('./sector9999.csv')


if __name__ == "__main__":
    unittest.main()
