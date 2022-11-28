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
        # positive:

        # negative:

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.save_sector, [10, 10], './')
        self.assertRaises(TypeError, ut.save_sector, 10, './')
        self.assertRaises(TypeError, ut.save_sector, ['a', 'a'], 10)

    def test_get_sector_tics(self):
        # positive:

        # negative:

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.get_sector_tics, [10, 10], './')
        self.assertRaises(TypeError, ut.get_sector_tics, 10, './')
        self.assertRaises(TypeError, ut.get_sector_tics, ['a', 'a'], 10)

    def test_build_names_from_sectors(self):
        # positive:

        # negative:

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.build_names_from_sectors, [1], './')
        self.assertRaises(TypeError, ut.build_names_from_sectors, 10, './')
        self.assertRaises(TypeError, ut.build_names_from_sectors, ['a'], 10)

    def build_all_stars_table(self):
        # positive:

        # negative:

        # assert: Check that incorrect types raise errors.
        self.assertRaises(TypeError, ut.build_all_stars_table, [10, 10], './')
        self.assertRaises(TypeError, ut.build_all_stars_table, 10, './')
        self.assertRaises(TypeError, ut.build_all_stars_table, ['a', 'a'], 10)

    def test_save_raw_lc(self):
        # positive:

        # negative:

        # assert: Check that incorrect types raise errors.
        #         Check that files that do not exist raise errors.
        self.assertRaises(TypeError, ut.save_raw_lc, 1, './', 2, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 'star', 5, 2, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 'star', './', 2.0, 2.0)
        self.assertRaises(TypeError, ut.save_raw_lc, 1, './', 2, 2)
        # self.assertRaises(FileNotFoundError, ut.save_raw_lc,
        #                   'ZZ top', './', 2, 2.0)

    def test_analyze_lc(self):
        # positive:

        # negative:

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
        intercept, slope, slope_err = ut.calculate_slope_powerlaw(self.re_x,
                                                                  self.re_y)

        # Positive Unit Tests
        # -------------------
        # Intercept exists and is a positive value.
        # The slope error is less than slope itself.
        # Note: Slope is always negative, so abs() is used to compare.
        # FFD slope should be between -0.4 and -5.0
        self.assertTrue(intercept)
        self.assertGreater(intercept, 0.0)
        self.assertGreater(abs(slope), slope_err)
        self.assertLessEqual(slope, -0.4)
        self.assertGreaterEqual(slope, -5.0)

        # Negative Unit Tests
        # -------------------
        # FFD slope is not beyond the bounds
        self.assertFalse(slope > -0.4 and slope < -5.0)

        # Error Assertions
        # ----------------
        # Check that incorrect types raise errors.
        # Check that empty lists do not pass.
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [5, 4], [2])
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 'ban', 'Y')
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 10.0, 20.0)
        self.assertRaises(ValueError, ut.calculate_slope_powerlaw, [], [])

    def test_get_time_and_energy(self):
        # positive:

        # negative:

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

    def generate_ffd(self):
        # positive:

        # negative:

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


if __name__ == "__main__":
    unittest.main()
