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
        # Generate fake FFD data, as is done in utils.py
        self.y = np.unique(np.random.choice(range(1, 5000), size=100))
        self.y.sort()
        self.lo_y, self.lo_x = ut.get_log_freq(self.y, 80.0)
        self.re_x, self.re_y = ut.get_middle_ffd_regime(self.lo_x, self.lo_y)

    def test_get_spectral_temp(self):

        # positive: check that a valid spectral type gives the right values
        self.assertTrue(ut.get_spectral_temp('M') == (2000, 3500))

        # negative: valid spectral type does not return random values
        self.assertFalse(ut.get_spectral_temp('M') != (2000, 3500))

        # asserts: non-strings give TypeError, non-recognized spectral
        #          types give ValueError
        self.assertRaises(ValueError, ut.get_spectral_temp, 'banana')
        self.assertRaises(ValueError, ut.get_spectral_temp, 'Y')
        self.assertRaises(TypeError, ut.get_spectral_temp, 10)

    def test_get_middle_ffd_regime(self):

        # pos: returned array is shorter than the input
        self.assertTrue(len(self.re_x) <= len(self.lo_x) and
                        len(self.re_y) <= len(self.lo_y))

        # neg: returned array is not longer than the input one
        self.assertFalse(len(self.re_x) > len(self.lo_x) and
                         len(self.re_y) > len(self.lo_y))

        # asserts: checking that non-numpy-arrays don't pass
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, [5, 4], [2, 2])
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 'banana', 'Y')
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 10.0, 20.0)

    def test_calculate_slope_powerlaw(self):
        intercept, slope, slope_err = ut.calculate_slope_powerlaw(self.re_x,
                                                                  self.re_y)

        # positive: the error of the slope is less than the slope itself
        # NOTE: slope is always negative (error always pos), so abs() is used
        self.assertTrue(abs(slope) > slope_err)

        # positive: FFD slope should be between -0.4 and -2.0
        self.assertTrue(slope <= -0.4 and slope >= -2.0)

        # neg: FFD slope is not beyond the bounds
        self.assertFalse(slope > -0.4 and slope < -2.0)

        # asserts: TypeError when incorrect types are passed
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [5, 4], [2])
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 'ban', 'Y')
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 10.0, 20.0)

    def test_get_log_freq(self):
        # pos: returned array is = to log10(original)
        self.assertTrue((self.lo_y == np.log10(self.y)).all())

        # neg: returned array is not different than log10(original)
        self.assertFalse((self.lo_y != np.log10(self.y)).all())

    def test_get_time_and_energy(self):
        # asserts: type check and can't input empty list
        self.assertRaises(ValueError, ut.get_time_and_energy, [])
        self.assertRaises(TypeError, ut.get_time_and_energy, [-1, 1.0])

    def tearDown(self):
        self.y = 0
        self.re_x = 0
        self.re_y = 0
        self.lo_x = 0
        self.lo_y = 0


if __name__ == "__main__":
    unittest.main()
