"""This tests the functions in utils.py to make sure they can handle different
situations where inputs cannot be found.

"""
import sys
import os
import unittest
import numpy as np
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(src_path)
import utils as ut  # nopep8


class TestDataProcessor(unittest.TestCase):

    def test_get_spectral_temp(self):
        self.assertTrue(ut.get_spectral_temp('M') == (2000, 3500))
        self.assertFalse(ut.get_spectral_temp('M') != (2000, 3500))
        self.assertRaises(ValueError, ut.get_spectral_temp, 'banana')
        self.assertRaises(ValueError, ut.get_spectral_temp, 'Y')
        self.assertRaises(ValueError, ut.get_spectral_temp, 10)

    def test_get_middle_ffd_regime(self):
        self.assertTrue(len(self.re_x) <= len(self.x) and
                        len(self.re_y) <= len(self.y))
        self.assertFalse(len(self.re_x) > len(self.x) and
                         len(self.re_y) > len(self.y))
        # self.assertRaises(TypeError, ut.get_middle_ffd_regime, [5,4], [2,2])
        # self.assertRaises(TypeError, ut.get_middle_ffd_regime, 'banana', 'Y')
        # self.assertRaises(TypeError, ut.get_middle_ffd_regime, 10.0, 20.0)

    def test_calculate_slope_powerlaw(self):
        intercept, slope, slope_err = ut.calculate_slope_powerlaw(self.re_x,
                                                                  self.re_y)
        # self.assertTrue(slope > slope_err)
        self.assertTrue(slope <= -0.4 and slope >= -2.0)
        self.assertFalse(slope > -0.4 and slope < -2.0)
        # self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [5,4], [2])
        # self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 'ban', 'Y')
        # self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 10.0, 20.0)

    def setUp(self):

        # Generate fake FFD data, as is done in utils.py
        self.y = np.unique(np.random.choice(range(1, 5000),
                                            size=100,
                                            replace=False))
        self.y.sort()
        self.x = np.arange(len(self.y))[::-1] + 1

        # Get only the middle section of the FFD
        self.re_x, self.re_y = ut.get_middle_ffd_regime(np.log10(self.x),
                                                        np.log10(self.y))

    def tearDown(self):
        self.x = 0
        self.y = 0
        self.re_x = 0
        self.re_y = 0


if __name__ == "__main__":
    unittest.main()
