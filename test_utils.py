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
        self.assertTrue(self.re_y.max() < -0.4 and
                        self.re_y.max() > -2.0)
        self.assertFalse(self.re_y.max() > -0.4 and
                         self.re_y.max() < -2.0)
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, [5,4], [2,2])
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 'banana', 'Y')
        self.assertRaises(TypeError, ut.get_middle_ffd_regime, 10.0, 20.0)

    def test_calculate_slope_powerlaw(self):
        intercept, slope, slope_err = calculate_slope_powerlaw(self.re_x,
                                                               self.re_y)
        self.assertTrue(slope > slope_err)
        self.assertTrue(slope =< -0.4 and slope > -2.0)
        self.assertFalse(slope > -0.4 and slope < -2.0)
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, [5,4], [2])
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 'bana', 'Y')
        self.assertRaises(TypeError, ut.calculate_slope_powerlaw, 10.0, 20.0)

    def setUp(self):
        self.x = np.linspace(0, 10, 100)
        self.y = np.random.uniform(low=-5, high=-0.1, size=(100,))
        self.re_x, self.re_y = ut.get_middle_ffd_regime(self.x, self.y)

    def tearDown(self):
        self.x = 0
        self.y = 0
        self.re_x = 0
        self.re_y = 0


if __name__ == "__main__":
    unittest.main()
