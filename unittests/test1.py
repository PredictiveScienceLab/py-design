"""
Unittests for best.design

Author:
    Ilias Bilionis

Date:
    8/31/2013
"""


import best.design
import numpy as np
import unittest


class TestDesign(unittest.TestCase):

    def test_latinize(self):
        x = np.random.rand(1000, 2)
        x_lhc = best.design.latinize(x)


if __name__ == '__main__':
    unittest.main()
