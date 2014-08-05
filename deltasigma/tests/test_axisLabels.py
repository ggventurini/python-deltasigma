import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises

class TestAxisLabels(unittest.TestCase):
    """Test function for axisLabels()"""
    def setUp(self):
        self.ran = np.arange(100)
    
    def test_axis_label_geneation_part1(self):
        """Test function for axisLabels() 1/3"""
        ss = ds.axisLabels(self.ran, incr=10)
        r = ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90']
        self.assertTrue(r == ss)
        
    def test_axis_label_generation_part2(self):
        """Test function for axisLabels() 2/3"""
        ss = ds.axisLabels(self.ran, incr=(15, 10))
        r = ['10', '25', '40', '55', '70', '85']
        self.assertTrue(r == ss)

    @raises(ValueError)
    def test_axis_label_incr_length(self):
        """Test function for axisLabels() 3/3"""
        ds.axisLabels([1., 2., 3., 4.], [1, 2, 3])

