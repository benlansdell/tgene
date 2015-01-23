#!/usr/bin/env python

from feature import *
import unittest

class TestFeatureFunctions(unittest.TestCase):

    def setUp(self):
        self.testFeature = Feature('ref', 'source', ['exon'], [[0,10]], [10], '+', None, 'Gene 1')

    def test_Feature(self): #(ref, source, etype, coords, score, strand, frame, name):
        self.assertRaises(
        #self.assertRaises(ValueError, SequenceDict, [0,1,2], 3)
        

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestFeatureFunctions)
	unittest.TextTestRunner(verbosity=2).run(suite)