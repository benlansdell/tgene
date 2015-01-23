#!/usr/bin/env python

from io import *
import unittest

class TestIoFunctions(unittest.TestCase):

    def setUp(self):
        self.matrixarray1 = array([[1,2,3],[1,2,3]])
        self.matrixarray2 = array([])
        self.matrixarray3 = array([[[1,0],[1,0]],[[0,1],[1,0]]])
        self.matrixarray4 = array([1,2,3])
   
    def test_writeMatrix(self): #(ref, source, etype, coords, score, strand, frame, name):
        self.assertRaises(ValueError, writeMatrix, self.matrixarray3, 'test.mat')
        self.assertEqual(writeMatrix(self.matrixarray1, 'test1.mat'), True)
        self.assertEqual(writeMatrix(self.matrixarray2, 'test2.mat'), True)
        self.assertEqual(writeMatrix(self.matrixarray4, 'test4.mat'), True)      
    def test_readMatrix(self): #(ref, source, etype, coords, score, strand, frame, name):
        #self.assertRaises(ValueError, writeMatrix, self.matrixarray3, 'test.mat')
        self.assert_((readMatrix('test1.mat') == self.matrixarray1).all())
        self.assert_((readMatrix('test2.mat') == self.matrixarray2).all())
        self.assert_((readMatrix('test4.mat') == self.matrixarray4).all())
        

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestIoFunctions)
	unittest.TextTestRunner(verbosity=2).run(suite)