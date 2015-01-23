#!/usr/bin/env python

from maths import *
import unittest

class TestMathsFunctions(unittest.TestCase):

    def setUp(self):
        self.matrixarray = array([[1,2,3],[1,2,3]])
        self.nullarray = array([])
        self.vectorarray = array([1,2,3,4])
        self.null = []
        self.seq = range(10)
        self.baddist = [0.01, 0.1, 0.9]
        self.gooddist = [0.1, 0.1, 0.8]
        self.negdist = [-.1, 0.3, 0.8]
        self.kernelest1 = [1,2,3]
        self.kernelest2 = [-1,2,3]
        self.kernelest3 = range(10) + [100]

    def test_arrayDimension(self):
        self.assertEqual(arrayDimension(self.matrixarray), 2)
        self.assertEqual(arrayDimension(self.vectorarray), 1)
        self.assertEqual(arrayDimension(self.nullarray), 0)

    def test_chooseRand(self):
        self.assertRaises(ValueError, chooseRand, self.baddist)
        self.assertRaises(ValueError, chooseRand, self.negdist)
        for i in range(10):
            self.assert_(chooseRand(self.gooddist) in range(len(self.gooddist)))
            
    def test_normalCdf(self):
        self.assertEqual(normalCdf(-11, 0, 1), 0)
        self.assertEqual(normalCdf(0,0,20), 0.5)
        self.assertEqual(normalCdf(11, 0, 1), 1)
        self.assertRaises(ZeroDivisionError, normalCdf, 0, 0, 0)
        
    def test_discreteNormal(self): #(mu, sigma, xmin, xmax)
        self.assert_(round(sum(discreteNormal(10, 10, 0, 100)),10)== 1)
        self.assert_(round(sum(discreteNormal(10, 10, 20, 100)),10)== 1)
        self.assertEqual(len(discreteNormal(10,10,20,100)), 80)
        self.assertEqual(len(discreteNormal(10,10,10,10)), 0)

    def test_ctsNormal(self): #(mu, sigma, xmin, xmax)
        self.assertEqual(len(ctsNormal(10,10,20,100)), 80)
        self.assertEqual(len(ctsNormal(10,10,10,10)), 0)

    def test_continuousKernelDensity(self):
        self.assertRaises(ValueError, continuousKernelDensity, [], 1)
        self.assertEqual(len(continuousKernelDensity(self.kernelest1, 1)[0]), 4)
        self.assertEqual(len(continuousKernelDensity(self.kernelest2, 1)[0]), 5)
        
    def test_kernelDensity(self):
        self.assertRaises(ValueError, kernelDensity, [], 1)
        self.assertEqual(len(kernelDensity(self.kernelest1, 1)[0]), 4)
        self.assertEqual(len(kernelDensity(self.kernelest2, 1)[0]), 5)
        self.assertEqual(sum(kernelDensity(self.kernelest1, 1, 0.05, 100)[0]), 1)

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestMathsFunctions)
	unittest.TextTestRunner(verbosity=2).run(suite)


