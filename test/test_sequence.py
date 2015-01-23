#!/usr/bin/env python

from sequence import *
from sequence import _decimal, _basefour, _genHigherOrder, _higherOrderDict
import unittest

class TestSequenceFunctions(unittest.TestCase):

    def test_seqToInt(self):
        self.assertEqual(seqToInt('ATGCGT'), [0,3,2,1,2,3])
        self.assertRaises(ValueError, seqToInt, 'ATGCB')
        self.assertRaises(ValueError, seqToInt, 'aTCG')
        
    def test_intToSeq(self):
        self.assertEqual(intToSeq([0,1,2,3]), 'ACGT')
        self.assertRaises(ValueError, intToSeq, [0,1,2,3,4])

    def test_translate(self):
        self.assertEqual(translate('ATGATGATG'), 'MMM')
        self.assertRaises(ValueError, translate, 'ATGATM')
        self.assertRaises(ValueError, translate, 'ATGATGAT')
        self.assertEqual(translate(''), '')
        self.assertEqual(translate([0,3,2]), 'M')
        self.assertEqual(translate('TAA'), '*')

    def test_reverse(self):
        self.assertEqual(reverse('ATG'), 'GTA')
        self.assertEqual(reverse([0,1,2]), [2,1,0])

    def test_reverseComlement(self):
        self.assertEqual(reverseComplement('ATG'), 'CAT')
        self.assertEqual(reverseComplement([0,1,2]), [1,2,3])
        self.assertRaises(ValueError, reverseComplement, 'ATGF')

    def test_replaceString(self):
        self.assertEqual(len(replaceString(range(10), range(5), 5)), 10)
        self.assertEqual(len(replaceString(range(10), range(5), 50)), 10)
        self.assertEqual(replaceString(range(10), range(5),5)[5], 0)
        self.assertEqual(replaceString(range(10), range(5),5)[4], 4)
        self.assertEqual(replaceString(range(10), range(5),5)[9], 4)

    def test_decimal(self):
        self.assertRaises(TypeError, _decimal, 'A')
        self.assertEqual(_decimal([0,1,2,3]), 228)
        self.assertRaises(ValueError, _decimal, [4])

    def test_basefour(self):#Deprecated
        self.assertEqual(_basefour(228, 3), [0,1,2,3])

    def test_genHigherOrder(self): #seq, order):
        self.assertRaises(ValueError, _genHigherOrder, [0,1,2], 3)
        self.assertRaises(TypeError, _genHigherOrder, 'ATGCF', 2)
        self.assertRaises(ValueError, _genHigherOrder, [0,4], 1)
        self.assertEqual(_genHigherOrder([0,1,2,3], 3), [None, None, None, 228])
        
    def test_genLowerOrder(self): #seq, order):#Deprecated
        pass
        
    def test_higherOrderDict(self): #(seq, maxorder):
        self.assertRaises(ValueError, _higherOrderDict, [0,1,2], 3)
        self.assertRaises(TypeError, _higherOrderDict, 'ATGCF', 2)
        self.assertRaises(ValueError, _higherOrderDict, [0,4], 1)
        self.assertEqual(_higherOrderDict([0,1,2], 2), {0:[0,1,2], 1:[None, 4, 9], 2:[None, None, 36]})
    
    def test_SequenceDict(self):
        self.assertRaises(ValueError, SequenceDict, [0,1,2], 3)
        self.assertRaises(ValueError, SequenceDict, 'ATGCF', 2)
        self.assertRaises(ValueError, SequenceDict, [0,1,2,3,4], 5)
        self.assertRaises(ValueError, SequenceDict, [0,4], 1)
        self.assertEqual(SequenceDict('AACGTA', 5)['+'][0], [0,0,1,2,3,0])
        self.assertEqual(SequenceDict('AACGTA', 5)['-'][0], [3,0,1,2,3,3])
        self.assertEqual(str(SequenceDict('AACGTA', 5)), 'AACGTA')
        self.assertEqual(len(SequenceDict('AACGTA', 5)), 6)
        for strand in SequenceDict('AACGTA', 5):
        	self.assert_(strand in ['+', '-'])
        

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
	unittest.TextTestRunner(verbosity=2).run(suite)
