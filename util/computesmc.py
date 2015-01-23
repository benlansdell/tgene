#!/usr/bin/env python

from common.tiling import *

def computeStats():
        print 'Loading data...'
        #a = ExpressionData(['../data/tiling/raw/gcrma/average/chr2L.ints1out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints6out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints9out.adh0',
        #                    '../data/tiling/raw/gcrma/average/chr2L.ints12out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints15out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints18out.adh0',
        #                    '../data/tiling/raw/gcrma/average/chr2L.ints21out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints24out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints27out.adh0',
        #                    '../data/tiling/raw/gcrma/average/chr2L.ints30out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints33out.adh0','../data/tiling/raw/gcrma/average/chr2L.ints36out.adh0'], '../results/transfrags/manaketal.combined.gff')
        a = ExpressionData(['../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_10_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_11_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_12_C01.sig.gr',
                            '../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_1_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_2_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_3_C01.sig.gr',
                            '../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_4_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_5_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_6_C01.sig.gr',
                            '../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_7_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_8_C01.sig.gr','../data/tiling/raw/pm-mm/sliding-window/chr2L/Dro_Total_AS_9_C01.sig.gr'], '../data/annotation/adh.flybase.2004.combined.gff.expressed')
        #a = ExpressionData(['../data/tiling/raw/pm-mm/Dro_Total_AS_10_T1.sig.gr0','../data/tiling/raw/pm-mm/Dro_Total_AS_11_T1.sig.gr3','../data/tiling/raw/pm-mm/Dro_Total_AS_12_T4.sig.gr6',
        #                    '../data/tiling/raw/pm-mm/Dro_Total_AS_1_T1.sig.gr9','../data/tiling/raw/pm-mm/Dro_Total_AS_2_T1.sig.gr12','../data/tiling/raw/pm-mm/Dro_Total_AS_3_T1.sig.gr15',
        #                    '../data/tiling/raw/pm-mm/Dro_Total_AS_4_T1.sig.gr18','../data/tiling/raw/pm-mm/Dro_Total_AS_5_T1.sig.gr21','../data/tiling/raw/pm-mm/Dro_Total_AS_6_T1.sig.gr24',
        #                    '../data/tiling/raw/pm-mm/Dro_Total_AS_7_T1.sig.gr27','../data/tiling/raw/pm-mm/Dro_Total_AS_8_T1.sig.gr30','../data/tiling/raw/pm-mm/Dro_Total_AS_9_T1.sig.gr33'], '../results/transfrags/manaketal.adh2.combined.gff')
        #a = ExpressionData(['../data/tiling/raw/gcrma/chr2L.ints1out.adh','../data/tiling/raw/gcrma/chr2L.ints2out.adh','../data/tiling/raw/gcrma/chr2L.ints3out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints4out.adh','../data/tiling/raw/gcrma/chr2L.ints5out.adh','../data/tiling/raw/gcrma/chr2L.ints6out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints7out.adh','../data/tiling/raw/gcrma/chr2L.ints8out.adh','../data/tiling/raw/gcrma/chr2L.ints9out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints10out.adh','../data/tiling/raw/gcrma/chr2L.ints11out.adh','../data/tiling/raw/gcrma/chr2L.ints12out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints13out.adh','../data/tiling/raw/gcrma/chr2L.ints14out.adh','../data/tiling/raw/gcrma/chr2L.ints15out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints16out.adh','../data/tiling/raw/gcrma/chr2L.ints17out.adh','../data/tiling/raw/gcrma/chr2L.ints18out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints19out.adh','../data/tiling/raw/gcrma/chr2L.ints20out.adh','../data/tiling/raw/gcrma/chr2L.ints21out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints22out.adh','../data/tiling/raw/gcrma/chr2L.ints23out.adh','../data/tiling/raw/gcrma/chr2L.ints24out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints25out.adh','../data/tiling/raw/gcrma/chr2L.ints26out.adh','../data/tiling/raw/gcrma/chr2L.ints27out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints28out.adh','../data/tiling/raw/gcrma/chr2L.ints29out.adh','../data/tiling/raw/gcrma/chr2L.ints30out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints31out.adh','../data/tiling/raw/gcrma/chr2L.ints32out.adh','../data/tiling/raw/gcrma/chr2L.ints33out.adh',
        #                    '../data/tiling/raw/gcrma/chr2L.ints34out.adh','../data/tiling/raw/gcrma/chr2L.ints35out.adh','../data/tiling/raw/gcrma/chr2L.ints36out.adh'], '../data/annotation/adh.flybase.2004.combined.gff.expressed')
        print 'Analysing data...'
        analyseData(a,'../data/tiling/raw/gcrma/correlation/mismc.sliding.adh.int.allexprssedgenes')

computeStats()