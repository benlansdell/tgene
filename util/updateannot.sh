chmod 755 ./gfftoucscgff.py
chmod 755 ./inttowig.py
#./inttowig.py ../data/tiling/adh/pm-mm/smc/fly.smc.pmmm.sum.combined.int
./inttowig.py ../data/tiling/adh/pm-mm/mismc/pmmm.mismc.fly.combined.int
#./inttowig.py ../data/tiling/adh/gcrma/averaged/chr2L.ints.added.fly.combined
#./inttowig.py ../data/tiling/adh/gcrma/averaged/gcrma.adh.int.expressed..fly.combined
#./inttowig.py ../data/tiling/adh/gcrma/averaged/gcrma.adh.int.mismc.expressed..fly.combined
#./inttowig2.py ../data/tiling/raw/gcrma/correlation/temp
#mv temp.wig adh2.mismc.int.wig
#./gfftoucscgff.py ../decode/drosophila1/flyA2.probelevel.fa.gff.combined
#./gfftoucscgff.py ../decode/drosophila3/flyTG10a.fa.gff
#./gfftoucscgff.py ../decode/drosophila3/flyTG1a.fa.gff
#./gfftoucscgff.py ../decode/drosophila3/fly42a.fa.gff
#./gfftoucscgff.py ../decode/drosophila1/flyA3.fa.gff
#./gfftoucscgff.py ../decode/drosophila1/flyA4.fa.gff
#./gfftoucscgff.py ../results/transfrags/manaketal.adh2.combined.gff
#./gfftoucscgff.py ../data/annotation/adh2.flybase.ucsc.2004.expressed.gff
#./gfftoucscgff.py ../decode/drosophila3/fly42.
#./combinefiles.py ./gcrma.adh.int.norand.expressed..fly.combined.wig ./gcrma.adh.int.mismc.expressed..fly.combined.wig ./fly.smc.pmmm.combined.int.wig ./fly.smc.pmmm.sum.combined.int.wig ./chr2L.ints.added.fly.combined.wig ./adh.ucsc.fa.scr.out ./adh.flybase.2004.better.expressed.gff ./fly80.fa.gff ./fly14a.fa.gff manaketal.combined.gff ./fly42.fa.gff
./combinefiles.py ./flyTG1a.fa.gff ./pmmm.mismc.fly.combined.int.wig 
scp *.ucsc unix22:/var/www/bioinf/folders/lansdell/dmel
