#!/usr/bin/env sh

rm -rf bin/;

mkdir bin/;

cd include
  echo 'building zlib'
  cd zlib-*/;
  make clean;
  ./configure --prefix=`pwd`
  make;
  cd ..;

  echo 'building samtools'
  cd samtools-*/;
    make clean;
    make;
  cd ..;
cd ..;

cd src/;
  echo 'building ATCGmapMerge'
  g++ ATCGmapMerge.cpp -o ../bin/ATCGmapMerge -lz -L ../include/zlib-1.2.8 

  echo 'building CGmapSelectByRegion'
  g++ CGmapSelectByRegion.cpp -o ../bin/CGmapSelectByRegion

  echo 'building CGmapMethInBed'
  g++ CGmapMethInBed.cpp  -o ../bin/CGmapMethInBed  

  echo 'building CGmapMethInFragReg'
  g++ CGmapMethInFragReg.cpp -o ../bin/CGmapMethInFragReg
  
  echo 'building CGmapFromBAM'
  gcc -o ../bin/CGmapFromBAM CGmapFromBAM.c -lz -L ../include/zlib-1.2.8 -lbam -L ../include/samtools-0.1.18;
  
  echo 'building CGmapToCGbz'
  gcc -o ../bin/CGmapToCGbz CGmapToCGbz.c  -L ../include/zlib-1.2.8  -lbam -L ../include/samtools-0.1.18 -lz
  
  echo 'building CGbzToCGmap'
  gcc -o ../bin/CGbzToCGmap CGbzToCGmap.c  -L ../include/zlib-1.2.8  -lbam -L ../include/samtools-0.1.18 -lz
  
  echo 'building ATCGmapToATCGbz'
  gcc -o ../bin/ATCGmapToATCGbz ATCGmapToATCGbz.c  -L ../include/zlib-1.2.8  -lbam -L ../include/samtools-0.1.18 -lz
  
  echo 'building ATCGbzToATCGmap'
  gcc -o ../bin/ATCGbzToATCGmap ATCGbzToATCGmap.c  -L ../include/zlib-1.2.8  -lbam -L ../include/samtools-0.1.18 -lz
  
  echo 'building CGbzFetchRegion'
  gcc -o ../bin/CGbzFetchRegion CGbzFetchRegion.c -L ../include/zlib-1.2.8 -lbam -L ../include/samtools-0.1.18 -lz
  
  echo 'building ATCGbzFetchRegion'
  gcc -o ../bin/ATCGbzFetchRegion ATCGbzFetchRegion.c -L ../include/zlib-1.2.8 -lbam -L ../include/samtools-0.1.18 -lz
cd ..

cd bin
  for pr in ../src/*.py ../src/*.pl ../src/*.R; do
    pn=`basename $pr | cut -d"." -f1`
    echo 'building '`basename $pn`
    ln -s $pr $pn
    chmod +x $pn
  done
cd ..

