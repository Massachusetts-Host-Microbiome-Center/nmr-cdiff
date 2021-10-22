#!/bin/csh

if ($# == 1) then
  set name = "$1_1D.fid"
else
  set name = "$1_$2_1D.fid"
endif

bruk2pipe -verb -in ./$1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 1386.66666666667 -dspfvs 20 -grpdly 67.9878387451172  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW        14423.077  \
  -xOBS         600.086  \
  -xCAR           4.657  \
  -xLAB              1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 1.22070e-01 \
  -out ./$1/$name -ov

sleep 5
