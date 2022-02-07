#!/bin/csh

if ($# == 1) then
  set name = "$1_1D.fid"
else
  set name = "$1_$2_1D.fid"
endif

bruk2pipe -verb -in ./$1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 1664 -dspfvs 20 -grpdly 67.9841613769531  \
  -xN             16384  \
  -xT              8192  \
  -xMODE            DQD  \
  -xSW        12019.231  \
  -xOBS         600.086  \
  -xCAR           4.657  \
  -xLAB              1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 1.22070e-01 \
  -out ./$1/$name -ov

sleep 5
