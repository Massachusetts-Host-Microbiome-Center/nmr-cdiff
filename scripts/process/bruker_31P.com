#!/bin/csh

if ($# == 1) then
  set name = "$1_1D.fid"
else
  set name = "$1_$2_1D.fid"
endif

bruk2pipe -verb -in ./$1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 1376 -dspfvs 20 -grpdly 67.9841156005859  \
  -xN              1024  \
  -xT               512  \
  -xMODE            DQD  \
  -xSW        14534.884  \
  -xOBS         242.918  \
  -xCAR          -0.679  \
  -xLAB             31P  \
  -ndim               1  \
| nmrPipe -fn MULT -c 7.62939e-03 \
  -out ./$1/$name -ov

sleep 5