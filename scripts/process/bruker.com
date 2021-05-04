#!/bin/csh

set name = "$1_1D.fid"

bruk2pipe -verb -in ./$1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 552 -dspfvs 20 -grpdly 68.0126800537109  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW        36231.884  \
  -xOBS         150.906  \
  -xCAR         100.000  \
  -xLAB             13C  \
  -ndim               1  \
| nmrPipe -fn MULT -c 3.05176e-02 \
  -out ./$1/$name -ov

sleep 5
