#!/bin/csh

bruk2pipe -verb -in $1/fid \
  -bad 0.0 -ext -noaswap -AMX -decim 24 -dspfvs 12 -grpdly 0  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW         7002.801  \
  -xOBS         499.842  \
  -xCAR           4.754  \
  -xLAB              1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 6.25000e+01 \
| nmrPipe -fn ZF -zf 2 \
| nmrPipe -fn FT \
