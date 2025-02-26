#!/bin/csh

bruk2pipe -in $1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 552.0 -dspfvs 20.0 -grpdly 68.0126800537109  \
  -xN              16384  \
  -xT              8192.0  \
  -xMODE            DQD  \
  -xSW          36231.884057971  \
  -xOBS         150.890990367  \
  -xCAR           99.41  \
  -xLAB             13C  \
  -ndim               1  \
| nmrPipe -fn MULT -c 0.0305176 \
| nmrPipe -fn EM -lb 1 -c 1.0 \
| nmrPipe -fn ZF -zf 3 \
| nmrPipe -fn FT \

sleep 1