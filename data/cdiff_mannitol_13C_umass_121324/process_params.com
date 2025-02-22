#!/bin/csh

bruk2pipe -in $1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 512 -dspfvs 20 -grpdly 67.9937744140625  \
  -xN              32768  \
  -xT              16384.0  \
  -xMODE            DQD  \
  -xSW          12019.2307692308  \
  -xOBS         600.13  \
  -xCAR           99.41  \
  -xLAB             1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 0.0305176 \
| nmrPipe -fn EM -lb 1 -c 1.0 \
| nmrPipe -fn ZF -zf 3 \
| nmrPipe -fn FT \

sleep 1