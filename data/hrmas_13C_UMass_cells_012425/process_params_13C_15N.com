#!/bin/csh

bruk2pipe -in $1/fid \
  -bad 0.0 -ext -aswap -AMX -decim 512.0 -dspfvs 20.0 -grpdly 67.983642578125  \
  -xN              16384  \
  -xT              8192.0  \
  -xMODE            DQD  \
  -xSW          39062.5  \
  -xOBS         150.902809  \
  -xCAR           99.41  \
  -xLAB             13C_15N  \
  -ndim               1  \
| nmrPipe -fn MULT -c 0.0305176 \
| nmrPipe -fn EM -lb 0.01 -c 1.0 \
| nmrPipe -fn ZF -zf 3 \
| nmrPipe -fn FT \

sleep 1