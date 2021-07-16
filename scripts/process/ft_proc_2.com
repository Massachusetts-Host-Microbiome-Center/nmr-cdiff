#!/bin/csh
#
set name = "$1_1D.fid"
set outname = "$1_1D.ft"

nmrPipe -in ./$name \
| nmrPipe -fn EM -lb 2.0 -c 1.0 \
| nmrPipe -fn ZF -zf 3 \
| nmrPipe -fn FT \
| nmrPipe -fn POLY -auto \
-out ./$outname -ov
