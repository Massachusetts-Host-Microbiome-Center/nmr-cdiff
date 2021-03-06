#!/bin/csh
#
set name = "$1_1D.ft.ps"
set outname = "$1_1D.ft.ps.bl"

nmrPipe -in ./$name \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
-out ./$outname -ov
