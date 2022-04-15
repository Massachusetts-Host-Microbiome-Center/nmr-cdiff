#!/bin/csh
#
nmrPipe \
| nmrPipe -fn EM -lb 1.0 -c 1.0 \
| nmrPipe -fn ZF -zf 3 \
| nmrPipe -fn FT \
| nmrPipe -fn POLY -auto
