# TRACMT
Robust transfer function analysis code for magnetotellurics.
In TRACMT, you can use the RRMS estimator developed by Usui et al. (2024).

_Yoshiya Usui, Makoto Uyeshima, Shin'ya Sakanaka, Tasuku Hashimoto, Masahiro Ichiki, Toshiki Kaida, Yusuke Yamaya, Yasuo Ogawa, Masataka Masuda, Takahiro Akiyama, New robust remote reference estimator using robust multivariate linear regression, Geophysical Journal International, 2024, ggae199, [https://doi.org/10.1093/gji/ggae19_9](https://doi.org/10.1093/gji/ggae199)_

We developed the RRMS estimator by applying the robust multivariate linear regression S-estimator to the two-input-multiple-output relationship between the local EM field and the reference magnetic field that leads to the same equation as by the original remote reference method.
The RRMS estimator can give an unbiased estimate of MT transfer function and suppress the influence of outliers in the local electric field and magnetic
field.

# How to compile TRACMT
1) Download all source files of TRACMT to a directory.
2) Download source files of CLAPACK (https://www.netlib.org/clapack/) to another directory and make librariy files.
3) Copy library files of CLAPACK (blas_LINUX.a, lapack_LINUX.a, and libf2c.a) to "lib" directory and copy header files (blaswrap.h, clapack.h, and f2c.h) to "include" directory"
4) If your compiler supports C++11, you can compile TRACMT by "make -f Makefile_C++11" command.
   If your compiler does NOT support C++11, download mt19937-64.tgz from http://math.sci.hiroshima-u.ac.jp/m-mat/MT/mt64.html and rename mt19937-64.c to mt19937-64.cpp. After copying mt19937-64.cpp and mt64.h to the source-file directory of TRACMT, you can compile TRACMT by make command (Make -f Makefile).

# Functional overview
**Input file format**: Text (Ascii) file / .ats files of Metronix instruments / .dat files of ELOG-MT

**Prewhitening**: Non-robust prewhitening / Robust prewhitening / Robust prewhitening with robust filter

**Transfer function estimation method**: OLS / Non-robust remote reference / Robust remote reference / RRMS estimator

**Error estimation method**: Parametric approach / Bootstrap method / Jackknife method

# Release note
_**v1.2.0**_ Jun. 13, 2024: Initial release

