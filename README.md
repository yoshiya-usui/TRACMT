# TRACMT
Robust transfer function analysis code for magnetotellurics.
In TRACMT, you can use the RRMS estimator developed by Usui et al. (2024).

_Yoshiya Usui, Makoto Uyeshima, Shin'ya Sakanaka, Tasuku Hashimoto, Masahiro Ichiki, Toshiki Kaida, Yusuke Yamaya, Yasuo Ogawa, Masataka Masuda, Takahiro Akiyama, New robust remote reference estimator using robust multivariate linear regression, Geophysical Journal International, 2024, ggae199, [https://doi.org/10.1093/gji/ggae19_9](https://doi.org/10.1093/gji/ggae199)_

We developed the RRMS estimator by applying the robust multivariate linear regression S-estimator to the two-input-multiple-output relationship between the local EM field and the reference magnetic field that leads to the same equation as by the original remote reference method.
The RRMS estimator can give an unbiased estimate of the MT transfer function and suppress the influence of outliers in the electric field and magnetic
field.

## How to compile TRACMT
1) Download all source files of TRACMT to a directory.
2) Download source files of CLAPACK (https://www.netlib.org/clapack/) to another directory and make library files.
3) Copy library files of CLAPACK (blas_LINUX.a, lapack_LINUX.a, and libf2c.a) to "lib" directory and copy header files (blaswrap.h, clapack.h, and f2c.h) to "include" directory.
4) If your compiler supports C++11, you can compile TRACMT by "make -f Makefile_C++11" command.
   If your compiler does NOT support C++11, download mt19937-64.tgz from http://math.sci.hiroshima-u.ac.jp/m-mat/MT/mt64.html and rename mt19937-64.c to mt19937-64.cpp. 
   After copying mt19937-64.cpp and mt64.h to the source-file directory of TRACMT, you can compile TRACMT by the make command (Make -f Makefile).
5) To read MTH5 files (Peacock et al., 2022), please use "Makefile_C++11_MTH5" or "Makefile_MTH5". Before compiling TRACMT, HDF5 library should be installed and path to HDF5 libraries should be added as follows.<br>
   export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/hdf5-1.10.5/lib"<br>
   export LIBRARY_PATH="${LIBRARY_PATH}:/usr/local/hdf5-1.10.5/lib"<br>
   export PATH="${PATH}:/usr/local/hdf5-1.10.5/lib"<br>

## Functional Overview
**Input file format**: Text (Ascii) file / .ats files of Metronix instruments / MTH5 files (Peacock et al., 2022) / .dat files of ELOG-MT

**Prewhitening**: Non-robust prewhitening / Robust prewhitening / Robust prewhitening with robust filter

**Transfer function estimation method**: Ordinary least square / Non-robust remote reference / Robust remote reference / RRMS estimator

**Error estimation method**: Parametric approach / Bootstrap method / Jackknife method

## Release note
_**v1.2.0**_ June. 13, 2024: Initial release.

_**v1.3.4**_ July. 14, 2024:  A parametric error estimation method was added for the RRMS estimator. Some bugs relating to the ELOG reading option and calibration were fixed.

_**v1.3.6**_ August. 15, 2024: A function to read ELOG-DUAL binary files was modified.

_**v2.0.0**_ January. 1, 2025: The fast and robust bootstrap method was added. Some equations of the parametric error estimation method were modified.

_**v2.1.0**_ January. 14, 2025: The equation of the AIC is modified.

_**v2.2.0**_ February. 10, 2025: The option for reading MTH5 files is added.


