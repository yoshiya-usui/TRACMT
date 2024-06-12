# TRACMT
Robust transfer function analysis code for magnetotellurics
In TRACMT, you can use the RRMS estimator developed by Usui et al. (2024).

_Yoshiya Usui, Makoto Uyeshima, Shin'ya Sakanaka, Tasuku Hashimoto, Masahiro Ichiki, Toshiki Kaida, Yusuke Yamaya, Yasuo Ogawa, Masataka Masuda, Takahiro Akiyama, New robust remote reference estimator using robust multivariate linear regression, Geophysical Journal International, 2024, ggae199, [https://doi.org/10.1093/gji/ggae19_9](https://doi.org/10.1093/gji/ggae199)_

We developed the RRMS estimator by applying the robust multivariate linear regression S-estimator to the two-input-multiple-output relationship between the local EM field and the reference magnetic field that leads to the same equation as by the original remote reference method.
The RRMS estimator can give an unbiased estimate of MT transfer function and suppress the influence of outliers in the local electric field and magnetic
field.
