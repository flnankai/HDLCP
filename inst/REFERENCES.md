# Method and implementation provenance

- Yan, Guowei; Zhuang, Dan; Feng, Long. Analytically Calibrated Upper L-Statistics
  for High-Dimensional Mean Change-Point Testing. Manuscript v37.
- Wang and Feng (2023), JRSS B, DOI: 10.1093/jrsssb/qkad048.
  DMS is an optional external GPL-3 dependency, https://github.com/ghwang-nk/DMS.
- Zhang, Wang and Shao (2022), JASA, DOI: 10.1080/01621459.2021.1884562.
  ZWS22's q=2,6 reference distributions come from DMS. The C++ scan included
  here was developed in this research project and implements the same statistics.
- Liu, Zhou, Zhang and Liu (2020), JRSS B, DOI: 10.1111/rssb.12375.
  The compiled multiplier-bootstrap norm calculation was developed in this
  research project; it is not represented as the original authors' software.
- UCI Machine Learning Repository, SECOM (2008), dataset 179.
  Dataset supplied separately by the user; not bundled in HDLCP.

Reusable numerical helpers and accelerated comparator kernels were adapted from
the research project's R/C++ files. Public inference and nuisance estimates are
specified by manuscript v37. The package uses GPL (>= 3).
