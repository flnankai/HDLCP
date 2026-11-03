# HDLCP

R functions for **Analytically Calibrated Upper L-Statistics for High-Dimensional
Mean Change-Point Testing**, by Guowei Yan, Dan Zhuang and Long Feng.
The implementation follows manuscript **v37**.

## Install

```r
install.packages("remotes")
remotes::install_github("flnankai/HDLCP")
# Required only for DMS and ZWS22 comparisons:
remotes::install_github("ghwang-nk/DMS")
library(HDLCP)
```

Source installation requires a C++17 compiler (Rtools on Windows). R >= 4.1,
Rcpp and RcppArmadillo are required. CVXR is used only if a fitted covariance
kernel requires the paper's max-entry-norm correlation-matrix projection.

## Test your data

Rows of `X` are observations in time order; columns are variables. Supply a finite
numeric matrix. The test does not accept the true covariance, sparsity, innovation
distribution or change point as inputs.

```r
reference <- hdlcp_calibrate(
  n = nrow(X), p = ncol(X), k = c(1, 5), q = c(.1, .9),
  B = 2000, workers = 25, seed = 20261103
)
fit <- hdlcp_test(X, reference, keep_paths = TRUE)
fit$tests
hdlcp_location(fit)
plot(fit)
```

Reuse `reference` for all samples with these dimensions and scan settings. The
fixed-order critical curves use an independent Gaussian finite-design anchor
plus a Poisson tangent increment, interpolated in the estimated clustering
parameter and made nonincreasing in the tail probability. Proportional scans
reuse the standard normal vector bank; its transformation uses the estimated
full covariance kernel of each sample. No observed-data bootstrap or permutation
is used by HDLCP. Reference generation takes time; increasing `B` and decreasing
`floor_level` improves numerical precision at a computational cost.

`F1-0`, `F5-0.5`, `P0.1-0`, etc. identify component tests. `CL0` and `CL0.5`
combine the requested component p-values with equal Cauchy weights. Their location
is the leftmost scan maximizer of the component with the smallest **global**
p-value. Exact p-value ties follow the supplied fixed-order list, then the
proportion list. A location is the **last row before the change**.

The paper's Cauchy validity theorem covers a matched fixed/proportional pair;
the four-component combination used in its numerical studies is an empirical
extension. To use a matched pair, provide one `k` and one `q`.

## Recompute the paper

These functions generate new data and recompute the statistics. They do not read
stored Monte Carlo outcomes. Defaults reproduce the designs and replication
counts in Sections 6-7; Monte Carlo estimates need not equal printed values.

```r
paper_designs()                         # complete parameter/signal grids
size <- reproduce_size(workers = 25)   # Table 1: 8 null designs, 1000 each
fixed <- reproduce_order_power("fixed", workers = 25)       # Figure 1
prop <- reproduce_order_power("proportional", workers = 25) # Figure 2
comparison <- reproduce_comparison(workers = 25)            # Figures 3-4
plot_paper_results(fixed, "figure1.pdf")
plot_paper_results(prop, "figure2.pdf")
plot_paper_results(comparison, "figure3.pdf")
plot_paper_results(comparison, "figure4.pdf", metric = "location")

# All simulation tables and figures, with final summaries only:
results <- reproduce_paper("results", workers = 25)
```

Order comparisons use 10,000 independent null and 10,000 alternative draws.
Adaptive comparisons use 1,000 null and 1,000 alternative draws. At tau = .3,
signals are multiplied by 1.2. The q=1 benchmark is evaluated directly and
calibrated empirically under H0; no q=1 asymptotic p-value is used. Table 1 uses
analytic p-values, not size-corrected thresholds. LZZL20 uses B=500 multiplier
draws, s0=p/2, powers 1:5 and infinity, and a 0.2n boundary.

All parallel work uses independent seeded tasks and is reproducible across worker
counts. Serial execution is available with `workers = 1`. No simulation runs on
package load, and no replication dumps or checkpoint files are written. Only
explicitly requested final summaries and figures are saved. A small run can be
requested with `reproduce_paper(mode="smoke")`; this is for users to exercise the
pipeline, not a precision assessment.

## Individual methods

```r
hdlcp_fixed_test(X, reference, k = 5, gamma = 0)
hdlcp_proportional_test(X, reference, q = .1, gamma = .5)
cauchy_combine(c(.02, .3))
dms_test(X, gamma = 0)
dms_test(X, gamma = .5)
zws22_test(X)
lzzl20_test(X, B = 500, seed = 1)
```

The comparators are explicitly separate from the proposed method. DMS calls the
published package; ZWS22 uses an algebraically accelerated scan with that package's
reference distribution; LZZL20 uses a compiled multiplier-bootstrap calculation.
All use minimum-component-p localization as specified in the numerical study.

## SECOM application

Download `secom.data` and `secom_labels.data` from
[UCI SECOM](https://archive.ics.uci.edu/dataset/179/secom), then run:

```r
dat <- read_secom("secom.data", "secom_labels.data")
application <- reproduce_secom("secom.data", "secom_labels.data",
                               output_dir = "secom-results", workers = 25)
```

The loader orders timestamps, applies the manuscript's missingness and flat-run
filters, and selects rows 493:892. The application checks the expected 400 x 361
window, evaluates the reference split at 229, and recomputes Tables 4-5 and
Figures 5-6. `reps = 0` omits the 500-repetition side-wise subsampling experiment.
The original dataset is not redistributed.

## 中文使用说明

`hdlcp_calibrate()` 预制与观测数据无关的参考对象；`hdlcp_test(X, reference)`
返回各固定阶、比例阶与柯西联合检验结果；`hdlcp_location()` 返回变点位置。
统计量、全部参数估计、偏度稳定化和定位规则均按 v37 实现。
`reproduce_size()`、`reproduce_order_power()`、`reproduce_comparison()` 和
`reproduce_secom()` 分别重新计算论文中的 size、阶数选择、方法比较和实证结果。
默认并行 25 个工作进程。复现函数需要使用者主动调用，不会在安装或加载包时运行。

## Validation and license

Package checks cover compilation, installation, documentation and deterministic
formula/interface tests. The development checks do not rerun the paper's Monte
Carlo experiments. GPL (>= 3). See `inst/REFERENCES.md` and `inst/METHODS.md`
for method attribution and the equation-to-function map.
