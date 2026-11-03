# v37 equation-to-function map

| Manuscript quantity | Implementation |
|---|---|
| All n-1 adjacent differences, Eq. (1) | `.prepare`: `(X[i+1,]-X[i,])/sqrt(2)`, `colMeans(D^2)` |
| Squared studentized CUSUM, upper ranks | `.prepare`, `hdlcp_test` |
| MOM correlations, Eq. (14) | `.summaries`: LR consecutive blocks separated by one unused difference row |
| Five-block trace estimate, Eq. (15) | `.summaries`: exclude equal and adjacent row pairs within every block |
| Gaussian MOM noise floor, Eq. (16) | `hdlcp_calibrate`: same deterministic blocks, independent Gaussian columns |
| Marginal tail/density, Eqs. (17)-(18) | `fixed_marginal_tail`, `fixed_marginal_density` |
| Rank-specific strong/weak pair clustering, Eq. (19) | `.theta`, Gaussian bivariate quadrature |
| Poisson tangent field, Eq. (9), Section 2.6 | `tangent_bank_cpp`: first-appearance roots, full Brownian trajectories, fixed geometric multiplicities; all requested k use the same paths |
| Monotone anchored critical curve, Eq. (20) | `.fixed_curve`: reverse cumulative maximum after interpolation |
| Fixed-order p-value/rejection, Eq. (21) | `.curve_p`, `hdlcp_test` |
| Full covariance and intensity | `.kernel`: Gamma plus trace term plus every selected full strong-pair correction |
| Gaussian scale center, Eq. (25) | `hdlcp_test`, `.order_scores` |
| Proportional statistics, Eqs. (26)-(27) | `hdlcp_test` |
| Ridge and clipped skewness, Eq. (28) | `hdlcp_test` |
| Gaussian skewness transport, Eq. (29) | `hdlcp_test`, reusable standard normal bank |
| Gamma=1/2 intensity/Gumbel calibration | `hdlcp_test` |
| Minimum-component-p localization, Eq. (32) | `hdlcp_test`, `hdlcp_location` |
| Endpoint/degenerate-scale conventions, Appendix A | `cauchy_combine`, `hdlcp_test` |

Numerical reference conventions: Gaussian reference tails evaluate `1-F_B(T)`
as the fraction of reference maxima strictly greater than T, as in Section 3.3.
There is no add-one adjustment. Fixed-order quantiles use type 8 on a tail grid
including 0 and 1, followed by the monotone envelope and generalized inverse.
The fixed-order rejection is evaluated directly as `L > c(alpha)` (Eq. 21),
including equality at a critical value. Values below/above the finite reference
range are 1/0, rather than being clipped to an interior probability.
Finite banks have finite tail resolution: numerical 0/1 values are flagged in
`diagnostics$reference_boundary` and must not be read as exact theoretical tails.
Analytic Gumbel tails retain log p-values for minimum-p selection without
floating-point underflow. Exact ties still use the supplied component order.

Tangent height truncation is a numerical integration tolerance; it does not
alter the observed upper-L statistic. At each grid time r, Poisson centers above
the cutoff are sampled and their backward Brownian paths are accepted only if
all earlier grid heights lie below that cutoff. These disjoint first-appearance
sets partition the restriction of the paper's Poisson path measure. An accepted
center is then followed at every later grid time, including below-cutoff
excursions and returns. Its geometric multiplicity remains attached to that
center. No independent replacement particles are introduced on return. The
sampler stops with an error if the cutoff does not leave k visible particles
at every grid time; it never silently conditions the reference on success. The
stored clustering grid is exactly the grid in Section 2.6. Bivariate pair-tail
interpolation uses rho=0,0.005,...,0.95 with the stated endpoint clipping.
Hermite covariance functions are integrated by Gaussian conditional moments and
Gauss-Legendre quadrature. The derivative is integrated directly through
`Gamma'(rho)=E[h'(Z1)h'(Z2)]`, where `h(z)=psi_q(z^2)`; it is odd in rho.
Gamma'(1)=4 mu is imposed analytically.

The raw estimated correlation kernel is retained when positive semidefinite.
Nontrivial indefiniteness is handled by the manuscript's max-entry-norm
semidefinite optimization using CVXR. Tiny negative eigenvalues within 1e-7
are treated as floating-point error when forming its square root.

The simulator knows Sigma only to generate data. Inference receives X and a
universal reference object; no AR coefficient or other true nuisance parameter
enters `hdlcp_test`.
