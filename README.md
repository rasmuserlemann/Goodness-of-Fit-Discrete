# Power study code in the paper *Conditional Goodness-of-Fit Tests for Discrete Distributions*
The code is written in R. It was written in collaboration with Prof. Bo Henry Lindqvist.

## Anderson-Darling, Cramer-von Mises and Kolmogorov Smirnov test statistics
Test statistics functions are AD, CM, KS. The first two are regular quadratic goodness-of-fit test statistics and the last one is a maximal type test. Only null hypothesis model maximum likelihood estimates are used. For the geometric distribution, this is calculated by n/(sum(d)+n).

## Likelihood based tests
Test statistics functions are CR, SB, SBabs, theta, SW, SWL, SWU. The first one is the so-called full likelihood test. Its main advantage is in detecting deviations from homogeneity in the sample. The next 3 tests use the beta-geometric distribution as the alternative in the likelihood ratio test. It detects deviations towards the beta-geometric distribution. Last 3 tests are also likelihood ratio tests with the Discrete Weibull type I as the alternative distribution.