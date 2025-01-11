# Multivariate Generalized Gaussian Distribution (MVGGD)
The MVGGD repository provides MATLAB code for calculating the probability density of a *d*-dimensional, Generalized Gaussian Distribution at realization vector $\boldsymbol{x}$ with mean vector $\boldsymbol{\mu}$, covariance matrix $\Sigma$, and shaping parameter $\beta$. The probability density function defined in [1] is:  <br>

$$
\begin{gather}
    p(\mathbf{x}|\boldsymbol{\mu},\Sigma,\beta) =  A(\beta,d)
    \exp(-[B(\beta,d)(\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1}(\mathbf{x}-\boldsymbol{\mu})]^{\beta}),
    \\ 
    \text{where}\quad A(\beta,d)= \Big(\frac{B(\beta,d)}{\pi}\Big)^{\frac{d}{2}}\cdot\frac{\Gamma(\frac{d}{2})\beta}{\Gamma(\frac{d}{2\beta})|\Sigma|^{\frac{1}{2}}}\quad\text{and}\quad B(\beta,d)= \frac{\Gamma(\frac{d+2}{2\beta})}{d\Gamma(\frac{d}{2\beta})}, 
\end{gather}
$$

When $\beta = 1$, *mvggd.m* becomes *mvnpdf.m*. Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1](https://www.tandfonline.com/doi/full/10.1080/03610926.2024.2439999#abstract) L. Hanson, B., Zhao, M., & Thomas, R. B. (2025). An extensible framework for the probabilistic search of stochastically-moving targets characterized by generalized Gaussian distributions or  experimentally-defined regions of interest. Communications in Statistics - Theory and Methods, 1–26. https://doi.org/10.1080/03610926.2024.2439999
