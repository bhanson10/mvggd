# Multivariate Generalized Gaussian Distribution (MVGGD)
The MVGGD repository provides MATLAB code for calculating the probability density of a *d*-dimensional, Generalized Gaussian Distribution at realization vector $\boldsymbol{x}$ with mean vector $\boldsymbol{\mu}$, covariance matrix $\Sigma$, and shaping parameter $\beta$. The full formulation can be found  <br><br>

$$
\begin{gather}
    p(\mathbf{x}|\boldsymbol{\mu},\Sigma,\beta) =  A(\beta,d)
    \exp\{-[B(\beta,d)(\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1}(\mathbf{x}-\boldsymbol{\mu})]^{\beta}\},
    \\ 
    \text{where}\quad A(\beta,d)= \Big(\frac{B(\beta,d)}{\pi}\Big)^{\frac{d}{2}}\cdot\frac{\Gamma(\frac{d}{2})\beta}{\Gamma(\frac{d}{2\beta})|\Sigma|^{\frac{1}{2}}}\quad\text{and}\quad B(\beta,d)= \frac{\Gamma(\frac{d+2}{2\beta})}{d\Gamma(\frac{d}{2\beta})}, 
\end{gather}
$$
