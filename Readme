
# Gaussian Curve Fitting Using Newton–Gauss Method

**Puspa Kamal Rai** · Ist M.Sc. Physics (25010204011)  
Department of Physics, Prashanti Nilayam Campus  
Sri Sathya Sai Institute of Higher Education

---

## Overview

This project implements **nonlinear least-squares Gaussian curve fitting** using the **Newton–Gauss (Levenberg–Marquardt) method**, coded in **Scilab**. A synthetic noisy dataset over $x \in [-5, 5]$ is fitted to the three-parameter Gaussian model:

$$f(x;\,A,\mu,\sigma) = A\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$$

where $A$ is the peak amplitude, $\mu$ is the centroid, and $\sigma$ is the standard deviation.

---

## Algorithm

The Newton–Gauss update rule minimises the sum of squared residuals:

$$U(\mathbf{a}) = \sum_{i=1}^{N} \left[y_i - f(x_i;\mathbf{a})\right]^2$$

The iterative update at each step is:

$$\mathbf{a}_{k+1} = \mathbf{a}_k - \left(\mathbf{J}^\top\mathbf{J} + \lambda\mathbf{I}\right)^{-1}\mathbf{J}^\top\mathbf{r}$$

where $\mathbf{J}$ is the Jacobian of the residuals and $\lambda$ is the Levenberg–Marquardt damping parameter.

**Damping strategy:**
- If the trial step reduces the cost → accept step, decrease $\lambda$ (trust Gauss–Newton)
- If the trial step increases the cost → reject step, increase $\lambda$ (fall back to gradient descent)

**Convergence:** the iteration stops when $\|\Delta\mathbf{a}\| < 10^{-6}$.

---

## Repository Structure

```
GaussianFit/
├── gaussian_fit.sce        # Main Scilab script
├── Gaussiafit.png          # Output plot (fitted curve)
├── gaussian_fit_report.pdf # Full LaTeX report
├── gaussian_fit_report.tex # LaTeX source
└── README.md
```

---

## Results

Starting from the initial guess $\mathbf{a}_0 = (105,\,-0.4,\,1.2)^\top$, the algorithm converges to:

| Parameter | Description       | Initial Guess | Fitted Value   |
|-----------|-------------------|---------------|----------------|
| $A$       | Peak amplitude    | 105.0         | ≈ 102.1        |
| $\mu$     | Centroid (mean)   | −0.4          | ≈ −0.02        |
| $\sigma$  | Standard deviation| 1.2           | ≈ 1.20         |

![Gaussian Fit](Gaussiafit.png)

---

## How to Run

1. Open **Scilab**.
2. Load the script:
   ```scilab
   exec('gaussian_fit.sce', -1);
   ```
3. The fitted parameters `[A, mu, sigma]` are printed to the console and the fitted curve is plotted.

---

## Dependencies

- [Scilab](https://www.scilab.org/) (any recent version)
- No external toolboxes required

---

## Report

The full written report (theory, methodology, results) is available as [`gaussian_fit_report.pdf`](gaussian_fit_report.pdf), typeset in LaTeX.

---

## References

1. K. Levenberg, "A method for the solution of certain non-linear problems in least squares," *Quarterly of Applied Mathematics*, vol. 2, no. 2, pp. 164–168, 1944.
2. D. W. Marquardt, "An algorithm for least-squares estimation of nonlinear parameters," *Journal of the Society for Industrial and Applied Mathematics*, vol. 11, no. 2, pp. 431–441, 1963.

---

## Author

**[Puspa Kamal Rai](https://puspa-opal.vercel.app/)**  
Ist M.Sc. Physics · Sri Sathya Sai Institute of Higher Education
