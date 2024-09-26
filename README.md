# Systems Identification

The assignments focus on analyzing stationary stochastic processes and applying parameter estimation techniques using Auto-Regressive (AR) models. All MATLAB codes for simulations are included.

## Objectives

This assignment covers the following objectives:
- Analyze a stationary stochastic process using its expected value, variance, and covariance functions.
- Demonstrate the identification of Auto-Regressive (AR) models of increasing order.
- Compare theoretical findings with empirical estimates obtained through simulations using Least Squares (LS) methods.
- Investigate the effects of colored noise on model accuracy.

## Task A: Theoretical Analysis

### A1. Stationary Process Analysis

Consider a stationary stochastic process \( y(t) \) generated by:

\[
y(t) = \frac{1 + c z^{-1}}{1 - a z^{-1}} e(t)
\]

Where \( e(t) \) is white noise with zero mean and variance \( \lambda^2 \), and constants \( a \in \mathbb{R}, |a| < 1 \) and \( c \in \mathbb{R}, |c| < 1 \).

- **Goal**: Determine the following:
  1. The expected value \( E[y(t)] \) of the process.
  2. The variance \( \text{var}[y(t)] = E\{[y(t) - E[y(t)]]^2\} \).
  3. The covariance function \( \gamma_y(\tau) = E\{(y(t) - E[y(t)])(y(t - \tau) - E[y(t)])\} \) for \( \tau = 1, 2 \).

### A2. AR Model Identification

Using the same process as in A1, estimate the parameters using Auto-Regressive models:
- **AR(1) Model**: \( y(t) = a_1 y(t - 1) + \xi(t) \)
- **AR(2) Model**: \( y(t) = a_1 y(t - 1) + a_2 y(t - 2) + \xi(t) \)

Where \( \xi(t) \) is white noise.

- **Goal**: Use a Parameter Estimation Method (PEM) to find the optimal parameters \( \theta_1 = a_1 \) and \( \theta_2 = [a_1, a_2]^T \) for the AR(1) and AR(2) models, respectively. Calculate the prediction errors \( \text{var}[\varepsilon_{\theta_1}(t)] \) and \( \text{var}[\varepsilon_{\theta_2}(t)] \).

### A3. Non-Zero Mean Noise

Consider that the white noise \( e(t) \) is replaced by \( e(t) + 1 \) (non-zero mean). Perform the same analysis as in A1 and A2, and discuss how the results differ due to the shift in the mean of the noise.

### A4. Numerical Analysis

Given the specific values \( a = -1/2 \), \( c = 1/4 \), and \( \lambda^2 = 4 \):
- **Goal**: Compute the expected value, variance, and covariance from A1 and the optimal parameters from A2. Compare the prediction errors of the AR(1) and AR(2) models and discuss your findings.

## Task B: Simulations

### B1. White Noise Simulations

Using the recursive equation:

\[
y(t) = ay(t - 1) + e(t) + ce(t - 1)
\]

Where \( a = -1/2 \), \( c = 1/4 \), and \( \lambda^2 = 4 \), perform the following simulations:
1. **Generate a dataset** of \( N = 1000 \) samples.
2. **Apply the Least Squares Algorithm** to estimate the AR(1) and AR(2) parameters.
3. **Compare** the estimates with the theoretical values from A2 and comment on the differences.
4. **Repeat** for \( N = 2000 \) samples and discuss the results.

Finally, generate 200 independent datasets for \( N = 1000 \) and \( N = 2000 \), apply the Least Squares Algorithm, and compute the empirical mean and variance of the estimates.

### B2. Coloured Noise Simulations

Using coloured noise generated by the equation:

\[
e(t) = \xi(t) + \frac{1}{2}\xi(t - 1) + \frac{1}{4}\xi(t - 2)
\]

Repeat the simulation steps from B1. Discuss how the presence of colored noise affects the estimation accuracy of the AR(1) and AR(2) models.

---

This README provides a structured and clear overview of the assignments, without any references to your university, adhering to your request for rephrasing. You can add this file directly to your GitHub repository.
