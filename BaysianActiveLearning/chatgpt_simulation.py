#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 17:52:28 2026

@author: vencov
"""

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Psychometric function (logistic)
# -----------------------------
def psychometric(x, alpha, beta=0.5):
    return 1 / (1 + np.exp(-beta * (x - alpha)))


# -----------------------------
# Simulated true listener
# -----------------------------
alpha_true = 20      # true threshold
beta_true = 0.5      # true slope

def simulate_response(x):
    p = psychometric(x, alpha_true, beta_true)
    return np.random.rand() < p


# -----------------------------
# Bayesian setup
# -----------------------------
alpha_grid = np.linspace(0, 40, 400)
posterior = np.ones_like(alpha_grid)
posterior /= np.sum(posterior)

stim_levels = np.linspace(0, 40, 50)
beta_model = 0.5

n_trials = 25

print("Starting Bayesian active learning...\n")

for trial in range(n_trials):

    entropies = []

    # Evaluate each possible stimulus level
    for x in stim_levels:
        p_detect = psychometric(x, alpha_grid, beta_model)

        p_yes = np.sum(p_detect * posterior)
        p_no = 1 - p_yes

        # Posterior if YES
        post_yes = posterior * p_detect
        post_yes /= np.sum(post_yes)

        # Posterior if NO
        post_no = posterior * (1 - p_detect)
        post_no /= np.sum(post_no)

        H_yes = -np.sum(post_yes * np.log(post_yes + 1e-12))
        H_no = -np.sum(post_no * np.log(post_no + 1e-12))

        H_exp = p_yes * H_yes + p_no * H_no
        entropies.append(H_exp)

    # Choose stimulus minimizing expected entropy
    x_next = stim_levels[np.argmin(entropies)]

    # Simulate listener response
    response = simulate_response(x_next)

    # Update posterior
    likelihood = psychometric(x_next, alpha_grid, beta_model)

    if response:
        posterior *= likelihood
        resp_text = "YES"
    else:
        posterior *= (1 - likelihood)
        resp_text = "NO"

    posterior /= np.sum(posterior)

    estimate = alpha_grid[np.argmax(posterior)]

    print(f"Trial {trial+1:02d}: Stimulus={x_next:.2f} dB | Response={resp_text} | Estimate={estimate:.2f}")

print("\nFinished.")
print(f"True threshold: {alpha_true}")
print(f"Final estimate: {estimate:.2f}")


# -----------------------------
# Plot final posterior
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(alpha_grid, posterior, label="Posterior")
plt.axvline(alpha_true, linestyle='--', label="True threshold")
plt.xlabel("Threshold (alpha)")
plt.ylabel("Probability")
plt.title("Posterior Distribution After Active Learning")
plt.legend()
plt.grid(True)
plt.show()