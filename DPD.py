import numpy as np
import matplotlib.pyplot as plt
import json

# --- Configuration Parameters ---
M       = 2         # Memory depth (number of past samples to include)
K       = 3         # Polynomial order (max power of amplitude term)
N_train = 5000      # Number of samples to generate for training
N_test  = 2000      # Number of samples for inference/testing
np.random.seed(0)   # Fix seed for reproducibility

# --- PA model coefficients (complex to simulate AM/AM + AM/PM) ---
# b[(m, k)] multiplies x[n-m] * |x[n-m]|^(k-1)
# --- PA model coefficients (complex to simulate AM/AM + AM/PM) ---
# b[(m, k)] multiplies x[n-m] * |x[n-m]|^(k-1)
b = {
    # m = 0 (instantaneous) terms
    (0, 1): 1.00   + 0.0j,    # linear gain
    (0, 2): 0.20   + 0.05j,   # 2nd‑order AM/AM + small AM/PM
    (0, 3): 0.10   - 0.02j,   # 3rd‑order

    # m = 1 (one‑sample delay) terms
    (1, 1): 0.80   + 0.10j,   # memory‑linear with phase twist
    (1, 2): 0.05   + 0.02j,   # memory 2nd‑order
    (1, 3): 0.02   + 0.01j,   # memory 3rd‑order

    # m = 2 (two‑sample delay) terms – now non‑zero to capture longer memory
    (2, 1): 0.60   + 0.08j,   # two‑tap linear
    (2, 2): 0.03   + 0.015j,  # two‑tap 2nd‑order
    (2, 3): 0.01   + 0.005j   # two‑tap 3rd‑order
}

def pa_model(x):
    """
    Simulate a memory‑polynomial Power Amplifier (PA).
    Inputs:
      x : complex numpy array of baseband samples
    Returns:
      y : complex numpy array, PA output same length as x
    Formula:
      y[n] = sum_{m=0..M} sum_{k=1..K} b[m,k] * x[n-m] * |x[n-m]|^(k-1)
    """
    # Pre‑allocate output array (complex)
    y = np.zeros_like(x, dtype=complex)
    # Loop over each output sample
    for n in range(len(x)):
        # Sum contributions from current + past M inputs
        for m in range(M+1):
            if n < m:
                continue  # skip if index negative
            xm = x[n-m]  # delayed input sample
            # Add each polynomial order term
            for k in range(1, K+1):
                y[n] += b[(m, k)] * xm * (np.abs(xm) ** (k-1))
    return y

# === TRAINING PHASE ===

# 1) Create random complex training waveform (e.g. QAM-like)
x_train = (np.random.randn(N_train) + 1j*np.random.randn(N_train)) / np.sqrt(2)

# 2) Pass through PA to get distorted training output
d_train = pa_model(x_train)

# 3) Build the predistorter basis matrix from PA output d_train
#    Each column = d_train delayed by m, multiplied by its own magnitude^(k-1)
cols_pred = []
for m in range(M+1):
    # Prepare delayed vector, zero-padding first m samples
    col = np.zeros_like(d_train)
    if m > 0:
        col[m:] = d_train[:-m]
    else:
        col[:] = d_train
    # For each polynomial order k, form basis term
    for k in range(1, K+1):
        cols_pred.append(col * (np.abs(col) ** (k-1)))

# Stack basis columns into a big matrix (N_train × ((M+1)*K))
X_pred = np.column_stack(cols_pred)

# 4) Solve least-squares: find DPD coefficients a_hat such that
#       X_pred @ a_hat ≈ x_train
a_hat, *_ = np.linalg.lstsq(X_pred, x_train, rcond=None)

# 5) Save the complex DPD coefficients
np.save('dpd_coeffs.npy', a_hat)  # numeric reload
# Also export to JSON as [real, imag] pairs for each (m,k)
coeff_dict = {}
for idx, (m,k) in enumerate([(m,k) for m in range(M+1) for k in range(1,K+1)]):
    c = a_hat[idx]
    coeff_dict[f"{m},{k}"] = [c.real, c.imag]
with open('dpd_coeffs.json', 'w') as f:
    json.dump(coeff_dict, f, indent=2)
print("Saved DPD coefficients to dpd_coeffs.npy and dpd_coeffs.json")

# === INFERENCE PHASE ===

# Reload DPD coefficients (complex)
coeffs = np.load('dpd_coeffs.npy')

# 1) Generate new random complex test waveform
x_test = (np.random.randn(N_test) + 1j*np.random.randn(N_test)) / np.sqrt(2)

# 2) PA output without any predistortion
y_no_dpd = pa_model(x_test)

# 3) Apply static DPD to x_test
def dpd_apply(x, coefs):
    """
    Build the same memory‑polynomial basis from x, then multiply by coefs.
    Returns predistorted x.
    """
    cols = []
    for m in range(M+1):
        col = np.zeros_like(x, dtype=complex)
        if m > 0:
            col[m:] = x[:-m]
        else:
            col[:] = x
        for k in range(1, K+1):
            cols.append(col * (np.abs(col) ** (k-1)))
    X_inf = np.column_stack(cols)
    return X_inf @ coefs

# Generate predistorted waveform and re‑amplify
x_predistorted = dpd_apply(x_test, coeffs)
y_dpd = pa_model(x_predistorted)

# === PERFORMANCE METRICS ===

# Complex‑baseband MSE (amplitude + phase)
mse_before = np.mean(np.abs(y_no_dpd - x_test)**2)
mse_after  = np.mean(np.abs(y_dpd    - x_test)**2)
print(f"Baseband MSE before DPD: {mse_before:.6f}")
print(f"Baseband MSE after  DPD: {mse_after:.6f}")

# Separate AM/AM and AM/PM MSE components
mse_am_before = np.mean((np.abs(y_no_dpd) - np.abs(x_test))**2)
mse_am_after  = np.mean((np.abs(y_dpd)    - np.abs(x_test))**2)
mse_ph_before = np.mean((np.angle(y_no_dpd) - np.angle(x_test))**2)
mse_ph_after  = np.mean((np.angle(y_dpd)    - np.angle(x_test))**2)
print(f"AM/AM MSE before: {mse_am_before:.6f}, after: {mse_am_after:.6f}")
print(f"AM/PM MSE before: {mse_ph_before:.6f}, after: {mse_ph_after:.6f}")

# ============================================== PLOTTING ===================================================

# --- SMOOTH PLOTS FOR AM/AM, AM/PM, AND POWER TRANSFER ---

import numpy as np
import matplotlib.pyplot as plt

# Compute instantaneous metrics
mag_in       = np.abs(x_test)
mag_out_no   = np.abs(y_no_dpd)
mag_out_dp   = np.abs(y_dpd)
phase_err_no = np.angle(y_no_dpd) - np.angle(x_test)
phase_err_dp = np.angle(y_dpd)    - np.angle(x_test)

# Convert to dB for power plots
P_in_dB     = 20*np.log10(mag_in)
P_out_no_dB = 20*np.log10(mag_out_no)
P_out_dp_dB = 20*np.log10(mag_out_dp)

# Define bins
n_bins      = 30
bins        = np.linspace(P_in_dB.min(), P_in_dB.max(), n_bins+1)
inds        = np.digitize(P_in_dB, bins)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Preallocate arrays for medians
am_am_no_med  = np.full(n_bins, np.nan)
am_am_dp_med  = np.full(n_bins, np.nan)
am_pm_no_med  = np.full(n_bins, np.nan)
am_pm_dp_med  = np.full(n_bins, np.nan)
pw_no_med     = np.full(n_bins, np.nan)
pw_dp_med     = np.full(n_bins, np.nan)

# Populate medians per bin
for i in range(1, n_bins+1):
    mask = (inds == i)
    if mask.sum() < 5:
        continue
    # AM/AM (voltage)
    am_am_no_med[i-1] = np.median(mag_out_no[mask])
    am_am_dp_med[i-1] = np.median(mag_out_dp[mask])
    # AM/PM (phase error)
    am_pm_no_med[i-1] = np.median(phase_err_no[mask])
    am_pm_dp_med[i-1] = np.median(phase_err_dp[mask])
    # Power transfer (dB)
    pw_no_med[i-1]    = np.median(P_out_no_dB[mask])
    pw_dp_med[i-1]    = np.median(P_out_dp_dB[mask])

# Plot
plt.figure(figsize=(15,4))

# 1) AM/AM: |v_out| vs P_in_dB
plt.subplot(1,3,1)
plt.plot(bin_centers, am_am_no_med, '-o', label='Before DPD')
plt.plot(bin_centers, am_am_dp_med, '-o', label='After  DPD')
plt.title('AM/AM vs Input Power')
plt.xlabel('Input power [dB]')
plt.ylabel('Median |v_out|')
plt.legend()

# 2) AM/PM: phase error vs P_in_dB
plt.subplot(1,3,2)
plt.plot(bin_centers, am_pm_no_med, '-o', label='Before DPD')
plt.plot(bin_centers, am_pm_dp_med, '-o', label='After  DPD')
plt.title('AM/PM vs Input Power')
plt.xlabel('Input power [dB]')
plt.ylabel('Median phase error [rad]')
plt.legend()

# 3) Power transfer: P_out_dB vs P_in_dB
plt.subplot(1,3,3)
plt.plot(bin_centers, pw_no_med, '-o', label='Before DPD')
plt.plot(bin_centers, pw_dp_med, '-o', label='After  DPD')
plt.title('Power Transfer')
plt.xlabel('Input power [dB]')
plt.ylabel('Median output power [dB]')
plt.legend()

plt.tight_layout()
plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import os
# import csv
# import json
# # --- Configuration Parameters ---
# M = 1        # Memory depth: number of past samples to include (0..M)
# K = 3        # Polynomial order: highest power of input amplitude
# N_train = 5000  # Number of samples to generate for training
# N_test = 2000   # Number of samples for inference/testing
# np.random.seed(0)  # Seed for reproducibility

# # --- Define the "PA" (Power Amplifier) model coefficients ---
# # Memory-polynomial coefficients: b[(m, k)] corresponds to the term for x[n-m] * |x[n-m]|^(k-1)
# b = {
#     (0, 1): 1.0,   # Linear gain term (no memory)
#     (0, 2): 0.2,   # AM/AM distortion term
#     (0, 3): 0.1,   # Higher-order distortion
#     (1, 1): 0.8,   # Memory term linear
#     (1, 2): 0.05,  # Memory term second order
#     (1, 3): 0.02   # Memory term third order
# }

# def pa_model(x):
#     """
#     Simulate the nonlinear PA behavior using a memory-polynomial model.
#     x: input waveform (1D numpy array)
#     returns: distorted output (same shape as x)
#     """
#     y = np.zeros_like(x)
#     # Iterate over each sample
#     for n in range(len(x)):
#         # Sum contributions from current and past M samples
#         for m in range(M + 1):
#             for k in range(1, K + 1):
#                 if n - m >= 0:
#                     # b[(m,k)] * x[n-m] * |x[n-m]|^(k-1)
#                     y[n] += b[(m, k)] * x[n - m] * (abs(x[n - m]) ** (k - 1))
#     return y

# # --- Training Phase: Generate data and fit DPD coefficients ---

# # 1) Generate and run PA
# x_train = np.random.randn(N_train)
# d_train = pa_model(x_train)

# # 2) Build basis from d_train (not x_train!)
# cols_pred = []
# for m in range(M+1):
#     col = d_train.copy()
#     col[:m] = 0
#     for k in range(1, K+1):
#         cols_pred.append(col * (np.abs(col) ** (k-1)))
# X_pred = np.column_stack(cols_pred)

# # 3) Solve X_pred · a_hat = x_train  (we want PA(DPD(x)) ≈ x)
# a_hat, *_ = np.linalg.lstsq(X_pred, x_train, rcond=None)

# # --- (Optionally) JSON dump for easy parsing in other languages ---
# json_file = 'dpd_coeffs.json'
# coeff_dict = {}
# idx = 0
# for m in range(M+1):
#     for k in range(1, K+1):
#         coeff_dict[f"{m},{k}"] = float(a_hat[idx])
#         idx += 1

# with open(json_file, 'w') as f:
#     json.dump(coeff_dict, f, indent=2)
# print(f"Saved JSON‑format coefficients to {json_file}")

# # 5) Save the fitted coefficients to disk for later use
# coeff_file = 'dpd_coeffs.npy'
# np.save(coeff_file, a_hat)
# print(f"DPD coefficients saved to: {coeff_file}")

# # --- Inference Phase: Apply DPD and measure performance ---

# # Load the saved coefficients
# coeffs = np.load(coeff_file)

# # Generate a fresh test waveform
# x_test = np.random.randn(N_test)

# # 1) PA output without any predistortion
# y_no_dpd = pa_model(x_test)

# # 2) Apply the predistorter: build inference matrix X_inf and multiply by coeffs
# def dpd_apply(x, coeffs):
#     """
#     Apply static DPD to input x using precomputed coefficients.
#     Returns the predistorted waveform.
#     """
#     cols_inf = []
#     for m in range(M + 1):
#         col = x.copy()
#         col[:m] = 0
#         for k in range(1, K + 1):
#             cols_inf.append(col * (np.abs(col) ** (k - 1)))
#     X_inf = np.column_stack(cols_inf)
#     return X_inf @ coeffs

# x_predistorted = dpd_apply(x_test, coeffs)

# # 3) Pass predistorted waveform through PA
# y_dpd = pa_model(x_predistorted)

# # --- Error Measurement ---
# # Compute error signals: difference between PA output and original input
# err_no_dpd = y_no_dpd - x_test
# err_dpd = y_dpd - x_test

# # Compute Mean Squared Error for before/after
# mse_no_dpd = np.mean(err_no_dpd**2)
# mse_dpd = np.mean(err_dpd**2)

# print(f"MSE before DPD: {mse_no_dpd:.6f}")
# print(f"MSE after  DPD: {mse_dpd:.6f}")


# print("Before DPD: min err, max err =", err_no_dpd.min(), err_no_dpd.max())
# print("After  DPD: min err, max err =", err_dpd.min(),    err_dpd.max())


# # --- Plot VIN vs VOUT before and after DPD ---

# # VIN vs VOUT before DPD
# plt.figure()
# plt.scatter(x_test, y_no_dpd, s=5)
# plt.title('PA Input vs Output BEFORE DPD')
# plt.xlabel('v_in')
# plt.ylabel('v_out')
# #plt.show()

# # VIN vs VOUT after DPD
# plt.figure()
# plt.scatter(x_test, y_dpd, s=5)
# plt.title('PA Input vs Output AFTER DPD')
# plt.xlabel('v_in')
# plt.ylabel('v_out')
# plt.show()