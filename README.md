## Memory-Polynomial DPD Model

This repository contains a Python implementation of a memory-polynomial Digital Predistortion (DPD) workflow for Power Amplifiers (PAs). The code covers:

* **PA Simulation**: A memory-polynomial PA model (`pa_model`) capturing AM/AM and AM/PM distortion up to user-defined memory depth (M) and polynomial order (K).
* **DPD Training**: Least-squares estimation to compute predistorter coefficients based on simulated PA output and known input signals.
* **DPD Inference**: Application of the computed predistorter to new signals, followed by PA re-amplification, to evaluate distortion compensation.
* **Performance Metrics & Visualization**: MSE evaluations (baseband, AM/AM, AM/PM) and binned median plots for AM/AM, AM/PM, and power-transfer characteristics.

---

### Features

* **Configurable Parameters**: Easily adjust memory depth `M`, polynomial order `K`, and dataset sizes (`N_train`, `N_test`).
* **Reproducible**: Fixed random seed ensures repeatable training/inference runs.
* **Complex Coefficients**: Supports complex-valued PA coefficients for joint AM/AM and AM/PM behavior.
* **JSON & NumPy Outputs**: Saves estimated DPD coefficients in both `.npy` and human-readable `.json` formats.

---

### Requirements

* Python 3.7+
* NumPy
* Matplotlib

Install dependencies with:

```bash
pip install numpy matplotlib
```

---

### Usage

1. **Configuration**

   * Open `dpd_workflow.py` (or your script file) and set:

     ```python
     M = 2        # Memory depth
     K = 3        # Polynomial order
     N_train = 5000
     N_test  = 2000
     ```
   * Adjust PA coefficient dictionary `b` for your specific amplifier model.

2. **Training Phase**

   ```bash
   python dpd_workflow.py --mode train
   ```

   * Generates `x_train`, simulates `d_train = pa_model(x_train)`, builds predistorter basis, solves least-squares for `a_hat`, and saves:

     * `dpd_coeffs.npy`
     * `dpd_coeffs.json`

3. **Inference Phase**

   ```bash
   python dpd_workflow.py --mode infer
   ```

   * Loads `dpd_coeffs.npy`, simulates a new `x_test`, compares PA output before/after DPD, computes MSE, and plots results.

---

### Script Breakdown

* **`pa_model(x)`**: Memory-polynomial PA simulator
* **Training**:

  * Generate random QAM-like waveform
  * Simulate PA distortion
  * Build basis `X_pred`
  * Solve `X_pred @ a_hat ≈ x_train`
  * Save DPD coefficients
* **Inference**:

  * Reload `a_hat`
  * Generate test waveform
  * Evaluate PA output before (`y_no_dpd`) and after (`y_dpd`) DPD
  * Compute MSE & AM/AM, AM/PM metrics
  * Plot binned medians for visual analysis

---

### Visualization

The final script produces a three-panel figure showing:

1. **AM/AM**: Median output magnitude vs input power
2. **AM/PM**: Median phase error vs input power
3. **Power Transfer**: Output power (dB) vs input power (dB)

These plots help assess predistortion performance across the input power range.

---

### Extending & Customization

* **Change PA Behavior**: Modify or replace `b` coefficients for different amplifier characteristics.
* **Add Noise**: Insert noise into training or inference to simulate realistic channel conditions.
* **Alternate Fitting**: Swap least-squares solver for regularized or weighted approaches.

---

### License

This project is provided under the MIT License. Feel free to use, modify, and distribute.
