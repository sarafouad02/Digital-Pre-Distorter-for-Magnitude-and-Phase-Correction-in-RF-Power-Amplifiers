# Digital Pre‑Distortion (DPD) Project

This repository contains two main phases of a Digital Pre‑Distortion system for RF Power Amplifiers:

1. **High‑Level Modeling Phase (Python)**
2. **RTL Implementation Phase (SystemVerilog)**

---

## 1. High‑Level Modeling Phase

### Overview

The high‑level model characterizes amplifier non‑linearity using a memory‑polynomial (Volterra‑series) approach, solves for complex correction coefficients via least‑squares, and exports them for hardware realization.

### Directory Structure

```
/high_level/
├── dpd_model.py         # Python script implementing PA model, regressor, LS solver
├── dpd_coeffs.json      # Output JSON with computed coefficients
├── dpd_real.txt         # (optional) golden real outputs for verification
├── dpd_imag.txt         # (optional) golden imag outputs for verification
└── requirements.txt     # Python dependencies (numpy, matplotlib)
```

### Prerequisites

* Python 3.8+
* `numpy`
* `matplotlib`
* `json`

Install dependencies:

```bash
pip install -r high_level/requirements.txt
```

### Usage

1. **Generate coefficients**:

   ```bash
   python high_level/dpd_model.py
   ```

   * Solves for DPD coefficients and writes `dpd_coeffs.json`.

2. **Visualize constellations** (optional):

   * The script plots three scatter plots: PA output without DPD, with Python‑DPD, and (if available) RTL‑DPD.

---

## 2. RTL Implementation Phase

### Overview

The RTL design ingests IQ samples, buffers past taps, applies the trained coefficients via a MAC‑array, and streams out pre‑distorted samples. Includes a self‑checking testbench against the Python golden‑model.

### Directory Structure

```
/rtl/
├── sample_buffer.sv    # Circular buffer for M+1 taps
├── dpd_mac_array.sv    # MAC array with CORDIC magnitude and polynomial terms
├── coeff_rom.sv        # Dual‑port ROM loading coeffs from .mem files
├── dpd_top.sv          # Top-level integration module
├── dpd_top_tb.sv       # Self‑checking SystemVerilog testbench
├── coeffs_real.mem     # Fixed‑point real coefficients (hex)
├── coeffs_imag.mem     # Fixed‑point imag coefficients (hex)
└── syn/                # Synthesis & constraint scripts
    ├── dpd_top.sdc     # SDC constraints
    ├── design_dc.tcl   # Synopsys Design Compiler script
    └── libraries/      # Standard‑cell .db files
```

### Prerequisites

* Synopsys Design Compiler
* Standard‑cell libraries (TSMC 130 nm)
* Vivado, QuestaSim or ModelSim for simulation

### Simulation

```bash
# Compile RTL
vlog rtl/*.sv
# Run testbench
vsim dpd_top_tb
# Observe PASS messages and waveform
```

### Coefficient Conversion

Convert JSON to `.mem` files (for coeff\_rom):

```bash
python scripts/json_to_mem.py high_level/dpd_coeffs.json rtl/coeffs_real.mem rtl/coeffs_imag.mem
```

### Synthesis & Timing

```bash
# Launch Design Compiler
dc_shell -f syn/design_dc.tcl
```

* Top module: `dpd_top`
* Clock period: 10 ns (100 MHz)
* Reports: `area.rpt`, `power.rpt`, `setup.rpt`, `hold.rpt`, `constraints.rpt`

---


