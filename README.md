# 🌀 Galerkin Reduced-Order Model (G-ROM) for 2D Navier–Stokes Equations

This repository provides a MATLAB implementation of a **Galerkin Reduced-Order Model (G-ROM)** for the two-dimensional incompressible Navier–Stokes equations

The code demonstrates how to construct and evolve a low-dimensional model using **POD (Proper Orthogonal Decomposition)** basis functions extracted from high-fidelity DNS data. It compares the G-ROM solution against the full-order model (FOM / DNS) in terms of kinetic energy and $L^2$ error.

---

## 📂 Repository Structure

Galerkin_ROM_NSE-main/
│
├── main_rom_galerkin.m # Main driver script for G-ROM simulation
│
├── ROM_FUN/ # ROM-related utility functions
│ ├── nse_g_rom_online_fun.m # Time integration of reduced model
│ └── other helper routines...
│
├── FEM_FUN/ # FEM-related utilities (e.g., assembly)
│
├── DATA/ # Precomputed data (DNS, POD, ROM matrices)
│ ├── nse_dns_t23_cylmesh35dt002Re_1eN3_T10.mat
│ ├── nse_pod_re1000_snap1500.mat
│ ├── nse_rom_mat_re1000_snap1500.mat
│ └── ke_dns_re1000.mat
│
├── save_required_rom_vars_reduced.m # Script to extract minimal ROM data
└── README.md




---

## ⚙️ Requirements

- **MATLAB R2020b** or later  
- Toolboxes: *none required beyond base MATLAB*  
- System: macOS, Windows, or Linux  

Optional: For large data handling, ensure MATLAB supports the `-v7.3` MAT-file format.

---

## 🚀 Running the G-ROM Simulation

1. **Clone this repository:**
   ```bash
   git clone https://github.com/<yourusername>/Galerkin_ROM_NSE.git
   cd Galerkin_ROM_NSE-main
