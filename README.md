# 🌀 Galerkin Reduced-Order Model (G-ROM) for 2D Navier–Stokes Equations

This repository provides a MATLAB implementation of a **Galerkin Reduced-Order Model (G-ROM)** for the two-dimensional incompressible Navier–Stokes equations

The code demonstrates how to construct and evolve a low-dimensional model using **POD (Proper Orthogonal Decomposition)** basis functions extracted from high-fidelity DNS data. It compares the G-ROM solution against the full-order model (FOM / DNS) in terms of kinetic energy and $L^2$ error.






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
