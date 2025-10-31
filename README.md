# Multiple Boilers on a Common Header

## 📖 Project Overview
This repository contains my final year project for the Department of Electrical Engineering, University of Cape Town.  
The project investigates Multiple boilers on a common header.  

While this configuration is widely used in industry for flexibility and redundancy, it introduces strong dynamic coupling: a change in one boiler’s operation directly affects the header pressure and the other boilers. Without proper coordination, this can lead to uneven load sharing, oscillations, and instability.

---

## 🎯 Objectives
- Develop a dynamic model of two drum boilers and a common header in MATLAB/Simulink.  
- Capture key thermohydraulic phenomena such as shrink and swell and pressure–flow coupling.  
- Design and implement a cascade control architecture that ensures:
  - Stable header pressure  
  - Safe drum water levels  
  - Fair load sharing between boilers  
- Validate the control system under different load conditions (100%, 80%, 20%).  
- Demonstrate feasibility for PLC implementation using IEC 61131 standards.

---

## 🛠️ Methodology
- **Plant Modelling:**  
  - First‑principles modelling of drum boilers and header dynamics.  
  - Inclusion of asymmetrical pipeline resistances to reflect real industrial layouts.  

- **Control Design:**  
  - Fast PI loop for header pressure regulation.  
  - Supervisory allocator for distributing total steam demand.  
  - Ratio PI controller for enforcing proportional load sharing.  
  - Dedicated PI controllers for drum water level control.  

- **Simulation Environment:**  
  - MATLAB/Simulink for plant and controller implementation.  
  - Closed‑loop tests for setpoint tracking and disturbance rejection.  

---

## 📊 Key Results
- At **100% and 80% load**, the cascade strategy achieved:
  - Accurate header pressure regulation  
  - Stable drum water levels  
  - Balanced load sharing  

- At **20% load**, the ratio controller struggled, leading to imbalance — highlighting the limitations of fixed‑gain controllers in nonlinear regimes.  

- Poorly tuned controllers reproduced **limit cycles**, confirming the system’s sensitivity and validating the need for careful tuning.  

---

## 🚀 Future Work
- Adaptive control strategies for low‑load operation.  
- State estimation to improve decision‑making.  
- Migration towards **Model Predictive Control (MPC)** for wide‑range robustness.  


---

How to run simulations:
The two important files found in each folder are the Init_two_Boilers.m/INIT_BOILER.m and Single_Boiler_Model.slx/Multi_Boiler_Common_header.slx.

You need to first open both the .m and the .slx files first.
The run the .m init file for that respective .slx 
e.g if I want to run the Single_Boiler_Model.slx simulation, i need to first run the INIT_BOILERtp initialise and assign the model parameters.
From there you should be good to go.
---

## 👤 Author
**Zwivhuya Ndou**  
Final Year Project – BSc Electrical & Computer Engineering  
University of Cape Town (2025)  

Supervisor: **Prof. Edward Boje**
