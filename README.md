# POLISOC: A Voltage-Deformation Hybrid SOC Estimation Algorithms Platform for Lithium-Ion Batteries

**POLISOC** is an open-source algorithms platform implemented in MATLAB for hybrid State of Charge (SOC) estimation and empirical modelling of lithium-ion batteries.

The key contribution of the platform is the inclusion of mechanical measurements (thickness change) in the SOC estimation, and the mechanical empirical modelling of lithium-ion batteries.  
The algorithm leverages the SOC-thickness change (SOC-THK) relationship, which is almost linear for most chemistries used in commercial lithium-ion batteries, and it is largely unaffected by the intensity of the applied current.  
These properties make deformation (thickness change) an interesting metric for SOC estimation.

The platform includes:
- A set of codes for SOC estimation and empirical modelling, contained in the `Codes` folder.
- A set of experimental data to run and test the algorithms and the empirical models, contained in the `Data` folder.
- Voltage and deformation interpolation functions, contained in the `Utility function` folder.

---

## 📁 Data Structures

The data provided were experimentally measured at the DIMEAS lab at Politecnico di Torino. Battery chemistries include **LCO**, **LFP**, and **NMC**.

### Characterization Tests
Tests such as constant current charge/discharge cycles and hybrid pulse power characterization (HPPC) have been carried out to measure:
- Capacity
- SOC-OCV and SOC-THK relationships
- Hysteresis parameters
- Resistance and capacitance values

Each battery sample's parameters are saved in the `Data` folder in a MATLAB structure named: param_"chemistry""#sample", e.g. param_NMC1.mat


#### Structure fields:
- `param.DthkD{num}`: Thickness change during full discharge cycle at 1C as a function of `param.SOC` [mm]
- `param.DthkC{num}`: Thickness change during full charge cycle at C/2 as a function of `param.SOC` [mm]
- `param.OCV`: Open-circuit voltage vs. SOC [V]
- `param.Q`: Cell capacity (1C constant current) [Ah]
- `param.R0`: Ohmic resistance [Ω]
- `param.R1`: First-order resistance [Ω]
- `param.C1`: First-order capacitance [F]
- `param.Me`: Voltage hysteresis amplitude vs. SOC [V]
- `param.Ge`: Voltage hysteresis decay rate [-]
- `param.Mm`: Deformation hysteresis amplitude vs. SOC [mm]
- `param.Gm`: Deformation hysteresis decay rate [-]

> The variable `{num}` refers to the reference performance test number (e.g., 1 → BOL, 44 → EOL).  
> Valid only for LCO batteries.

### Test Profiles

Testing profiles such as **Dynamic Stress Test (DST)** and **drive cycles** are also provided.

These are saved as MATLAB structures named: Meas_"chemistry""#sample""TestType""#test", e.g. Meas_NMC1_DST1.mat



#### Structure fields:
- `Meas.Time`: Time [s]
- `Meas.Voltage`: Recorded voltage [V]
- `Meas.TrueCurrent`: Actual battery current [A]
- `Meas.Current`: Sensor-measured current (may include bias) [A]
- `Meas.Deformation`: Battery thickness change [mm]
- `Meas.Temperature`: Battery temperature [°C]
- `Meas.Strain`: Thickness strain (if available) [mm/mm]

---

## ⚙️ SOC Estimation Algorithms

The platform includes two SOC estimation algorithms: **Deformation Inversion** and **Hybrid EKF**.

### 1. Deformation Inversion
This method estimates SOC by inverting the SOC-THK relationship.

To run:
1. Open `RunInversion.m` in the `Deformation inversion` folder.
2. Set:
   - Battery chemistry with the `BatType` variable.
   - Battery sample parameters (`param` structure)
   - Test profile (`Meas` structure)
3. Run the code.

> ⚠️ **Note**: For **LFP batteries**, the non-monotonic SOC-THK curve (especially between 30%-60% SOC) makes estimation unreliable in that range.

### 2. Hybrid EKF
An **Extended Kalman Filter (EKF)**-based SOC estimator that supports:

- **Voltage-based mode**: Traditional SOC estimation using voltage.
- **Deformation-based mode**: SOC estimation using deformation (thickness change).
- **Hybrid mode**: Combines both voltage and deformation measurements.

To run:
- Open `main.m` in the `Hybrid ekf` folder.
- Set:
  - Battery chemistry with the `BatType` variable.
  - Battery sample parameters (`param` structure)
  - Test profile (`Meas` structure)
  - Choose the mode (voltage-based, deformation-based, or hybrid)

> ⚡ You can also use the version in `Hybrid EKF_No Mechanical Hysteresis`, which is faster and excludes mechanical hysteresis (has minimal effect on results).

---

## 🔋 Equivalent Circuit Battery Models

The platform includes **empirical equivalent circuit models** for both electrical and mechanical responses:

### 1. Electrical Circuit Model (ECM)
The model calculates the **voltage** given a current profile.

To run:
- Open `ECM.m` in the `Electrical equivalent circuit` folder.
- Set:
  - Battery chemistry with the `BatType` variable.
  - Battery sample parameters (`param` structure)
  - Test profile (`Meas` structure)
  - Choose the mode (voltage-based, deformation-based, or hybrid)


### 2. Mechanical Circuit Model (MCM)
The model calculates the **deformation** (thickness change) given a current profile.

To run:
- Open `MCM.m` in the `Mechanical equivalent circuit` folder.
- Set:
  - Battery chemistry with the `BatType` variable.
  - Battery sample parameters (`param` structure)
  - Test profile (`Meas` structure)
  - Choose the mode (voltage-based, deformation-based, or hybrid)

---

## Summary

POLISOC integrates mechanical measurements into battery state estimation in a novel and efficient way.  
By providing access to both the models and experimental datasets, it enables users to experiment with advanced hybrid SOC estimation methods that go beyond traditional voltage-based approaches.

---

## Requirements
MATLAB 2006a and later

---

## References
More details on the implementation are included in the relevant paper "POLISOC: A Voltage-Deformation Hybrid SOC Estimation Algorithm for Lithium-Ion Batteries", currently in peer review. 


