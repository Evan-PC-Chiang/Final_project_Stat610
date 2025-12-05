# Final Project

## Description

This model uses a Markov Chain Monte Carlo (MCMC) inversion to estimate the dip-slip rates and locking ratio along the fault, based on a given 2D fault geometry and target surface motion.

## How to Run

### **1. Requirements**

- **MATLAB** (R2020 or newer recommended)

No external libraries are needed beyond standard MATLAB functions.

### **2. Prepare Input Files**

The model requires the following input files:

| **Input**          | **Description**                                             |
| ------------------ | ----------------------------------------------------------- |
| **Surface rates**  | Observed surface motion.                                    |
| **Fault geometry** | 2D fault geometry, discretized into segments, with fault ID |

> Make sure the geometry and rates use consistent coordinate systems and units (km).

### **3. Run the Model**

The workflow consists of two main steps.

#### 1. Build_Greens.m

```matlab
filename = '../data/2RampGeometries.dat'; % assign files
h = 20; % set thickness of the crust
SegPair = [15 16; 16 1; 17 18; 20 3; 22 4; 4 5; 5 6; 6 7; 8 9; 9 10; 10 11; 11 12; 12 13; 23 13; 14 29; 21 22; 26 27]; % add thrust faults by adjacent fault numbers.
```

#### 2. Do_mcmc_[file_name].m

Choose which kind of model you wish to run.

```matlab
model = 'Elastic'; % Elastic or Plate
kind = 'Shear';  % NoShear or Shear
```

#### 3. Output

A figure of dip-slip rates and locking ratio on each fault segment.

## Model setup

We consider two different lithospheric configurations:

1. **Two-layer model**

   An elastic crust overlying an inviscid mantle, similar to our 3D model. This configuration allows the system to flex in response to loading and fault deformation.

2. **Elastic half-space model**

   A simpler configuration used for comparison, where deformation is confined to an elastic medium without underlying mantle flow.

For the kinematic framework, we evaluate two scenarios:

- **Dislocation-only**

  Slip is accommodated entirely by fault dislocation.

- **Dislocation + shear zones**

  Shear zones are incorporated to simulate the geometry and kinematics of a fault-bend-fold system.

<div class="row mt-3 text-center" style="margin-bottom: 2rem">
<div class="col-sm mt-3 mt-md-0">
<table class="tg" style="margin: auto"><thead>
  <tr>
    <th class="tg-0pky" style="border-bottom: 1px solid white;border-right: 1px solid white;"></th>
    <th class="tg-0pky" style="border-bottom: 1px solid white;border-right: 1px solid white;padding-left: 1rem">Elastic half-space</th>
    <th class="tg-0pky" style="border-bottom: 1px solid white;padding-left: 1rem">Two-layer Plate</th>
  </tr></thead>
<tbody>
  <tr>
    <td class="tg-0pky" style="border-right: 1px solid white;"><b>Fault Dislocation</b></td>
    <td class="tg-0pky" style="border: 1px solid white;">Scenario 1</td>
    <td class="tg-0pky">Scenario 3</td>
  </tr>
  <tr>
    <td class="tg-0pky" style="border-right: 1px solid white;border-top: 1px solid white;"><b>Fault-bend-fold</b></td>
    <td class="tg-0pky">Scenario 2</td>
    <td class="tg-0pky"  style="border-top: 1px solid white;border-left: 1px solid white;">Scenario 4</td>
  </tr>
</tbody>
</table>
</div>
</div>

We use a **Markov Chain Monte Carlo (MCMC)** inversion to estimate slip rates and locking ratios along faults. The inversion is performed using both **geodetic surface rates** and **geological uplift/erosion rates**, derived from fault geometries obtained through thermokinematic modeling.

Our goal is to determine the set of slip rates and locking patterns that best reproduce the observed geodetic and geomorphic rates across Taiwan, while remaining consistent with the underlying fault geometry.

## Example

Scenario 1, Please ignore Invicid Mantle part

![https://evan-pc-chiang.github.io/assets/img/post/2D_model/2RampGeometries_Elastic_NoShear-1400.webp](https://evan-pc-chiang.github.io/assets/img/post/2D_model/2RampGeometries_Elastic_NoShear-1400.webp)

Scenario 2, Please ignore Invicid Mantle part

![](https://evan-pc-chiang.github.io/assets/img/post/2D_model/2RampGeometries_Elastic_Shear-1400.webp)

Scenario 3

![](https://evan-pc-chiang.github.io/assets/img/post/2D_model/2RampGeometries_Plate_NoShear-1400.webp)

Scenario 4

![](https://evan-pc-chiang.github.io/assets/img/post/2D_model/2RampGeometries_Plate_Shear-1400.webp)
