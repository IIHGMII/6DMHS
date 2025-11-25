# Code Analysis Document for 6DMHS System

## Overview

This document provides a detailed analysis of the MATLAB code implementation related to the paper "6D Movable Holographic Surface Assisted Integrated Data and Energy Transfer: A Sensing Enhanced Approach". The system implements a 6D Movable Holographic Surface (6DMHS) based Integrated Data and Energy Transfer (IDET) system.

## Core Concepts from the Paper

The paper proposes a three-stage protocol:
1. **Uplink Sensing Stage**: Acquire angular information of IDET receivers using holographic sensing technology
2. **Orientation Adjustment Stage**: Optimize the spatial positions and holographic beamforming of 6DMHS
3. **Downlink Transmission Stage**: Optimize digital beamforming and power allocation for IDET transmission

## Quick Reference Tables

### Core Files Quick Reference

| File Name | Category | Paper Section | Main Function |
|-----------|----------|--------------|--------------|
| `RUN_OPT_Protocol.m` | Main Flow | Section III-B | Complete three-stage protocol execution |
| `Sensing_Algorithm.m` | Sensing | Section IV-A | Holographic sensing and angle estimation |
| `optimize_system.m` | Optimization | Section IV-B | Joint optimization of position and beamforming |
| `calculate_EH_power.m` | Transmission | Section II-C | Calculate energy harvesting power |
| `Channel_Generation_init.m` | Channel | Section II-B | Rician channel generation |
| `Orientation_Initial.m` | Initialization | Section II-A | 6DMHS orientation initialization |

### Testing and Plotting Files Quick Reference

| File Name | Corresponding Figure | Function |
|-----------|---------------------|----------|
| `plot_EH_vs_Pt.m` | Fig. 5 | EH power vs transmit power |
| `Test_Sensing.m` | Fig. 7-8 | Sensing performance testing |
| `Test_6DMA.m` | Fig. 5 | Performance comparison of different configurations |

## Code File Classification and Functions

### I. Main Control Flow Files

#### 1. `RUN_OPT_Protocol.m`
**Function**: Main protocol execution file implementing the complete three-stage protocol

**Corresponding Paper Section**: Section III-B (Sensing-Enhanced 6DMHS Protocol)

**Key Features**:
- Initialize system parameters (frequency fc=30GHz, antenna array 32×32, K=4 users, B=4 RHS)
- Execute the complete three-stage protocol flow
- Call various sub-modules for sensing, optimization, and transmission

**Key Implementation**:
```matlab
% Parameter initialization
fc = 30e9;  % Carrier frequency
Mx = 32; My = 32; M = Mx * My;  % Antenna array size
K = 4;  % Number of users
B = 4;  % Number of holographic surfaces
```

**Corresponding Paper Equations**: Implements the system model in Section II

#### 2. `RUN_OPT_Protocol_with_Pt.m`
**Function**: Similar to `RUN_OPT_Protocol.m` but adds performance analysis under different transmit powers Pt

**Corresponding Paper Section**: Section VI (Simulation Results), especially Fig. 5

**Key Features**:
- Run the protocol at different transmit power levels
- Generate curves showing the relationship between transmit power and EH power

### II. Sensing Module (Stage I: Uplink Sensing)

#### 3. `Sensing_Algorithm.m`
**Function**: Core sensing algorithm implementation based on holographic sensing method for estimating user angle information

**Corresponding Paper Section**: Section IV-A (Holographic-based Sensing Method)

**Key Features**:
- Implement FFT-based holographic image processing
- Extract angular information (azimuth θ and elevation φ) of IDET receivers
- Process sensing signals in local coordinate system

**Key Implementation**:
```matlab
% FFT matrix for angle estimation
% N: Number of FFT sampling points (e.g., 1024)
% delta_n: Sample index vector from -(N-1)/2 to (N-1)/2
delta_n = [(-N + 1) / 2:1:(N - 1) / 2].';
% F: FFT transformation matrix for 2D angle extraction
F = exp(-1j * 2*pi * (delta_n * delta_n.') / N);
```

**Corresponding Paper Equations**:
- Eq. (27)-(29): Holographic sensing received signal model
- Eq. (30): FFT-based detection method
- Theorem 1: Sensing accuracy analysis

#### 4. `Sensing_Holographic_HalfWave.m` / `Sensing_Holographic_HalfWave1.m`
**Function**: Implement holographic sensing with half-wavelength spacing

**Corresponding Paper Section**: Lemma 2 (on the sufficiency of half-wavelength spacing)

**Key Features**:
- Deploy sensing elements at half-wavelength (λ/2) intervals
- Verify that half-wavelength spacing is sufficient for high-precision sensing

**Corresponding Paper Equation**: Eq. (94): d_S ≤ λ/2

#### 5. `Run_Sensing_HalfWavelength.m`
**Function**: Run and test half-wavelength sensing method

#### 6. `Test_Sensing.m`
**Function**: Test script for sensing algorithms

**Key Features**:
- Test performance of different sensing methods
- Compare holographic sensing with traditional LS methods

**Corresponding Paper Content**: Sensing performance comparison in Fig. 7-8

#### 7. `Benchmark_Sensing.m`
**Function**: Benchmark sensing method implementation (e.g., LS-based sensing)

**Corresponding Paper Content**: Baseline methods for performance comparison

### III. Optimization Module (Stage II: Orientation Adjustment)

#### 8. `optimize_system.m`
**Function**: Joint optimization of spatial positions and holographic beamforming

**Corresponding Paper Section**: Section IV-B (Orientation Optimization)

**Key Features**:
- Optimize spatial positions (translations) of 6DMHS
- Optimize holographic beamforming
- Adjust RHS orientation based on sensing information

**Key Implementation**:
```matlab
% Generate discrete spatial points on a sphere
% B1: Number of discrete spatial positions (e.g., 40)
[X, Y, Z] = Orientation_uniformSpherePoints(B1);

% Calculate steering vectors for each position
% k0: RHS index, b0: position index, k1: user index
% lambda: wavelength, a1: element positions after rotation
% e(:, k1): direction vector to user k1
for k0 = 1:K
    for b0 = 1:B1
        for k1 = 1:K
            a1 = q(b0, :) + (Rot_Matrix{k0} * Coor_Ele_init.').';
            a_st_set{k0, b0, k1} = exp(1j * 2*pi/lambda * (a1 * e(:, k1)));
        end
    end
end
```

**Corresponding Paper Equations**:
- Problem (P2): Orientation optimization problem
- Eq. (34)-(35): Holographic beamforming design
- Algorithm 1 & 2: Alternating optimization algorithms

#### 9. `optimize_spatial_position.m`
**Function**: Specifically optimize spatial positions (translation part) of 6DMHS

**Corresponding Paper Content**: Position optimization sub-problem

#### 10. `Orientation_Initial.m`
**Function**: Initialize spatial orientation of 6DMHS

**Key Features**:
- Uniformly distribute B RHS on a spherical surface
- Calculate initial rotation matrices (Rodrigues rotation formula)
- Set initial normal vectors

**Corresponding Paper Equations**:
- Eq. (1)-(6): 6DMHS coordinate transformation
- Eq. (2): Rodrigues rotation matrix

#### 11. `Orientation_uniformSpherePoints.m`
**Function**: Uniformly generate discrete position points on a sphere

**Key Features**:
- Generate uniformly distributed spherical coordinates
- Discretize the feasible position set of 6DMHS

#### 12. `Orientation_calculateRotationAngles.m`
**Function**: Calculate rotation angles

**Corresponding Paper Equations**: Eq. (4)-(5): Rotation vector and rotation angle calculation

#### 13. `Orientation_Rot.m`
**Function**: Execute rotation operations

### IV. Transmission and Performance Evaluation Module (Stage III: Downlink Transmission)

#### 14. `calculate_EH_power.m`
**Function**: Calculate Energy Harvesting (EH) power

**Corresponding Paper Section**: Section II-C (Downlink Transmission Signal Model)

**Key Features**:
- Calculate EH power based on optimized configuration
- Support parameterized transmit power input

**Corresponding Paper Equations**:
- Eq. (20): RF energy harvesting power
- Eq. (21): DC energy harvesting power (nonlinear EH model)

#### 15. `generate_effective_channel.m`
**Function**: Generate equivalent channel

**Key Features**:
- Calculate equivalent channel after holographic beamforming
- Equivalent channel = wireless channel × holographic beamforming × electromagnetic response

**Corresponding Paper Content**: "equivalent CSI" mentioned in Stage III

#### 16. `Channel_Generation_init.m`
**Function**: Initialize wireless channel generation

**Corresponding Paper Section**: Section II-B (Channel Model)

**Key Features**:
- Generate Rician channel (LoS + NLoS)
- Calculate path loss
- Generate scatterer parameters

**Corresponding Paper Equations**:
- Eq. (10)-(13): Wireless channel model
- Eq. (11): Rician channel gain

### V. Testing and Result Generation Files

#### 17. `plot_EH_vs_Pt.m` / `plot_EH_vs_Pt_simple.m` / `plot_EH_vs_Pt_final.m` / `plot_EH_vs_Pt_global.m`
**Function**: Generate different versions of EH power vs. transmit power plots

**Corresponding Paper Content**: Fig. 5 - EH power versus transmit power

#### 18. `fig5.m`
**Function**: Specifically generate Figure 5 from the paper

#### 19. `Test_6DMA.m` / `test_6DMA_run.m`
**Function**: Test overall performance of 6DMA system

**Key Features**:
- Compare different configurations (FPA, rotation only, translation only, full 6DMHS)
- Verify performance gains of 6DMHS

**Corresponding Paper Content**: Comparison of different configurations in Fig. 5

#### 20. `test_main.m` / `test1.m` / `test2.m`
**Function**: Various test scripts

#### 21. `main1.m`
**Function**: One of the main test programs, possibly for complete system testing

### VI. Utility Tool Files

#### 22. `fun_run_R00.m`
**Function**: Utility function, possibly related to throughput threshold R0

**Corresponding Paper Content**: Throughput constraints in Problem (P1)-(P2)

#### 23. `sig.m` / `p.m` / `vec.m`
**Function**: Simple utility functions

#### 24. `untitled.m` / `untitled.fig`
**Function**: Unnamed test files or figure files

## Code-to-Paper Mapping

### Three-Stage Protocol Implementation Mapping

| Protocol Stage | Paper Section | Main Code Files |
|---------------|--------------|----------------|
| Stage I: Uplink Sensing | Section IV-A | `Sensing_Algorithm.m`, `Sensing_Holographic_HalfWave.m` |
| Stage II: Orientation Adjustment | Section IV-B | `optimize_system.m`, `optimize_spatial_position.m`, `Orientation_Initial.m` |
| Stage III: Downlink Transmission | Section V | `calculate_EH_power.m`, `generate_effective_channel.m` |

### Key Technology Implementation Mapping

| Technology | Paper Section | Code Implementation |
|-----------|--------------|-------------------|
| Holographic Sensing | Section IV-A, Eq. (27)-(32) | `Sensing_Algorithm.m` |
| FFT Detection | Theorem 1 | FFT matrix in `Sensing_Algorithm.m` |
| Half-wavelength Spacing | Lemma 2 | `Sensing_Holographic_HalfWave.m` |
| Rodrigues Rotation | Eq. (2)-(5) | `Orientation_Initial.m` |
| Holographic Beamforming | Eq. (22), (34)-(35) | `optimize_system.m` |
| Position Optimization | Problem (P3) | `optimize_spatial_position.m` |
| EH Power Calculation | Eq. (20)-(21) | `calculate_EH_power.m` |

### Simulation Results Generation Mapping

| Figure | Paper Location | Code Files |
|--------|---------------|-----------|
| Fig. 5: EH vs Pt | Section VI | `plot_EH_vs_Pt.m`, `RUN_OPT_Protocol_with_Pt.m` |
| Fig. 6: EH vs R0 | Section VI | Related to `fun_run_R00.m` |
| Fig. 7-8: Sensing performance | Section VI | `Test_Sensing.m`, `Benchmark_Sensing.m` |

## System Parameter Settings

Key parameters in the code are consistent with the paper:

| Parameter | Code Value | Paper Value | Description |
|-----------|-----------|------------|-------------|
| fc | 30e9 | 30 GHz | Carrier frequency |
| M | 32×32 = 1024 | 32×32 | Number of RHS elements |
| K | 4 | 4 | Number of users |
| B | 4 | 4 | Number of RHS |
| Q | 1 | 1 | Number of feeds |
| K_rice | 10 | 10 | Rician factor |
| Pt | 40 dBm | 40 dBm | Transmit power |
| d | λ/4 | λ/4 | Element spacing |

## Code Execution Flow

Complete code execution flow:

1. **Initialization** (`Orientation_Initial.m`)
   - Set system parameters
   - Initialize 6DMHS spatial orientation
   - Generate initial coordinates

2. **Channel Generation** (`Channel_Generation_init.m`)
   - Generate Rician wireless channel
   - Calculate path loss
   - Generate scatterer information

3. **Sensing Stage** (`Sensing_Algorithm.m`)
   - Receive uplink sensing signals
   - Process holographic image with FFT
   - Estimate user angular information

4. **Optimization Stage** (`optimize_system.m`)
   - Optimize RHS positions based on sensing information
   - Design holographic beamforming
   - Adjust 6DMHS orientation

5. **Transmission Stage** (`calculate_EH_power.m`, `generate_effective_channel.m`)
   - Generate equivalent channel
   - Optimize digital beamforming
   - Calculate IDET performance metrics

6. **Performance Evaluation** (various test and plot files)
   - Compare different configurations
   - Generate simulation results
   - Plot performance curves

## Code Implementation of Core Innovations

### 1. Holographic Sensing Technology
**Paper Innovation**: Propose parallel sensing method based on holographic images, avoiding signal aliasing in traditional serial methods

**Code Implementation**: `Sensing_Algorithm.m` uses FFT matrix to directly extract angle information from holographic images
```matlab
F = exp(-1j * 2*pi * (delta_n * delta_n.') / N);
% Obtain data directly from sensing elements in parallel, avoiding aliasing
```

### 2. Utilization of Directional Beamforming Gain
**Paper Innovation**: Identify and exploit the directional beamforming gain property of RHS

**Code Implementation**: `optimize_system.m` aligns the maximum beamforming gain direction with users

### 3. Sensing-Enhanced Three-Stage Protocol
**Paper Innovation**: Use sensing information to reduce algorithm complexity from O(M²N²L²)

**Code Implementation**: Two-step method of sensing then optimization, `RUN_OPT_Protocol.m` coordinates the entire process

## Summary

This codebase completely implements the 6DMHS-IDET system proposed in the paper, including:

1. **Complete three-stage protocol**: Sensing, optimization, transmission
2. **Innovative holographic sensing method**: FFT-based parallel sensing
3. **Joint optimization algorithms**: Joint design of position, rotation, and beamforming
4. **Performance evaluation framework**: Comparison with multiple baseline schemes

The code structure is clear and well-modularized, with each file corresponding to specific sections or algorithms from the paper, making it easy to understand and reproduce the paper's results.

## Recommended Code Reading Order

1. **Start with protocol flow**: `RUN_OPT_Protocol.m` - Understand the overall process
2. **Then initialization**: `Orientation_Initial.m`, `Channel_Generation_init.m` - Understand system model
3. **Deep dive into sensing**: `Sensing_Algorithm.m` - Understand core innovation
4. **Study optimization**: `optimize_system.m` - Understand joint optimization methods
5. **Finally testing**: Various test and plot files - Understand performance verification

Following this order allows for a systematic understanding of the correspondence between the paper's theory and code implementation.
