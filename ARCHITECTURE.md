# System Architecture and Code Flow Diagram

## High-Level System Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    6DMHS-IDET System                        │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌────────────┐      ┌──────────┐      ┌──────────────┐   │
│  │  Stage I   │  →   │ Stage II │  →   │  Stage III   │   │
│  │  Sensing   │      │ Optimize │      │ Transmission │   │
│  └────────────┘      └──────────┘      └──────────────┘   │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Three-Stage Protocol Flow

### Stage I: Uplink Sensing
```
IDET Receivers (K users)
    │
    ├─── User 1 ─── Sensing Signal ──┐
    ├─── User 2 ─── Sensing Signal ──┤
    ├─── User 3 ─── Sensing Signal ──┼──→ 6DMHS (B RHS)
    └─── User K ─── Sensing Signal ──┘         │
                                               │
                                        ┌──────▼──────┐
                                        │   Sensing   │
                                        │  Algorithm  │
                                        └──────┬──────┘
                                               │
                                        Angle Info (θ, φ)
```

**Key Files**: 
- `Sensing_Algorithm.m` - FFT-based holographic sensing
- `Sensing_Holographic_HalfWave.m` - Half-wavelength sensing
- `Test_Sensing.m` - Sensing performance evaluation

### Stage II: Orientation Adjustment
```
Angle Information (θ, φ)
         │
         ▼
┌─────────────────────┐
│ Joint Optimization  │
│  - Position (q)     │
│  - Rotation (R)     │
│  - Beamforming (Ψ)  │
└─────────┬───────────┘
          │
    ┌─────┴─────┐
    │           │
┌───▼───┐   ┌───▼────┐
│ Space │   │ Holo-  │
│  Pos  │   │ graphic│
│ Opt.  │   │  BF    │
└───┬───┘   └───┬────┘
    │           │
    └─────┬─────┘
          ▼
    6DMHS Orientation
```

**Key Files**:
- `optimize_system.m` - Joint optimization
- `optimize_spatial_position.m` - Position optimization
- `Orientation_Initial.m` - Initialize 6DMHS orientation
- `Orientation_uniformSpherePoints.m` - Generate discrete positions

### Stage III: Downlink Transmission
```
Optimized 6DMHS Configuration
         │
         ▼
┌──────────────────┐
│ Equivalent CSI   │
│    Estimation    │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Digital BF & PS  │
│   Optimization   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ IDET Transmission│
│  - Data (ID)     │
│  - Energy (EH)   │
└──────────────────┘
```

**Key Files**:
- `generate_effective_channel.m` - Equivalent channel
- `calculate_EH_power.m` - Energy harvesting calculation
- `Channel_Generation_init.m` - Channel model

## Module Dependency Diagram

```
                    RUN_OPT_Protocol.m
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
        ▼                   ▼                   ▼
┌───────────────┐   ┌──────────────┐   ┌──────────────┐
│   Sensing     │   │ Optimization │   │ Transmission │
│    Module     │   │    Module    │   │    Module    │
└───────┬───────┘   └──────┬───────┘   └──────┬───────┘
        │                  │                   │
        ▼                  ▼                   ▼
  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
  │ Sensing_    │    │ optimize_   │    │ calculate_  │
  │ Algorithm.m │    │ system.m    │    │ EH_power.m  │
  └─────┬───────┘    └─────┬───────┘    └─────┬───────┘
        │                  │                   │
        ▼                  ▼                   ▼
  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
  │ Sensing_    │    │ optimize_   │    │ generate_   │
  │ Holographic │    │ spatial_    │    │ effective_  │
  │ _HalfWave.m │    │ position.m  │    │ channel.m   │
  └─────────────┘    └─────────────┘    └─────────────┘
```

## Initialization and Support Modules

```
┌──────────────────────────────────────────────────────┐
│              Initialization Layer                     │
├──────────────────────────────────────────────────────┤
│                                                       │
│  ┌─────────────────┐      ┌──────────────────┐     │
│  │ Orientation_    │      │ Channel_         │     │
│  │ Initial.m       │      │ Generation_      │     │
│  │                 │      │ init.m           │     │
│  └────────┬────────┘      └─────────┬────────┘     │
│           │                         │               │
│           ▼                         ▼               │
│  ┌─────────────────┐      ┌──────────────────┐     │
│  │ Orientation_    │      │ Rician Channel   │     │
│  │ uniform         │      │ Model            │     │
│  │ SpherePoints.m  │      │ (LoS + NLoS)     │     │
│  └─────────────────┘      └──────────────────┘     │
│                                                      │
└──────────────────────────────────────────────────────┘
```

## Data Flow Through the System

```
1. System Parameters
   ├─ fc = 30 GHz
   ├─ M = 32×32 = 1024
   ├─ K = 4 users
   ├─ B = 4 RHS
   └─ Pt = 40 dBm
         │
         ▼
2. Initialization
   ├─ 6DMHS positions (q)
   ├─ Initial rotations (R)
   └─ Channel generation (h)
         │
         ▼
3. Sensing Stage
   ├─ Receive uplink signals
   ├─ Holographic image processing
   └─ Extract angles (θ, φ)
         │
         ▼
4. Optimization Stage
   ├─ Optimize positions
   ├─ Design beamforming
   └─ Align max gain direction
         │
         ▼
5. Transmission Stage
   ├─ Estimate equivalent CSI
   ├─ Optimize digital BF
   └─ Calculate IDET metrics
         │
         ▼
6. Performance Evaluation
   ├─ EH power
   ├─ Throughput
   └─ Compare benchmarks
```

## Testing and Visualization Pipeline

```
Main Protocol (RUN_OPT_Protocol.m)
         │
         ├─→ RUN_OPT_Protocol_with_Pt.m
         │        │
         │        └─→ plot_EH_vs_Pt.m ──→ Fig. 5
         │
         ├─→ Test_Sensing.m ──→ Fig. 7-8
         │
         ├─→ Test_6DMA.m ──→ Performance Comparison
         │
         └─→ fig5.m ──→ Generate Paper Figures
```

## File Categories

### 📁 Core Protocol Files (3 files)
```
RUN_OPT_Protocol.m              - Main protocol
RUN_OPT_Protocol_with_Pt.m      - Protocol with Pt sweep
calculate_EH_power.m            - EH power calculation
```

### 🔍 Sensing Files (7 files)
```
Sensing_Algorithm.m             - Core sensing algorithm
Sensing_Holographic.m           - Basic holographic sensing
Sensing_Holographic_HalfWave.m  - Half-wavelength sensing
Sensing_Holographic_HalfWave1.m - Variant implementation
Run_Sensing_HalfWavelength.m    - Run half-wavelength test
Test_Sensing.m                  - Sensing tests
Benchmark_Sensing.m             - Benchmark methods
```

### ⚙️ Optimization Files (6 files)
```
optimize_system.m               - Joint optimization
optimize_spatial_position.m     - Position optimization
Orientation_Initial.m           - Initialize orientation
Orientation_Rot.m               - Rotation operations
Orientation_calculateRotationAngles.m - Angle calculation
Orientation_uniformSpherePoints.m - Generate positions
```

### 📡 Channel & Transmission Files (3 files)
```
Channel_Generation_init.m       - Channel generation
generate_effective_channel.m    - Equivalent channel
fun_run_R00.m                   - Throughput related
```

### 📊 Testing & Plotting Files (9 files)
```
Test_6DMA.m                     - System tests
test_6DMA_run.m                 - Test runner
test_main.m                     - Main tests
test1.m, test2.m, test.m        - Various tests
plot_EH_vs_Pt.m                 - Plot EH vs Pt
plot_EH_vs_Pt_simple.m          - Simple plot
plot_EH_vs_Pt_final.m           - Final plot
plot_EH_vs_Pt_global.m          - Global plot
fig5.m                          - Generate Fig. 5
main1.m                         - Main test program
```

### 🔧 Utility Files (4 files)
```
sig.m, p.m, vec.m              - Simple utilities
untitled.m                      - Unnamed test
```

## Paper-to-Code Mapping Summary

| Paper Section | Code Module | Main Files |
|--------------|-------------|-----------|
| Section II-A | System Model | `Orientation_Initial.m` |
| Section II-B | Channel Model | `Channel_Generation_init.m` |
| Section II-C | Signal Model | `calculate_EH_power.m` |
| Section III-B | Three-Stage Protocol | `RUN_OPT_Protocol.m` |
| Section IV-A | Holographic Sensing | `Sensing_Algorithm.m` |
| Section IV-B | Orientation Optimization | `optimize_system.m` |
| Section V | Transmission Optimization | `generate_effective_channel.m` |
| Section VI | Simulation Results | Test & Plot files |

## Execution Flow (Simplified)

```
START
  │
  ├─→ [1] Initialize System Parameters
  │      (Orientation_Initial.m)
  │
  ├─→ [2] Generate Wireless Channels
  │      (Channel_Generation_init.m)
  │
  ├─→ [3] Perform Holographic Sensing
  │      (Sensing_Algorithm.m)
  │
  ├─→ [4] Optimize 6DMHS Configuration
  │      (optimize_system.m)
  │
  ├─→ [5] Calculate IDET Performance
  │      (calculate_EH_power.m)
  │
  └─→ [6] Evaluate & Visualize Results
         (plot_EH_vs_Pt.m)
  │
END
```

## Key Algorithms Implementation

### Algorithm 1: Holographic Sensing
```
File: Sensing_Algorithm.m
Input: Uplink sensing signals
Process:
  1. Receive holographic image
  2. Apply FFT transformation
  3. Peak detection
Output: Estimated angles (θ, φ)
```

### Algorithm 2: Alternating Optimization
```
File: optimize_system.m
Input: Sensing angles, channel info
Process:
  1. Fix Ψ, optimize (q, R)
  2. Fix (q, R), optimize Ψ
  3. Iterate until convergence
Output: Optimal configuration
```

### Algorithm 3: Fractional Programming
```
File: calculate_EH_power.m
Input: Optimized 6DMHS config
Process:
  1. Fix ρ, optimize X
  2. Fix X, optimize ρ
  3. Iterate until convergence
Output: Optimal EH power
```

## Performance Metrics Calculated

```
┌────────────────────────────────┐
│     Performance Metrics         │
├────────────────────────────────┤
│                                 │
│  • EH Power (P_EH)             │
│    - RF power (Eq. 20)         │
│    - DC power (Eq. 21)         │
│                                 │
│  • Throughput (R)              │
│    - SINR (Eq. 18)             │
│    - Rate (Eq. 19)             │
│                                 │
│  • Sensing Accuracy            │
│    - RMSE of angles            │
│    - CDF analysis              │
│                                 │
│  • Beamforming Gain            │
│    - Directional gain          │
│    - Alignment accuracy        │
│                                 │
└────────────────────────────────┘
```

## Quick Navigation Guide

**Want to understand sensing?** 
→ Start with `Sensing_Algorithm.m`

**Want to understand optimization?** 
→ Start with `optimize_system.m`

**Want to run the system?** 
→ Start with `RUN_OPT_Protocol.m`

**Want to see results?** 
→ Check `plot_EH_vs_Pt.m` and test files

**Want to reproduce paper figures?** 
→ Run corresponding test scripts (see README.md)

---

*This diagram document complements the detailed code analysis documents (CODE_ANALYSIS.md and 代码分析文档.md) by providing visual representations of the system architecture and code relationships.*
