# 6DMHS: 6D Movable Holographic Surface for IDET

## Repository Overview

This repository contains the MATLAB implementation of the research paper:

**"6D Movable Holographic Surface Assisted Integrated Data and Energy Transfer: A Sensing Enhanced Approach"**

*Authors: Zhonglun Wang, Yizhe Zhao, Gangming Hu, Yali Zheng, and Kun Yang*

## About the Project

The project implements a novel 6D Movable Holographic Surface (6DMHS) based Integrated Data and Energy Transfer (IDET) system. The system leverages:

- **Holographic-based sensing** for acquiring user angular information
- **6D movable antenna technology** (3 rotations + 3 translations) 
- **Directional beamforming gain optimization** of Reconfigurable Holographic Surface (RHS)
- **Three-stage protocol** for efficient IDET transmission

## Key Innovations

1. **Holographic Sensing Technology**: Novel FFT-based parallel sensing method that avoids signal aliasing
2. **Directional Beamforming Gain Utilization**: Exploits the amplitude-controlled beamforming characteristics of RHS
3. **Sensing-Enhanced Protocol**: Reduces optimization complexity from O(M²N²L²) to practical levels
4. **Joint Optimization**: Simultaneous optimization of position, rotation, and beamforming

## Repository Structure

```
6DMHS/
├── 6D Movable Holographic Surface ... .pdf    # Research paper (PDF)
├── 6D Movable Holographic Surface ... .md     # Research paper (Markdown)
├── CODE_ANALYSIS.md                           # Detailed code analysis (English)
├── 代码分析文档.md                             # Detailed code analysis (Chinese)
├── README.md                                  # This file
└── Proposed-v2/                               # MATLAB implementation
    ├── RUN_OPT_Protocol.m                     # Main protocol execution
    ├── Sensing_Algorithm.m                    # Holographic sensing
    ├── optimize_system.m                      # Joint optimization
    ├── calculate_EH_power.m                   # Energy harvesting calculation
    ├── Channel_Generation_init.m              # Channel model
    ├── Orientation_Initial.m                  # 6DMHS initialization
    └── ... (see detailed analysis documents for complete file descriptions)
```

## Getting Started

### Prerequisites

- MATLAB R2020a or later (recommended)
- Signal Processing Toolbox
- Optimization Toolbox (for CVX, if used)

### Quick Start

1. Clone the repository:
```bash
git clone https://github.com/IIHGMII/6DMHS.git
cd 6DMHS/Proposed-v2
```

2. Run the main protocol:
```matlab
% In MATLAB
[P_out, Beam_Gain_Ite] = RUN_OPT_Protocol(K_rice);
```

3. Generate performance plots:
```matlab
% Plot EH power vs transmit power
plot_EH_vs_Pt();
```

### Main Scripts

- **`RUN_OPT_Protocol.m`**: Complete three-stage protocol execution
- **`RUN_OPT_Protocol_with_Pt.m`**: Performance vs. transmit power analysis
- **`Test_Sensing.m`**: Sensing performance evaluation
- **`Test_6DMA.m`**: System performance comparison
- **`plot_EH_vs_Pt.m`**: Generate performance figures

## Three-Stage Protocol

### Stage I: Uplink Sensing
- Users transmit sensing signals sequentially (TDMA)
- Base station performs holographic sensing
- Extract angular information (θ, φ) using FFT

**Key files**: `Sensing_Algorithm.m`, `Sensing_Holographic_HalfWave.m`

### Stage II: Orientation Adjustment
- Optimize 6DMHS positions and rotations
- Design holographic beamforming
- Align maximum beamforming gain direction with users

**Key files**: `optimize_system.m`, `optimize_spatial_position.m`

### Stage III: Downlink Transmission
- Estimate equivalent CSI
- Optimize digital beamforming and power splitting
- Perform IDET transmission

**Key files**: `calculate_EH_power.m`, `generate_effective_channel.m`

## System Parameters

Default parameters (consistent with paper):

| Parameter | Value | Description |
|-----------|-------|-------------|
| fc | 30 GHz | Carrier frequency |
| M | 32×32 | RHS element array size |
| K | 4 | Number of users |
| B | 4 | Number of RHS |
| Q | 1 | Number of feeds per RHS |
| K_rice | 10 | Rician factor |
| Pt | 40 dBm | Transmit power |
| d | λ/4 | Element spacing |

## Documentation

For detailed analysis of each code file and its relationship to the paper:

- **English**: See [CODE_ANALYSIS.md](CODE_ANALYSIS.md)
- **中文**: See [代码分析文档.md](代码分析文档.md)

These documents provide:
- Function and purpose of each MATLAB file
- Mapping between code and paper sections
- Correspondence to equations and algorithms
- Execution flow and recommended reading order

## Performance Benchmarks

The implementation supports comparison with several baseline schemes:

1. **FPA (Fixed Position Antenna)**: Traditional fixed RHS configuration
2. **6DMHS with rotation only**: Only rotation optimization
3. **6DMHS with translation only**: Only position optimization
4. **Full 6DMHS**: Complete position and rotation optimization
5. **LS-based sensing**: Traditional least-squares sensing method

## Reproducing Paper Results

To reproduce figures from the paper:

```matlab
% Figure 5: EH power vs transmit power
RUN_OPT_Protocol_with_Pt();
plot_EH_vs_Pt_final();

% Figure 7-8: Sensing performance
Test_Sensing();

% Figure 9: Feed position optimization
% (Configure parameters in RUN_OPT_Protocol.m)
```

## Code Structure Overview

### Modules

1. **Sensing Module**
   - Holographic sensing implementation
   - FFT-based angle detection
   - Half-wavelength sensing validation

2. **Optimization Module**
   - Joint position-rotation-beamforming optimization
   - Rodrigues rotation matrix calculation
   - Discrete position generation

3. **Transmission Module**
   - Channel generation (Rician model)
   - EH power calculation
   - Equivalent channel generation

4. **Evaluation Module**
   - Performance comparison
   - Result plotting
   - Benchmark testing

## Citation

If you use this code in your research, please cite:

```bibtex
@article{wang2024sixd,
  title={6D Movable Holographic Surface Assisted Integrated Data and Energy Transfer: A Sensing Enhanced Approach},
  author={Wang, Zhonglun and Zhao, Yizhe and Hu, Gangming and Zheng, Yali and Yang, Kun},
  journal={To be published},
  note={Preprint available},
  year={2024}
}
```

*Note: Journal name will be updated once the paper is published.*

## Key Results

The proposed 6DMHS-IDET system achieves:

- **Superior sensing accuracy**: Holographic sensing outperforms LS-based methods
- **Higher EH power**: 6DMHS significantly improves over fixed RHS
- **Reduced complexity**: Sensing-enhanced protocol reduces computational burden
- **Better IDET performance**: Outperforms perfect CSI schemes by reducing pilot overhead

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Authors

- **Zhonglun Wang** - University of Electronic Science and Technology of China
- **Yizhe Zhao** - University of Electronic Science and Technology of China
- **Gangming Hu** - University of Electronic Science and Technology of China
- **Yali Zheng** - University of Essex
- **Kun Yang** - University of Essex

## License

This project is for academic research purposes only. The code is provided as-is for reproducing the results from the research paper. For commercial use or redistribution, please contact the authors. Please refer to the paper for detailed methodology and results.

## Contact

For questions or issues, please:
- Open an issue on GitHub
- Contact the authors via email (see paper)

## Acknowledgments

Please refer to the research paper for funding information and acknowledgments.

---

**Keywords**: 6D Movable Antenna (6DMA), Reconfigurable Holographic Surface (RHS), Integrated Data and Energy Transfer (IDET), Holographic Sensing, Beamforming Optimization
