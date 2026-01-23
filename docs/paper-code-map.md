# 论文与 `Proposed-v2/` 代码对应关系（基于 `RHS_WET_NPA.tex`）

## 总体对应：三阶段协议（Stage I/II/III）

论文在协议设计部分给出三阶段流程，代码中可按下列链路理解：

- Stage I：上行感知（Uplink Sensing）
  - 初始姿态/几何：`Proposed-v2/Orientation_Initial.m`、`Proposed-v2/Orientation_uniformSpherePoints.m`
  - 信道生成（含感知相关子信道输出）：`Proposed-v2/Channel_Generation_init.m`
  - 全息感知 + FFT 检测：`Proposed-v2/Sensing_Algorithm.m`
- Stage II：姿态调整（Orientation Adjustment）
  - 旋转对准：示例在 `Proposed-v2/fig5.m` 中用罗德里格斯旋转更新 `Rot_Matrix{k}`
  - 平移离散选择 + 全息/数字波束（示例实现）：`Proposed-v2/optimize_system.m`
- Stage III：下行等效信道 + 数能优化（Downlink Transmission）
  - 等效信道生成：`Proposed-v2/generate_effective_channel.m`
  - 数字波束 + 功率分割（FP/CVX 交替优化）：`Proposed-v2/fun_run_R00.m`

## 按论文章节映射

- System Model
  - 坐标/姿态与离散位置采样：`Proposed-v2/Orientation_Initial.m`
  - Rician 信道与可见性判定（法向量与入射方向点积）：`Proposed-v2/Channel_Generation_init.m`
- Directional Beamforming Gain & Protocol
  - 端到端串联示例：`Proposed-v2/fig5.m`、`Proposed-v2/main1.m`
- Holographic-based Sensing & Orientation Adjustment
  - 全息图构造 + 2D FFT 峰值检测（对应论文 Theorem 1 的 FFT-based detection 叙述）：`Proposed-v2/Sensing_Algorithm.m`
  - 旋转矩阵更新（最大增益方向对准估计方向）：`Proposed-v2/fig5.m`
- Joint Beamforming and Power Splitter Design
  - 最小 EH 功率最大化的 CVX/FP 交替更新：`Proposed-v2/fun_run_R00.m`

## 关键变量对照（论文符号 → 代码变量）

- RHS 数量/用户数：论文常设 `K=B` → 代码 `K`、`B`
- 元素规模：`M=Mx*My` → 代码 `Mx`、`My`、`M`
- 旋转矩阵：论文 `R_b` → 代码 `Rot_Matrix{b}`
- 平移离散选择：论文 `s_b`（one-hot 选择向量）→ 代码 `s(:,b)`
- 等效信道：论文 `\\bar{h}` → 代码 `H_eff`
- 功率分割：论文 `\\rho_k` → 代码 `rho(k)`

## 推荐复现实验

- “EH vs Pt”曲线：直接运行 `Proposed-v2/fig5.m`
- 逐点计算：调用 `Proposed-v2/RUN_OPT_Protocol_with_Pt.m`（注意需要 CVX）

## 备注（易踩坑）

- `Proposed-v2/plot_EH_vs_Pt_final.m` 里的 “LoS-Only / Perfect CSI” 当前是按比例缩放的占位写法，不是严格实现。
- `Proposed-v2/plot_EH_vs_Pt_global.m` 引用了 `RUN_OPT_Protocol_modified`（仓库内未提供），建议优先用 `Proposed-v2/fig5.m`。
- 仓库里部分旧的 `test*.m`/`Test_*.m` 脚本可能仍使用旧接口（参数数量不匹配）；若报错，请以 `Proposed-v2/fig5.m` 的调用链为准。
