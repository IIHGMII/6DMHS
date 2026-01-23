# 6DMHS（6D Movable Holographic Surface）论文与仿真代码

本仓库包含：

- 论文 TeX：`RHS_WET_NPA.tex`
- MATLAB 仿真代码：`Proposed-v2/`

## 推荐入口

- 复现/跑通三阶段流程并绘制“EH vs Pt”：`Proposed-v2/fig5.m`
- 逐点计算（函数形式）：`Proposed-v2/RUN_OPT_Protocol_with_Pt.m`

## 依赖

- MATLAB
- `Sensing_Algorithm.m` 中使用 `normrnd`，通常需要 Statistics and Machine Learning Toolbox
- `fun_run_R00.m` 等优化脚本依赖 CVX（若安装了 MOSEK 会优先使用，否则使用 CVX 默认求解器）

## 文档

- 论文与代码对应关系：`docs/paper-code-map.md`

## 运行提示

- 若未安装 CVX，`fun_run_R00.m` 会报错；请先按 CVX 官方流程安装并 `cvx_setup`。
- `RUN_OPT_Protocol_with_Pt.m` 默认使用较小阵列（`Mx=My=6`）以便快速跑通；如需论文同规模可改为 `32`。
