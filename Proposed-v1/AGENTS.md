# Repository Guidelines

## 项目结构与模块组织
- 根目录以 MATLAB 脚本/函数为主（`*.m`）。核心流程通常由脚本驱动，配套函数完成建模与优化。
- 典型入口脚本：`reproduce_fig5.m`（复现实验图）、`main1.m`、`RUN_OPT_Protocol.m`、`test_6DMA_run.m`、`Test_Sensing.m`、`Test_6DMA.m`。
- 关键模块示例：
  - 信道与场景：`Channel_Generation_init.m`
  - 姿态/几何：`Orientation_*.m`
  - 感知算法：`Sensing_Algorithm.m`、`Sensing_Holographic*.m`
  - 优化（FP/CVX）：`fun_run_R00.m`、`Benchmark_Sensing.m`
- 结果与资产：`Data_record/`（如 `fig5_results.mat`、`fig5_reproduce.png`）。

## 构建、测试与本地运行命令
本仓库无传统“构建”步骤，主要通过 MATLAB 运行脚本：
- MATLAB 命令行（推荐，便于 CI/复现）：
  - `matlab -batch "reproduce_fig5"`：生成 Fig.5 曲线图（依赖 CVX，建议 MOSEK）。
  - `matlab -batch "test_main"`：快速画图/ sanity check（脚本型测试）。
- MATLAB GUI：在根目录直接运行同名脚本（例如在 Command Window 输入 `reproduce_fig5`）。

## 代码风格与命名约定
- 缩进：保持与现有代码一致（建议 4 空格）；`end` 对齐；避免混用 tab。
- 命名：文件名沿用当前风格（如 `Channel_Generation_init.m`、`Orientation_Initial.m`）；计数/维度常用大写（`K,B,M`），坐标常用前缀 `Coor_`。
- 可复现性：新增仿真脚本请设置随机种子（例如 `rng(1,'twister')`），并把可调参数集中在脚本顶部。
- 目录与路径：脚本默认在仓库根目录运行；如需显式加路径，使用 `addpath(genpath(pwd))`，避免写死本机绝对路径。
- 文件模式示例：算法族用前缀聚类（`Sensing_*`、`Orientation_*`），实验复现用动词（`reproduce_*.m`、`Run_*.m`、`Test_*.m`）。

## 测试指南
- 当前以“脚本回归”为主：优先保证 `reproduce_fig5.m` 与 `test_main.m` 可直接运行、不报错、输出图/结果合理。
- 运行较慢的脚本（如含 `parfor`、大规模 CVX）建议提供降阶参数（迭代次数、候选点数）以便快速检查。
- 若你新增函数文件（`function ...`），请同时提供一个可运行的最小示例（mini demo）脚本，说明输入/输出维度与单位。

## 提交与 Pull Request 指南
- 提交信息：历史上中英文均可，倾向“简短主题行”（例如 `初始化`、`fig5重置版`、`Add ...`）。建议可选加前缀：`fig5:`、`sensing:`、`opt:`。
- PR 要求：
  - 描述改动目的、关键参数变化、如何复现（给出运行命令与预期产物路径，如 `Data_record/`）。
  - 若改动影响图形结果，请附截图或导出的 `*.png`。
  - 避免提交 MATLAB 自动保存文件（`*.asv`）和运行日志（如 `replay_pid*.log`），除非用于复现问题。
  - 若引入/移除依赖（例如 CVX solver、toolbox），请在 PR 里明确“可选/必需”、安装方式与替代方案。

## 依赖与配置提示
- 优化相关脚本使用 CVX，并在多处指定 `cvx_solver mosek`；若本机无 MOSEK，请在本地改为已安装的求解器并在 PR 中说明。
- 含 `parfor` 的脚本需要 Parallel Computing Toolbox；如不具备，请提供可串行运行的替代路径或参数开关。

## 数据与产物管理
- 默认把可复现产物放到 `Data_record/`：图用 `*.png`，中间结果用 `*.mat`，可交互图用 `*.fig`（但 PR 更推荐同时导出 `png`）。
- 不要提交临时/本地文件：`*.asv`、`replay_pid*.log`、以及个人实验的超大 `*.mat`（除非它是论文复现必需且已在 PR 中说明来源与用途）。
- 变更会影响数值结果时，请在 PR 中说明“随机性来源”（seed、`parfor` 调度）与“单位/量纲”（例如 `dBm`、`W`、`Hz`）。

## 建议工作流（贡献者/Agent）
1. 先跑通基线：`matlab -batch "test_main"`，确认环境与路径无误。
2. 修改后做最小回归：与改动相关的脚本（例如 `reproduce_fig5.m`）至少运行一次，并记录关键参数（`Pt_dBm_list`、迭代次数、候选点数等）。
3. 输出与落盘：
   - 图：优先导出为 `Data_record/*.png`（不要只保存 `*.fig`）。
   - 数值：保存为 `Data_record/*.mat`，并在 PR 描述里写清字段含义/单位。
4. 性能注意：默认参数可能较慢；提交时请保留“快速模式”参数（如 `max_iter`、`outer_iter`）便于 reviewer 复现。
5. 提交前自查：`git diff` 确认未混入 `*.asv`/日志/临时产物，且脚本在干净工作区可直接运行。

## 常见问题排查
- CVX 报错/不可行：先检查求解器可用性（`cvx_solver`）、噪声/功率量纲（`dBm -> W` 转换），再尝试降低约束强度或减少迭代次数定位问题。
- 结果不稳定：确认是否固定 `rng(...)`，以及是否在循环中意外重置随机数或覆盖全局变量（`clear`/`clc` 的位置）。
