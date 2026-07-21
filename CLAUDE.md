# 20251222_Zesheng_Xuanwu_brain

脑组织代谢通量分析（MFA）项目，用于阿尔茨海默病（AD）与司美格鲁肽（Semaglutide）相关的代谢研究。

## 项目背景

本项目研究阿尔茨海默病模型小鼠在司美格鲁肽干预后的脑组织代谢变化，通过稳定同位素示踪和代谢通量分析技术，揭示 AD 相关的代谢重编程机制。

## 项目结构

```
20251222_Zesheng_Xuanwu_brain/
├── main.py                 # 主入口脚本
├── scripts/
│   ├── config.py           # 脚本配置
│   ├── main.py             # 脚本入口
│   ├── src/                # 核心源代码
│   │   ├── core/           # MFA核心算法（求解器、模型、采样器）
│   │   ├── common/         # 通用配置和工具
│   │   ├── common_parallel_solver/  # 并行求解器
│   │   ├── mfa_analysis/   # MFA分析主流程
│   │   ├── io_functions/   # 输入输出函数
│   │   └── model_data_combination/  # 数据与模型组合配置
│   ├── figures/            # 图表生成
│   │   ├── common/         # 图表通用配置
│   │   ├── figure_content/ # 图表内容定义
│   │   ├── figure_elements/  # 图表元素
│   │   └── figure_plotting_package/  # 图表绑制工具包
│   ├── data/               # 数据文件
│   │   ├── data_20251115_ad_sema/  # AD Semaglutide 实验数据
│   │   └── config.py       # 数据配置
│   ├── model/              # 代谢网络模型
│   │   ├── complete_model/ # 完整模型
│   │   └── multi_tissue_compartmentalized_network/  # 多组织区室化网络
│   └── output/             # 输出结果
├── docker_commands.sh      # Docker 命令
└── python_main_execute.sh  # Python 执行脚本
```

## 核心模块说明

### scripts/src/core/
**MFA 核心算法库**，包含代谢通量分析的核心实现：
- `solver/` - 优化求解器（SLSQP、Numba 加速）
- `model/` - EMU 模型和代谢网络模型
- `sampler/` - 采样器（用于初始值生成）
- `common/` - 通用类和配置
- `data/` - 数据结构定义

### scripts/src/common/
**通用工具模块**，包含项目级别的配置和辅助函数。

### scripts/src/common_parallel_solver/
**并行求解器**，实现多核并行优化：
- `common_parallel_solver.py` - 并行求解主类
- `common_parallel_predictor.py` - 并行 MID 预测
- `feasible_solution_generator.py` - 可行解生成器
- `result_class.py` - 结果数据结构

### scripts/src/mfa_analysis/
**MFA 分析主流程**，包含完整的分析管道：
- `mfa_main_functions.py` - 主分析函数
- `argument_and_entrance.py` - 命令行参数解析
- `mfa_config_generator.py` - 配置生成
- `result_processing_and_plotting.py` - 结果处理和绑图

### scripts/src/io_functions/
**输入输出函数**，处理数据读取和结果输出：
- `result_processing_functions.py` - 结果处理
- `result_output_functions.py` - 结果输出
- `plotting_functions.py` - 基础绑图函数

### scripts/src/model_data_combination/
**数据与模型组合配置**，定义具体实验的数据和参数：
- `data_20251115_ad_sema.py` - AD Semaglutide 实验配置
- `common_data_model_loader.py` - 数据模型加载器

### scripts/figures/common/
**图表通用模块**，包含图表生成的通用配置和素材。

### scripts/figures/figure_content/
**图表内容定义**，定义具体图表的内容结构。

### scripts/figures/figure_elements/
**图表元素**，各种具体图表类型的实现：
- `data_figure_plotting.py` - 数据图表
- `mid_comparison_figure.py` - MID 对比图
- `flux_comparison_figure.py` - 通量对比图
- `violin_box_data_figure.py` - 小提琴/箱线图

### scripts/figures/figure_plotting_package/
**图表绑制工具包**，面向科研论文的矢量图表生成库（详见该目录下的 CLAUDE.md）。

### scripts/data/
**数据目录**，包含实验数据和处理配置：
- `data_20251115_ad_sema/` - AD Semaglutide 实验的原始数据
- `config.py` - 数据路径和加载配置
- `data_loader.py` - 数据加载函数

### scripts/model/
**代谢网络模型目录**，包含代谢网络定义：
- `complete_model/` - 完整代谢网络模型
- `multi_tissue_compartmentalized_network/` - 多组织区室化网络模型
- `model_loader.py` - 模型加载函数

---

## ⚠️ 重要：核心模块保护

以下目录是项目的**核心代码库**，具有复杂的依赖关系：

```
scripts/src/core/
scripts/src/common/
scripts/src/common_parallel_solver/
scripts/src/mfa_analysis/
scripts/src/io_functions/
scripts/figures/figure_elements/
scripts/figures/figure_plotting_package/
```

**请勿修改这些目录**，除非你明确知道自己在做什么。这些模块之间存在复杂的依赖关系，修改可能导致：
- 破坏已有的测试用例
- 引入难以追踪的 bug
- 影响分析结果的正确性

如需修改，请：
1. 先与项目负责人确认
2. 确保修改向后兼容
3. 进行完整测试

## 代码规范

### 使用 ParameterName 常量

构建 `figure_data_parameter_dict` 或 `figure_config_dict` 时，**必须使用 `ParameterName` 类中的字符串常量作为 key**，禁止使用裸字符串。

```python
# 正确
from ...common.config import ParameterName
params = {
    ParameterName.horizontal: True,
    ParameterName.x_tick_labels_list: labels,
    ParameterName.figure_config_dict: {
        ParameterName.column_width: 0.15,
        ParameterName.edge: 0.2,
        ParameterName.x_label_format_dict: {'font_size': 12},
        ParameterName.supplementary_text_format_dict: {
            ParameterName.font_size: 10,
            ParameterName.width: 0.1,
        },
    },
}

# 错误 — 裸字符串
params = {
    'horizontal': True,
    'x_tick_labels_list': labels,
    'figure_config_dict': {'column_width': 0.15},
}
```

例外：`metabolomics_style` / `common_figure_materials` 中的非框架函数（如 `recs`、`all_pct`、`x_label`、`y_label` 等作为 `default_parameter_extract` 的自定义 key）可使用裸字符串，因为它们是用户层参数，不属于框架 `ParameterName`。

### Subfigure 统一模式

所有 Subfigure 使用 `element = Elements.XxxFigure(**{ParameterName.figure_data_parameter_dict: params, ParameterName.scale: 0.XX, ...})` + `move_center(element, 0.XX, subfigure_size)` 模式。scale 值手动硬编码，不使用 probe 动态计算。

### 文件换行规范

所有代码文件（`.py`、`.md`、`.json` 等）使用 Unix 风格换行（LF），禁止 CRLF。文件末尾必须以 `\n` 结尾。

## 使用方法

```bash
# 设置 Python 路径
export PYTHONPATH=$PYTHONPATH:`pwd`

# 实验数据绑图
python main.py computation experimental_mfa raw_experimental_data_plotting

# 通量分析（使用 20 核并行）
python main.py computation experimental_mfa flux_analysis -p 20

# 结果处理
python main.py computation experimental_mfa result_process
```

## Docker 支持

项目支持 Docker 运行，详见 `docker_commands.sh`：

```bash
# 启动容器
bash docker_commands.sh -p 20 start

# 查看容器日志
bash docker_commands.sh log

# 停止容器
bash docker_commands.sh stop

# 处理结果
bash docker_commands.sh result

# 构建镜像
bash docker_commands.sh build
```

## 依赖

- Python 3.12+
- NumPy
- Matplotlib
- SciPy（SLSQP 求解器）
- Numba（可选加速）