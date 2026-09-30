# FEM2D `problem_settings` 参考示例

本文只提供配置模板，不会被 `main.py` 导入或自动执行。使用时，在 `main.py` 的 `problem_settings` 中修改对应条目，然后设置 `active_problem = 'Plate_with_hole'` 或 `'Cantilever_beam'`。当前求解器只识别这两个问题名称；不要仅靠添加第三个名称来创建新问题。

## 每个设置控制什么

| 设置 | 含义 | 修改时要注意 |
| --- | --- | --- |
| `model_characteristics` | 带孔板为 `[孔半径 R, 远场 x 向拉力 Tx]`；悬臂梁为 `[长度 L, 高度 H, 端部竖向合力 P]` | `R`、`L`、`H` 必须与网格几何一致。`Tx`、`P` 是载荷幅值，可改变大小或符号。 |
| `mesh_file` | 主计算使用的 `.k` 网格 | 路径相对于 `src/FEM2D/Inputs/`。 |
| `mesh_files` | 收敛分析依次使用的网格 | 通常从粗到细；仅 `run_convergence_study = True` 时使用。 |
| `material_properties` | `[E, nu]`，弹性模量和泊松比 | 改变材料后，刚度和位移会改变；应力边界载荷仍由解析应力场计算。 |
| `integration_order` | Q4 域积分和边界积分的一维高斯点数 | 当前示例使用 `2`，修改它不会改变网格密度。 |
| `dirichlet` | 按网格节点集施加位移约束 | 每个方向的标志：`0` 不约束，`1` 规定为 0，`99` 规定为该节点的解析位移。 |
| `traction` | 按网格边界线单元集施加面力 | 名称必须存在于 `.k` 文件；面力由当前问题的解析应力乘外法向量计算。 |
| `2d_problem_type` | `plane_stress` 或 `plane_strain` | 改变本构矩阵；悬臂梁使用 `99` 时也改变解析位移。 |
| `title` | 网格图标题 | 不影响数值计算。 |

`mesh_file` 与 `mesh_files` 可以选不同密度，但必须属于同一几何问题。输出目录名称只包含问题名和主网格文件名；仅改变载荷、材料或边界条件后再次运行，会覆盖同名目录中的结果，请先保存需要比较的旧结果。

## 示例 A：带孔板，基准拉伸

这是当前的基准设置。将下列条目放入 `problem_settings`，选择 `active_problem = 'Plate_with_hole'`：

```python
'Plate_with_hole': {
    'model_characteristics': [1.0, 1.0],  # [R, Tx]
    'mesh_file': '2d_plate_with_hole/model_plate_with_hole_level_2.k',
    'mesh_files': [
        '2d_plate_with_hole/model_plate_with_hole_level_1.k',
        '2d_plate_with_hole/model_plate_with_hole_level_2.k',
        '2d_plate_with_hole/model_plate_with_hole_level_3.k',
    ],
    'material_properties': [1.0e3, 0.30],
    'integration_order': 2,
    'dirichlet': {'FIXED_LEFT': (1, 0), 'FIXED_BOTTOM': (0, 1)},
    'traction': ['TRACTION_RIGHT', 'TRACTION_TOP'],
    '2d_problem_type': 'plane_strain',
    'title': 'Plate with Hole',
},
```

`FIXED_LEFT` 限制 (u_x=0)，`FIXED_BOTTOM` 限制 (u_y=0)，代表四分之一模型的对称边界；`TRACTION_RIGHT` 与 `TRACTION_TOP` 是外边界的解析面力。计算结果应显示孔周围的非均匀应力；三个网格可用于观察能量范数误差变化。

### 示例 A1：拉力加倍

以上配置保持不变，仅改为：

```python
'model_characteristics': [1.0, 2.0],  # R 不变，Tx: 1 → 2
```

在线弹性范围内，位移、应力以及未归一化的能量误差约加倍；网格形状和单元数不变。`Tx` 是由解析应力施加在外边界的拉力参数，并不是在某个节点直接加力。

### 示例 A2：比较平面应力

在示例 A 的基础上仅改为：

```python
'2d_problem_type': 'plane_stress',
```

相同网格和边界面力下，本构矩阵变为平面应力矩阵，位移通常与平面应变结果不同。当前能量误差仍以解析应力场为参照。此设置适合比较两种二维材料假设，不表示改变了板厚；程序按单位厚度计算。

## 示例 B：悬臂梁，解析左端位移与向上剪力

这是与当前 Notebook 一致的基准设置。选择 `active_problem = 'Cantilever_beam'`：

```python
'Cantilever_beam': {
    'model_characteristics': [10.0, 2.0, 5.0],  # [L, H, P]
    'mesh_file': '2d_cantilever_beam/model_cantilever_level_2.k',
    'mesh_files': [
        '2d_cantilever_beam/model_cantilever_level_1.k',
        '2d_cantilever_beam/model_cantilever_level_2.k',
        '2d_cantilever_beam/model_cantilever_level_3.k',
    ],
    'material_properties': [1.0e3, 0.30],
    'integration_order': 2,
    'dirichlet': {'FIXED_LEFT': (99, 99)},
    'traction': ['TRACTION_RIGHT'],
    '2d_problem_type': 'plane_stress',
    'title': 'Cantilever Beam',
},
```

`P > 0` 表示右端向上剪力的合力。右边界实际施加的是由解析应力生成的**分布面力**，不是单个端点上的集中力。`(99, 99)` 把二维解析解在左边界每个节点的 (u_x,u_y) 作为规定位移；它不等同于整条左边界刚性夹持为零。当前三个网格的能量误差大致随网格细化下降。

### 示例 B1：改为向下剪力

在示例 B 的基础上仅改为：

```python
'model_characteristics': [10.0, 2.0, -5.0],  # P 改为负值
```

线弹性响应的位移和应力符号反转，位移绝对值及能量误差大小理论上与 `P = 5.0` 相同。这里仍是右边界的分布剪力，不是集中节点力。

### 示例 B2：比较传统零位移支承

在示例 B 的基础上仅改为：

```python
'dirichlet': {'FIXED_LEFT': (1, 0)},
```

这表示整条左边界 (u_x=0)、(u_y) 暂不约束；程序还会在左边界中点额外施加 (u_y=0)，防止刚体平移。它**不是**整条左边界 (u_x=u_y=0) 的刚性夹持。此模型与二维解析位移边界不完全一致，因此左端附近的误差及收敛结果可能明显不同。

如果要测试整条左边界刚性夹持，可改成 `(1, 1)`；这同样是与示例 B 的解析边界不同的物理问题。当前程序的额外中心点约束仍会重复规定该点的 (u_y=0)，因此不推荐直接用 `(1, 1)` 做正式比较，除非同时调整该辅助约束逻辑。

### 示例 B3：比较平面应变

在示例 B 的基础上仅改为：

```python
'2d_problem_type': 'plane_strain',
```

本构矩阵和 `(99, 99)` 使用的解析位移都会同步切换为平面应变版本。相同载荷下，平面应变一般更刚、位移较小。网格、几何尺寸和右端载荷集合不变。

### 示例 B4：悬臂梁中 `0` 和 `1` 的边界组合

`FIXED_LEFT` 的两个数字依次对应左边界节点的 `(u_x, u_y)`：

| 设置 | 数学含义 | 结果 |
| --- | --- | --- |
| `(0, 0)` | (u_x,u_y) 都不约束 | 左端完全自由，存在刚体运动，通常会导致刚度矩阵奇异。 |
| `(1, 0)` | (u_x=0)，(u_y) 不约束 | 程序会额外固定左边界中心节点的 (u_y=0)，用于去除竖向刚体平移；这是当前的传统支承测试。 |
| `(0, 1)` | (u_x) 不约束，(u_y=0) | 通常仍有水平方向刚体运动，不能作为当前悬臂梁的完整约束。 |
| `(1, 1)` | (u_x=u_y=0) | 左边界全固定，是传统刚性夹持模型；与 `(99,99)` 的解析位移边界不是同一个问题。 |
| `(99, 99)` | (u_x,u_y) 都采用解析位移 | 当前推荐的二维解析解比较设置。 |

例如，测试传统全固定左端时，可以把悬臂梁条目改成：

```python
'dirichlet': {'FIXED_LEFT': (1, 1)},
```

此时左端的两个位移分量都被规定为零。测试只限制竖向位移时可以写：

```python
'dirichlet': {'FIXED_LEFT': (0, 1)},
```

但这个设置没有约束水平方向刚体运动，直接求解很可能出现 `Singular matrix`。因此它只适合作为理解边界标志的示例，不建议作为正式悬臂梁算例。

需要注意，`1` 表示“规定为零”，并不是“施加单位位移”；`0` 表示“该方向不施加位移条件”。如果想让某个方向使用解析位移，必须使用 `99`，不能用 `1` 代替。

## 自己设计新示例时

1. 先确定仍是带孔板还是悬臂梁，并保持 `active_problem`、`problem_settings` 的键与现有解析应力/载荷函数对应。
2. 若只研究载荷大小、材料参数、平面应力/应变或网格密度，可在对应现有条目中修改。改变孔半径、梁长或梁高时，必须先准备相匹配的新网格文件。
3. 不要只改 `traction` 集合来声称创建了新的载荷方向。当前程序的面力方向由 `stress_plate_with_hole` 或 `stress_cantilever_beam` 的应力场决定；新的载荷类型需要相应的新应力/面力定义。
4. 检查 `dirichlet` 和 `traction` 的集合名称是否存在于新 `.k` 网格中，并确保约束足以去除刚体运动。
5. 查看 `Results_2d_.../` 下的位移、应力、网格与收敛结果；不要把某一网格上的结果直接当作解析解。
