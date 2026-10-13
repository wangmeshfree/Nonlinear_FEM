# FEM2D 使用说明 — MAE 5036

二维小变形、线弹性静力分析的教学程序。使用四节点双线性四边形（Q4）单元，每个节点有两个自由度 `(u_x, u_y)`，支持平面应力和平面应变，按单位厚度计算。

## 1. 四个 Python 文件的作用

| 文件 | 作用 | 同学们通常需要修改什么 |
| --- | --- | --- |
| `main.py` | 程序入口。选择算例，建立材料矩阵，读取网格，组装刚度与载荷，施加边界条件，求解 `K U = F`，如果问题有解析解，则做收敛分析。 | 修改 `active_problem` 选择要运行的算例；理解这里的求解流程。 |
| `fem2d_problem_lib.py` | 是求解问题库，里面的 `problem_settings`，集中定义读取哪些网格、材料、高斯积分点数、边界集合、载荷。 | 可以修改现有算例参数，或添加自己的算例，通常被称为Input deck。 |
| `fem2d_exact_solutions.py` | 收集了弹性力学问题的解析解，已有 `PlateWithHole` 和 `CantileverBeam` 两个问题。每个类提供 `displacement(x,y)`、`stress(x,y)`、`body_force(x,y)` 和 `traction(x,y,normal)`。 |
| `fem2d_lib.py` | 是所有2d有限元程序共用的函数库：网格读取、Q4 形函数、Jacobian、单元积分、整体组装、位移和面力边界处理、应力计算、绘图、误差与收敛分析。 | 学习和实现有限元算法；仅改变算例参数时通常不用修改。 |

文件之间的关系：`main.py` 读取 `fem2d_problem_lib.py` 的配置并调用 `fem2d_lib.py`进行计算；如果问题有解析解，配置中的 `solution` 指向 `fem2d_exact_solutions.py`。`Inputs/` 用于存放已经生成好的网格文件，我们使用Gmsh去生成.k文件让程序进行读取。

## 2. 安装与运行

运行脚本需要 Python、NumPy 和 Matplotlib：

```sh
python -m pip install numpy matplotlib
python main.py
```

在本目录运行上述命令，或在 IDE 中运行 `main.py`。已有 `.k` 网格可以直接使用；只有重新生成网格时才需要 Gmsh。

## 3. 选择算例

在 `main.py` 中设置一个名称：

```python
active_problem = 'Cantilever_beam'
# active_problem = 'Plate_with_hole'
# active_problem = 'Beam_with_tip_load'
```

名称必须是 `fem2d_problem_lib.py` 中 `problem_settings` 的键。当前三个算例的区别如下；以代码中的最新配置为准。

| 算例 | 当前材料假设 | 边界和载荷 | 解析比较 |
| --- | --- | --- | --- |
| `Plate_with_hole` | 平面应力 | 四分之一圆孔板；左边 `u_x=0`、底边 `u_y=0`；上边和右边施加 Kirsch 解析牵引。 | 有，默认进行收敛分析。 |
| `Cantilever_beam` | 平面应变 | 左边施加解析位移；右端施加解析应力对应的抛物线分布剪力。 | 有，默认进行收敛分析。 |
| `Beam_with_tip_load` | 平面应力 | 左端两个方向固定；右边施加均匀竖向面力。 | 无，不进行解析误差分析。 |

### 两个解析基准的含义

带孔板的参数为 `[孔半径 a, 远场 x 向拉应力 sigma0]`。Kirsch 解属于无限板；有限计算区域的外边界使用解析牵引 `sigma @ normal`，而不是简单地将上边设为自由、右边设为均匀拉力。孔边不施加载荷，表示自由边界。

悬臂梁的参数为 `[长度 L, 高度 H, 右端剪力合力 P]`，坐标范围为 `0 <= x <= L`、`-H/2 <= y <= H/2`。`P > 0` 表示向上加载。解析位移在左端一般不为零，因此 `(99,99)` 与左边整条边完全固定的 `(1,1)` 不是同一个边界条件。

`Beam_with_tip_load` 虽然名称含有 tip load，当前实现施加的是右端均匀分布面力，不是单个节点的集中力。默认高度为 2，竖向面力为 2.5，单位厚度下竖向合力为 5。

## 4. 配置参数

主要修改 `fem2d_problem_lib.py`：

| 配置项 | 含义 |
| --- | --- |
| `solution` | 可选的解析解类。存在时，主程序用 `solution(*model_characteristics, E, nu, plane_problem_type)` 创建实例。 |
| `model_characteristics` | 解析类所需的几何和载荷参数。几何参数必须与网格一致；无解析解时该字段不会自动生成载荷。 |
| `mesh_file` | 主计算网格，路径相对于 `Inputs/`。 |
| `mesh_files` | 收敛分析使用的网格列表，应从粗到细，且几何、材料和边界问题一致。 |
| `material_properties` | `[E, nu]`。 |
| `integration_order` | 每个自然坐标方向的高斯点数；`2` 表示 Q4 使用 `2 x 2` 积分。也用于边界积分及当前误差积分，不改变网格密度。 |
| `2d_problem_type` | `'plane_stress'` 或 `'plane_strain'`。 |
| `dirichlet` / `dirichlet_values` | 节点集合上的位移标志和指定值。 |
| `traction` / `traction_values` | 边界线单元集合上的面力标志和指定值。 |
| `body_force` | 无 `solution` 时使用的常量 `(bx, by)`，默认 `(0,0)`；有解析解时调用解析类的 `body_force`。 |
| `convergence` | `True` 或 `False`；当前收敛分析还要求存在 `solution`。 |
| `title` | 网格图标题。 |

### 边界标志：0、1、99

每组两个标志分别对应全局 x、y 分量。

| 标志 | `dirichlet` 中的含义 | `traction` 中的含义 |
| --- | --- | --- |
| `0` | 不规定该位移分量。 | 该集合不贡献该方向的面力。 |
| `1` | 取 `dirichlet_values` 对应值，默认零。 | 取 `traction_values` 对应值，默认零。 |
| `99` | 取解析解的位移分量。 | 取解析牵引的分量。 |

`1` 不是单位位移或单位力；只有未提供对应值时，它才意味着指定零值。`99` 必须有 `solution`。

例如，左端固定、右端施加均匀竖向面力：

```python
'dirichlet': {'FIXED_LEFT': (1, 1)},
'dirichlet_values': {'FIXED_LEFT': (0.0, 0.0)},
'traction': {'TRACTION_RIGHT': (0, 1)},
'traction_values': {'TRACTION_RIGHT': (0.0, 2.5)},
```

面力使用全局坐标，按单位长度、单位厚度输入，程序通过边界积分换算为节点力。集合名称必须存在于 `.k` 网格。当前版本不会自动增加中心点位移约束，因此自定义问题需要自行提供足以去除刚体运动的边界条件。

## 5. 输出结果与收敛分析

程序在本目录建立 `Results_2d_<问题名>_<主网格名>/`，主计算输出：

| 文件 | 内容 |
| --- | --- |
| `nodal_displacements.txt` | `x, y, u_x, u_y`。 |
| `nodal_stresses.txt` | `x, y, sigma_11, sigma_22, sigma_12`。 |
| `mesh.png` | 单元网格。 |
| `displacement.png` | 位移场；有解析解时同时显示解析比较。 |
| `stress.png` | 节点应力场；有解析解时同时显示解析比较。 |
| `convergence.png` | 启用收敛分析时生成的相对位移 L2 与相对能量误差图。 |

节点应力是各单元角点应力的平均值，属于后处理结果。Q4 位移较准确并不意味着粗网格应力也同样准确。高斯点应力计算与绘图函数仍在库中，但当前 `main.py` 没有调用它们，也不会自动生成高斯点应力文件。

启用收敛分析：

```python
'convergence': True,
```

还需配置 `solution` 和 `mesh_files`。程序重新求解每个网格，控制台打印相对位移 L2 误差、相对能量误差，以及平均收敛阶。两种误差均通过域积分计算，不是节点误差平均。典型光滑解下，Q4 位移 L2 误差可接近二阶、能量误差接近一阶；实际结果取决于网格和解的性质。

当前版本只打印收敛数据并保存图，不自动生成 `convergence.txt`。改变参数后再次运行会覆盖同名结果目录中的文件，需要比较时请先保存旧结果。

## 6. 添加自己的算例

1. 准备 Q4 的 LS-DYNA `.k` 网格，包括节点、四节点单元，以及载荷边界的线单元与命名集合。
2. 在 `problem_settings` 中新增条目，配置网格、材料和边界条件，然后在 `main.py` 中选择该名称。
3. 没有解析解时，省略 `solution`，用标志 `0/1`、指定值和 `body_force` 定义问题；设置 `'convergence': False`。无需仅为新问题名称修改 FEM 求解器。
4. 若要解析误差比较，在 `fem2d_exact_solutions.py` 中增加相同接口的解析解类，导入配置文件并赋给 `solution`。
5. 检查边界条件是否与解析解一致、集合名称是否正确，以及单元和边界线的节点顺序。解析牵引的外法向由边界线方向计算，反向连接会改变载荷方向。



建议阅读顺序：先运行现有算例，再看 `main.py` 的流程和配置，最后进入单元矩阵、边界积分及组装函数。
