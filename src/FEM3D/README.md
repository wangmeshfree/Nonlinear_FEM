# FEM3D — MAE 5036 三维线弹性有限元

在 FEM2D 教学结构基础上建立的独立三维版本。采用八节点三线性六面体 H8，每节点三个平移自由度。求解小变形、各向同性线弹性静力问题。二维圆孔板和悬臂梁解析解没有被当作三维实体精确解复用。

## 四个主要文件

| 文件 | 作用 |
| --- | --- |
| `main.py` | 选择算例、组装、求解、保存位移/应力/反力、生成三维变形图和收敛图。 |
| `fem3d_input_decks.py` | `problem_settings`：材料、网格、体力、边界标志、解析解及收敛开关。 |
| `fem3d_exact_solutions.py` | 均匀拉伸及三维制造解的位移、应力、体力、牵引；三维本构矩阵。 |
| `fem3d_lib.py` | H8 形函数、Jacobian、B 矩阵、体/面积分、整体组装、约束、应力、误差和 VTK 输出。 |

另有 `test_fem3d.py` 验证程序及 `Inputs/generate_meshes.py` 网格生成工具。现成 `.k` 网格已经提供。

## 安装与运行

```sh
python -m pip install numpy matplotlib scipy
python main.py
python test_fem3d.py
```

在 `main.py` 顶部修改 `active_problem` 为 `Manufactured_cube`、`Uniform_tension` 或 `Cantilever_solid`。默认运行 `Manufactured_cube`。在对应输入配置中设置 `convergence=False` 可只算主网格。入口按教学步骤从上到下直接执行，没有命令行解析或 `main()` 包装。图像直接保存，不弹出阻塞窗口。

## 三个示例

- **Uniform_tension**：单位立方体的 x 向均匀拉伸，左面施加解析位移、右面解析牵引，其余面无牵引。精确位移为 `(sigma*x/E, -nu*sigma*y/E, -nu*sigma*z/E)`。用于 affine patch test。
- **Manufactured_cube**：单位立方体，精确位移 `u=(A*sin(pi*x)*sin(pi*y)*sin(pi*z),0,0)`，六面施加解析位移，体力由 `b=-div(sigma)` 推导。默认 2、4、6 个单元/方向的网格，检验位移 L2 和能量收敛。
- **Cantilever_solid**：`10 x 2 x 1` 实体梁，左端完全固定，右面均匀 y 向牵引 `0.01`，合力 `0.02`。无精确解，不输出解析误差；默认载荷保持位移较小。不是节点集中力。

三维没有 `plane_stress`/`plane_strain` 参数；本构是完整三维 Hooke 定律。

## 标志及配置

`dirichlet`、`traction` 的标志均是 `(x,y,z)` 三个分量：

- `0`：不约束位移，或该集合不贡献该方向的牵引。
- `1`：使用 `dirichlet_values` 或 `traction_values` 对应分量，缺省为零。
- `99`：使用 `solution` 对象的解析位移/牵引。

```python
'dirichlet': {'XMIN': (1, 1, 1)},
'dirichlet_values': {'XMIN': (0., 0., 0.)},
 'traction': {'XMAX': (0, 1, 0)},
 'traction_values': {'XMAX': (0., .01, 0.)},
```

牵引是全局坐标的力/面积，程序通过面积分转换成节点力；体力是力/体积。没有 `solution` 时由常量 `body_force` 给出，有 `solution` 时使用解析类方法。模型全部采用一致单位，三维不再乘二维单位厚度。

`material_properties=[E,nu]`，要求 `E>0,-1<nu<0.5`；`integration_order=2` 为 `2x2x2` 体积分和 `2x2` 面积分。一点积分不支持，以免引入未处理的 hourglass 模式。误差积分默认每方向三点，并逐点使用正确权重。

## H8 顺序与网格

自然坐标局部节点顺序：

```text
1 (-1,-1,-1)   2 (+1,-1,-1)   3 (+1,+1,-1)   4 (-1,+1,-1)
5 (-1,-1,+1)   6 (+1,-1,+1)   7 (+1,+1,+1)   8 (-1,+1,+1)
```

支持逗号分隔 LS-DYNA 教学网格：`*NODE`、单行八节点 `*ELEMENT_SOLID`、四节点边界 `*ELEMENT_SHELL`、命名 `*SET_NODE_LIST` 和 `*SET_SHELL_LIST`。此 reader 不是完整 LS-DYNA 解析器，不支持固定宽度、两行连接表、tetrahedral 或退化实体单元。

集合名称由对应 set keyword 后的 `$ NAME` 注释给出。同名节点与面集合合并。提供 `XMIN/XMAX/YMIN/YMAX/ZMIN/ZMAX` 六面集合。边界 shell 仅作面力积分，不增加壳刚度；其节点顺序应使两条参数切向的叉积指向外法向。程序检查体单元 Jacobian 必须为正，重复位移约束合并，冲突约束报错。

重新生成输入：

```sh
python Inputs/generate_meshes.py
```

修改几何或密度后需要同步修改配置与网格；程序不会因为改变解析参数而自动重建网格。

## 输出

每个问题保存到 `Results_3d_<problem>_<mesh>/`：

- `nodal_displacements.txt`：`x y z ux uy uz`。
- `nodal_stresses.txt`：角点应力邻接平均，顺序 `xx yy zz yz xz xy`。
- `gauss_stresses.txt`：高斯点坐标与应力，同一分量顺序。
- `reactions.txt`：各节点三个反力分量；自由节点仅有数值残差。
- `deformation.png`：原网格和放大变形网格，标明放大倍数。
- `solution.vtu`：ParaView 可打开；包含三维位移与六分量节点应力，可使用 Warp By Vector。
- `summary.json`：残差、载荷和反力合力、误差与收敛阶。
- `convergence.txt`、`convergence.png`：有精确解且启用收敛时生成。

位移相对 L2 范数为 `sqrt(integral |u_h-u|^2 / integral |u|^2)`；能量误差为 `sqrt(integral (eps_h-eps)^T D (eps_h-eps) / integral eps^T D eps)`。应变顺序为 `xx,yy,zz,gamma_yz,gamma_xz,gamma_xy`，剪应变使用工程剪应变，剪应力不乘 2。精确解范数为零时相对误差返回 NaN。节点应力平均属于后处理，不代替高斯点精度检验。

## 验证与适用范围

`test_fem3d.py` 检查形函数完备性、畸变单元线性位移再现、刚度对称性及六个刚体模态、三点积分体力合力、均匀拉伸及反力平衡、制造解体力与应力散度一致性、网格误差下降、实体悬臂梁合力平衡。

这是教学版：使用 SciPy COO 组装、CSR 存储及稀疏直接求解。稀疏分解仍有额外内存和计算成本。细长梁粗 H8 网格可能剪切锁定，接近不可压缩材料可能体积锁定。没有实现塑性、几何非线性、接触、动力学或稳定化。新增三维孔洞问题需要真实三维网格和对应边界/解析基准，不能直接把二维 `.k` 壳网格交给本程序。

## 稀疏刚度矩阵与内存

整体 K 使用数值 COO 三元组数组收集贡献，再转 CSR，自动合并重复位置。F、U 保留 NumPy 数组。约束采用自由/固定 DOF 分块，求解使用 `scipy.sparse.linalg.spsolve`，不转换成密集矩阵。

`3d_cube_16.k` 实測：4,913 节点、4,096 H8、14,739 DOF。旧 float64 密集 K 的存储需求按 N*N*8 精确计算为 1.619 GiB；新 CSR 的 data/indices/indptr 实测合计 12.05 MiB，约减少 138 倍。COO 组装三元组数组为 36 MiB。测试进程峰值 RSS 约 338 MiB，包含 Python、库、子矩阵及直接分解，不是单独 K 的大小。为避免大量内存分配，16 层测试没有实际建立密集 K；仅在小网格验证密集与稀疏结果一致。

在自己的 Python 环境中安装 `scipy`。测试环境安装不等于所有 IDE/Jupyter 环境已安装。
