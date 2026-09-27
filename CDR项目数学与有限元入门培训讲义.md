---
status: reviewing
date: 2026-09-27
last_updated: 2026-09-27
course: CDR 求解器团队入门培训
audience: 具备 Python 与微积分基础、尚未系统学习有限元的新成员
tags:
  - 类型/培训讲义
  - 有限元
  - FEALPy
  - CDR
---

# CDR 项目数学与有限元入门培训讲义

## 0. 培训目标

这份讲义帮助新成员建立一条完整主线：

```text
物理收支
  → 偏微分方程
  → 弱形式
  → 网格与基函数
  → 有限元矩阵
  → 多分量块系统
  → 时间推进
  → 强对流稳定化
  → 数值验证
  → FEALPy 实现
```

完成学习后，应能回答：

1. CDR 方程的每一项在描述什么；
2. 为什么有限元要做分部积分；
3. 一个函数如何变成一组节点系数；
4. $M,K,B,R,F$ 分别来自哪里；
5. 为什么多个物质必须组成块矩阵联立求解；
6. 为什么扩散很小时普通 Galerkin 会振荡；
7. SUPG 解决什么问题，又不解决什么问题；
8. 怎样证明程序结果值得信任。

本讲义允许在数学推导中直接引用 FEALPy 对象，因为项目的目标不是孤立地学习公式，而是理解公式如何进入程序。

## 1. 第一层：一个位置上的收支

### 1.1 从常微分方程开始

先只观察一个位置上的浓度 $u(t)$。最简单的局部模型是

$$
d u_t+cu=f.
$$

移项后：

$$
u_t=\frac{f-cu}{d}.
$$

- $d>0$：容量，决定同样净作用下变化有多快；
- $c>0$：当地消耗系数；
- $f$：外部加入量。

例如 $d=1,c=2,f=30$：

$$
u_t=30-2u.
$$

当 $u=10$ 时，$u_t=10>0$；当 $u=20$ 时，$u_t=-10<0$。平衡值由 $u_t=0$ 得到：

$$
u_*=\frac fc=15.
$$

这个 ODE 是 CDR 方程的局部骨架。空间中的每一个位置都有自己的局部收支，扩散和对流再把不同位置连接起来。

### 1.2 常见困惑：左边的正号是否表示增加

不是。方程写成 $d u_t+cu=f$ 时，反应项移回变化率一侧是

$$
d u_t=f-cu.
$$

判断增加或减少，应先把 $u_t$ 单独写出，再看右边的符号。

## 2. 第二层：空间中的扩散与对流

### 2.1 扩散通量与净流入

一维扩散通量为

$$
J_{\mathrm{diff}}=-a u_x.
$$

$u_x$ 指向浓度增加方向，负号保证扩散从高浓度流向低浓度。一个小区域中的扩散增加量是左边流入减去右边流出：

$$
-\frac{\partial J_{\mathrm{diff}}}{\partial x}
=\frac{\partial}{\partial x}(a u_x).
$$

多维写成

$$
-\nabla\cdot\boldsymbol J_{\mathrm{diff}}
=\nabla\cdot(a\nabla u).
$$

因此项目方程中的 $\nabla\cdot(a\nabla u)$ 是扩散净流入。

在三个等距点上，

$$
u_{xx}\approx\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}.
$$

- $1,5,1$：二阶差分为负，中间高峰下降；
- $5,1,5$：二阶差分为正，中间低谷上升；
- $20,10,10$：左侧有扩散流入，右侧暂时无浓度差，中间上升。

### 2.2 梯度为零与净流出并不矛盾

对

$$
u(x,y)=10-x^2-y^2,
$$

中心处 $\nabla u=0$，但附近扩散通量

$$
\boldsymbol J_{\mathrm{diff}}=(2x,2y)
$$

向四周张开。其散度为 $4>0$，表示中心周围小区域存在净流出，所以中心浓度下降。

梯度回答“这一点向哪里流”，散度回答“这一点周围总体流入还是流出”。

### 2.3 变扩散系数

$$
\nabla\cdot(a\nabla u)
=a\Delta u+\nabla a\cdot\nabla u.
$$

一维中

$$
(a u_x)_x=a u_{xx}+a_xu_x.
$$

因此变系数扩散不能简化为 $a\Delta u$。系数的空间变化也会改变通量收支。

### 2.4 对流与固定观察点

恒定速度的一维纯对流方程为

$$
u_t+b u_x=0,
$$

解为

$$
u(x,t)=u_0(x-bt).
$$

跟着水团移动时，水团的浓度保持不变；站在固定位置时，不同上游水团依次到达，读数会变化。若浓度向右增大：

- 水向右流，左侧较低浓度水团到达，固定点读数下降；
- 水向左流，右侧较高浓度水团到达，固定点读数上升。

项目使用 $\boldsymbol b\cdot\nabla u$。守恒通量形式与它的关系是

$$
\nabla\cdot(\boldsymbol b u)
=\boldsymbol b\cdot\nabla u+(\nabla\cdot\boldsymbol b)u.
$$

## 3. 完整 CDR 方程与数据

先写收支形式：

$$
d u_t
=\nabla\cdot(a\nabla u)
-\boldsymbol b\cdot\nabla u
-cu+f.
$$

再统一移项：

$$
\boxed{
d u_t-\nabla\cdot(a\nabla u)
+\boldsymbol b\cdot\nabla u+cu=f.
}
$$

含时问题还需要：

$$
u(\boldsymbol x,0)=u_0(\boldsymbol x),
$$

$$
u|_{\partial\Omega}=g(\boldsymbol x,t).
$$

初值说明开始时内部状态，Dirichlet 边界直接规定边缘值。稳态问题令 $u_t=0$；瞬态问题研究到达稳态前的演化。

## 4. 从强形式到弱形式

### 4.1 为什么不能直接代入分片直线

P1 有限元函数在每个单元内是直线。它在单元内部二阶导数为零，在单元连接处斜率跳跃，普通二阶导数不存在。强形式包含二阶扩散算子，不能直接逐点要求分片直线满足它。

### 4.2 分部积分降阶

以

$$
-u''=f,\qquad u(0)=u(1)=0
$$

为例。乘测试函数 $v$ 并积分：

$$
\int_0^1(-u'')v=\int_0^1fv.
$$

分部积分：

$$
\int_0^1u'v'-[u'v]_0^1=\int_0^1fv.
$$

测试函数在 Dirichlet 边界为零，故

$$
\int_0^1u'v'=\int_0^1fv.
$$

完整标量 CDR 弱形式为

$$
\int_\Omega d u_t v
+\int_\Omega a\nabla u\cdot\nabla v
+\int_\Omega(\boldsymbol b\cdot\nabla u)v
+\int_\Omega cuv
=\int_\Omega fv.
$$

弱形式没有把折线变光滑，而是降低了对解的导数要求。测试函数检查其支撑区域内的加权平衡。

## 5. 网格、基函数与节点系数

### 5.1 一个单元不等于一个基函数

一维 P1 线段有左右两个节点，因此有两个局部基函数：

$$
\phi_1(x)=1-x,\qquad\phi_2(x)=x
$$

（这里单元取 $[0,1]$）。一般标量近似为

$$
u_h=U_1\phi_1+U_2\phi_2.
$$

两个系数就是左右节点值。二维 P1 三角形有 3 个局部基函数，三维 P1 四面体有 4 个；P2 还会增加边中点自由度。

### 5.2 FEALPy 对应关系

```python
mesh = IntervalMesh.from_interval_domain([0, 1], nx=2)
space = LagrangeFESpace(mesh, p=1)
```

- `mesh`：节点、单元及几何关系；
- `space`：自由度编号、基函数及其导数；
- `space.basis(...)`：计算 $\phi_j$；
- `space.grad_basis(...)`：计算 $\nabla\phi_j$。

局部自由度在相邻单元共享时会映射到同一个全局自由度，单元小矩阵的贡献在这些位置相加。

## 6. 标量矩阵如何产生

写

$$
u_h=\sum_jU_j(t)\phi_j(x),
$$

依次取 $v=\phi_i$。于是：

$$
M_{ij}=\int_\Omega d\phi_j\phi_i,
$$

$$
K_{ij}=\int_\Omega a\nabla\phi_j\cdot\nabla\phi_i,
$$

$$
B_{ij}=\int_\Omega(\boldsymbol b\cdot\nabla\phi_j)\phi_i,
$$

$$
R_{ij}=\int_\Omega c\phi_j\phi_i,
$$

$$
F_i=\int_\Omega f\phi_i.
$$

最终得到

$$
M\dot{\boldsymbol U}+(K+B+R)\boldsymbol U=\boldsymbol F.
$$

矩阵行由测试函数编号决定，矩阵列由未知系数编号决定。对流矩阵通常不对称。

FEALPy 的主要积分器对应为：

```python
ScalarMassIntegrator(coef=d)       # M
ScalarDiffusionIntegrator(coef=a)  # K
ScalarConvectionIntegrator(coef=b) # B
ScalarMassIntegrator(coef=c)       # R
ScalarSourceIntegrator(source=f)   # F
```

## 7. 时间推进

空间离散后，连续空间变量已经进入矩阵条目，剩下

$$
M(t)\dot{\boldsymbol U}(t)+A(t)\boldsymbol U(t)=\boldsymbol F(t).
$$

后向欧拉：

$$
(M^{n+1}+\Delta tA^{n+1})\boldsymbol U^{n+1}
=M^{n+1}\boldsymbol U^n+\Delta t\boldsymbol F^{n+1}.
$$

BDF2：

$$
\left(\frac32M^{n+1}+\Delta tA^{n+1}\right)\boldsymbol U^{n+1}
=M^{n+1}\left(2\boldsymbol U^n-\frac12\boldsymbol U^{n-1}\right)
+\Delta t\boldsymbol F^{n+1}.
$$

BDF2 需要两个历史值，所以项目先用一次后向欧拉得到 $\boldsymbol U^1$，之后再使用 BDF2。

## 8. 多物质耦合：从函数到块矩阵

### 8.1 小写函数与大写节点向量

两种物质是两个函数：

$$
\boldsymbol u_h(x)=
\begin{pmatrix}
u_{1,h}(x)\\u_{2,h}(x)
\end{pmatrix}.
$$

一个 P1 线段有两个空间基函数，所以

$$
u_{1,h}=U_{1,1}\phi_1+U_{1,2}\phi_2,
$$

$$
u_{2,h}=U_{2,1}\phi_1+U_{2,2}\phi_2.
$$

长向量为

$$
\boldsymbol U=
\begin{pmatrix}
U_{1,1}\\U_{1,2}\\U_{2,1}\\U_{2,2}
\end{pmatrix}
=
\begin{pmatrix}\boldsymbol U_1\\\boldsymbol U_2\end{pmatrix}.
$$

一般地，$m$ 个物质、每个物质 $N$ 个自由度，共有 $mN$ 个未知数。

### 8.2 $C$ 为什么变成 $R$

反应矩阵

$$
C=
\begin{pmatrix}
2&-1\\-1&2
\end{pmatrix}
$$

描述同一位置上的

$$
\begin{pmatrix}
2u_1-u_2\\-u_1+2u_2
\end{pmatrix}.
$$

有限元还要对空间基函数积分。令

$$
S_{ij}=\int_\Omega\phi_j\phi_i.
$$

那么数字 $2$ 扩展成空间块 $2S$，数字 $-1$ 扩展成 $-S$：

$$
R=
\begin{pmatrix}
2S&-S\\-S&2S
\end{pmatrix}.
$$

第一块行是物质 1 的节点方程组，第二块行是物质 2 的节点方程组：

$$
\begin{cases}
2S\boldsymbol U_1-S\boldsymbol U_2=\boldsymbol F_1,\\
-S\boldsymbol U_1+2S\boldsymbol U_2=\boldsymbol F_2.
\end{cases}
$$

两个块同时包含 $\boldsymbol U_1,\boldsymbol U_2$，所以必须联立求解。只有非对角反应块为零时，各物质才可能完全分开。

### 8.3 一般 $m$ 分量

$$
(C\boldsymbol u)_r=\sum_{s=1}^mC_{rs}u_s.
$$

离散反应块为

$$
R^{rs}_{ij}=\int_\Omega C_{rs}\phi_j\phi_i.
$$

第 $r$ 条离散方程：

$$
M\dot{\boldsymbol U}_r
+(K+B)\boldsymbol U_r
+\sum_sR^{rs}\boldsymbol U_s
=\boldsymbol F_r.
$$

记忆方式：块行是方程，块列是未知分量。`BlockForm` 拼接所有块，稀疏求解器一次返回完整长向量。

## 9. 适定性与数值稳定不是同一件事

适定性包括存在、唯一和连续依赖。稳态零边界标量问题的能量为

$$
A(u,u)=\int_\Omega a|\nabla u|^2
+\int_\Omega\left(c-\frac12\nabla\cdot\boldsymbol b\right)u^2.
$$

若 $a\ge a_0>0$，有效反应不太负，则扩散能够控制解。多分量时看

$$
C_s=\frac{C+C^T}{2},
\qquad
Q=C_s-\frac12(\nabla\cdot\boldsymbol b)I.
$$

需要分清：

- 连续方程适定性；
- 离散矩阵条件数；
- 空间离散是否振荡；
- 时间格式稳定域。

大常速度的散度为零，不一定破坏连续适定性，却可能让 $Pe$ 很大，使粗网格离散振荡。

## 10. 强对流与 SUPG

### 10.1 为什么 $a$ 小会成为强对流问题

$$
Pe=\frac{|b|h}{2a}.
$$

$a$ 小并没有增大 $b$，而是削弱了扩散，使对流在单元尺度上占主导。边界层厚度约为 $a/b$；若 $h\gg a/b$，网格无法解析薄层。

一维 P1 Galerkin 的内部方程为

$$
\left(-\frac ah-\frac b2\right)U_{i-1}
+\frac{2a}{h}U_i
+\left(-\frac ah+\frac b2\right)U_{i+1}=0.
$$

当 $Pe>1$，相邻增量比例

$$
\frac{U_{i+1}-U_i}{U_i-U_{i-1}}
=\frac{1+Pe}{1-Pe}
$$

为负，节点增量交替换号，产生非物理锯齿。精确解仍可保持单调，因此这不是物理振荡。

### 10.2 SUPG 做了什么

SUPG 增加

$$
\sum_K\int_K\tau(\boldsymbol b\cdot\nabla v_h)R(u_h).
$$

简单一维 P1 情形中，它等价于把流线方向的有效扩散增加到

$$
a_{\mathrm{eff}}=a+\tau b^2,
$$

降低有效 Péclet 数和振荡风险。完整项目使用包含时间、变系数扩散、反应和源项的残差，以保持一致性。

SUPG 不会增加网格自由度。稳定的节点曲线仍可能无法准确解析比网格更薄的边界层。

## 11. 验证：为什么已知答案仍要计算机再算一次

制造解相当于用标准砝码校准电子秤。先选择真解，再反推出源项、边界和初值；程序只看到题目数据，最后比较数值解与真解。

小代数残差

$$
\frac{\|A\boldsymbol U-\boldsymbol F\|_2}{\max(1,\|\boldsymbol F\|_2)}
$$

只说明求解器准确解出了已装配系统。若扩散符号写反，错误系统仍可能有很小残差。

误差范数为

$$
E_0=\left[\int_\Omega\sum_r(u_r-u_{h,r})^2\right]^{1/2},
$$

$$
E_1=\left[\int_\Omega\sum_r|\nabla u_r-\nabla u_{h,r}|^2\right]^{1/2}.
$$

观测收敛阶：

$$
p\approx\log_2\frac{E(h)}{E(h/2)}.
$$

验证链应包括：残差、边界误差、制造解、空间收敛、时间收敛、多分量耦合、强对流对比和输出读回。

## 12. FEALPy 阅读顺序

建议按责任而不是按文件长度阅读：

1. `cdr/pde.py`：系数、源项、边界、初值和精确解的数据约定；
2. `cdr/solver.py` 的网格与空间：`box_mesh`、`LagrangeFESpace`、`TensorFunctionSpace`；
3. `assemble`：双层分量循环、标量积分器与 `BlockForm`；
4. `DirichletBC`：非零边界如何进入右端；
5. `solve_time`：BE 起步与 BDF2；
6. `cdr/stabilization.py`：SUPG 的求积点数组、残差和局部矩阵；
7. `errors` 与 `export_vtu`：误差计算和结果输出。

阅读每段代码时回答四个问题：输入是什么、数组形状是什么、对应哪一项数学公式、输出进入下一步哪里。

## 13. 常见问题

### 一个单元为什么不只对应一个基函数？

有限元基函数通常对应自由度。P1 线段有两个端点自由度，所以有两个局部基函数；P1 三角形有三个。

### 为什么 $C$ 是小矩阵，$R$ 却很大？

$C$ 只描述物质分量关系；$R$ 还要作用于每个分量的全部空间自由度。$m$ 个分量、每个分量 $N$ 个自由度时，完整矩阵是 $mN\times mN$。

### 能否先求完物质 1，再求物质 2？

非对角反应块不为零时不能，因为两组新时刻未知量互相出现在对方方程中，必须联立求解。

### 线性系统残差很小，是否证明 PDE 正确？

不能。残差只验证求解器是否准确满足装配后的系统，不验证弱形式、系数、符号和边界是否正确。

### SUPG 后没有振荡，是否已经得到精确边界层？

不能这样推断。SUPG 控制稳定性；解析薄层仍需要足够细的网格或更合适的自适应策略。

## 14. 入门自检

1. 在 $d u_t+2u=30$ 中，当前 $u=10$ 时，$u_t$ 是多少？
2. 为什么 $\nabla u=0$ 不能单独推出扩散净变化为零？
3. 变系数扩散为什么不能把 $\nabla\cdot(a\nabla u)$ 写成 $a\Delta u$？
4. P1 线段有几个局部基函数？两种物质时有几个局部向量自由度？
5. $B_{ij}=\int(\boldsymbol b\cdot\nabla\phi_j)\phi_i$ 中，导数作用在哪个基函数上？
6. BDF2 为什么不能直接计算第一个新时刻？
7. 块 $R^{23}$ 表示哪个分量影响哪条方程？
8. $a$ 变小、$b,h$ 不变时，$Pe$ 怎样变化？
9. 小残差能够证明什么，不能证明什么？
10. 为什么制造解测试不等于真实任务本身？

### 参考答案

1. $u_t=10$。
2. 梯度是点上的通量方向，净变化由通量散度决定。
3. 还存在 $\nabla a\cdot\nabla u$。
4. 两个标量基函数；两种物质共四个局部向量自由度。
5. 表示未知解的试探基函数 $\phi_j$。
6. 它需要 $U^n$ 和 $U^{n-1}$ 两个历史值。
7. 第 3 个分量影响第 2 条方程。
8. 增大。
9. 证明已装配线性系统被准确求解，不能单独证明 PDE 实现正确。
10. 制造解用于校准求解器；真实任务通常没有解析答案。
