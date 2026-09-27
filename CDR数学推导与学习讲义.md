---
status: reviewing
session_state: completed
date: 2026-09-26
last_updated: 2026-09-26
course: CDR 求解器数学原理
current_section: M9
current_stage: 第一轮数学主线完成
next_action: 需要加深时从适定性证明、SUPG 完整残差或误差范数计算中选择一个专题
last_archived_scope: CDR 数学主线 M1—M9 及与 FEALPy 装配的第一轮衔接
tags:
  - 类型/学习讲义
  - 有限元
  - CDR
---

# 从局部收支到有限元矩阵：CDR 数学学习讲义

这份讲义服务于本项目的实际求解流程。主线不是先罗列符号，而是回答一个问题：一个位置上的量为什么会改变，计算机又怎样把这种变化算出来？

## 学习进度

| ID | 主题 | 当前状态 | 已有理解证据 |
| --- | --- | --- | --- |
| M1 | 局部收支与四种作用 | 已讨论 | 能区分扩散、对流、反应和外部源，并解释容量系数 $d$ |
| M2 | 梯度、通量与散度 | 已讨论 | 能用“净流入”解释扩散，能处理变系数扩散的乘法求导 |
| M3 | 初值、边界与稳态 | 已讨论 | 能判断热量从高温边界向内部扩散，并得到稳态解 $1-x$ |
| M4 | 弱形式与有限元空间 | 已讨论 | 知道分部积分降阶、帽子函数、测试函数和逐行装配 |
| M5 | 标量矩阵与时间推进 | 已讨论 | 能区分 $M,K,B,R$，理解由 $U^n$ 求 $U^{n+1}$ |
| M6 | 多分量反应耦合 | 已讨论 | 能说明函数耦合如何变成节点向量的块矩阵，并理解整体联立求解 |
| M7 | 适定性与能量估计 | 已讨论 | 能区分存在、唯一、连续依赖、条件数与离散振荡；详细证明未检验 |
| M8 | 强对流、振荡与 SUPG | 已讨论 | 能用 $a$ 小、$Pe$ 大、边界层薄解释振荡风险，知道 SUPG 的处理作用 |
| M9 | 误差、收敛阶与可信度 | 已讨论 | 理解残差不能单独证明 PDE 正确，理解制造解是求解器的标准砝码 |

阅读正文不等于已经掌握。表中的“已讨论”表示已经通过问答形成理解，但尚未做完整独立复述或综合练习。

## 1. 从一个位置的收支开始

先不考虑空间，只观察一个位置上的浓度 $u(t)$：

$$
d u_t=-cu+f,
$$

也就是

$$
d u_t+cu=f.
$$

- $u_t$ 是浓度随时间的变化速度；
- $d>0$ 是容量系数，同样的净作用下，$d$ 越大，变化越慢；
- $c>0$ 时，$-cu$ 表示当地消耗；
- $f$ 是外部加入。

例如 $d=1,c=2,f=30$：

$$
u_t=30-2u.
$$

当 $u=10$ 时，$u_t=10>0$，外部加入多于消耗，浓度上升。令 $u_t=0$ 可得平衡值 $u=15$。

这个一阶常微分方程是整个 CDR 模型的骨架。空间中的每一个位置都有一笔这样的账；扩散和对流再把不同位置联系起来。

## 2. 空间中的两种联系

### 2.1 扩散：从高处流向低处

一维扩散通量为

$$
J_{\mathrm{diff}}=-a u_x.
$$

负号表示物质沿浓度降低的方向移动。一个小区域中的扩散净增加量是“流入减去流出”：

$$
-\frac{\partial J_{\mathrm{diff}}}{\partial x}
=\frac{\partial}{\partial x}(a u_x).
$$

多维中写成

$$
-\nabla\cdot\boldsymbol J_{\mathrm{diff}}
=\nabla\cdot(a\nabla u).
$$

因此，$\nabla\cdot(a\nabla u)$ 表示扩散造成的净流入。

当 $a$ 为常数时，

$$
d u_t=a u_{xx}.
$$

在等距的三个位置上，

$$
u_{xx}\approx\frac{u_{\mathrm L}-2u_{\mathrm C}+u_{\mathrm R}}{h^2}.
$$

- 数据 $1,5,1$ 给出负值：中间是高峰，浓度下降；
- 数据 $5,1,5$ 给出正值：中间是低谷，浓度上升。

扩散因此具有“削峰填谷”的效果。

### 2.2 梯度为零不代表没有扩散变化

考虑二维山峰

$$
u(x,y)=10-x^2-y^2.
$$

中心处 $\nabla u=0$，但中心附近的通量箭头整体向外张开：

$$
\boldsymbol J_{\mathrm{diff}}=-\nabla u=(2x,2y),
$$

$$
\nabla\cdot\boldsymbol J_{\mathrm{diff}}=4>0.
$$

这表示中心附近存在净流出，所以中心浓度下降。梯度描述一个点上的流动方向，散度描述周围流量的收支。

### 2.3 变扩散系数

当 $a=a(x)$ 时，

$$
(a u_x)_x=a u_{xx}+a_xu_x.
$$

例如 $a=1+x^2,u=x^2$，则

$$
(a u_x)_x=2+6x^2,
$$

而 $a u_{xx}=2+2x^2$。遗漏的 $4x^2$ 来自 $a_xu_x$。

### 2.4 对流：固定观察点接收上游水团

只有恒定速度 $b$ 的一维对流满足

$$
u_t+b u_x=0,
$$

其解为

$$
u(x,t)=u_0(x-bt).
$$

跟随水团观察时，水团浓度不变；站在固定位置观察时，新的上游水团不断到达，读数会改变。

多维对流项是 $\boldsymbol b\cdot\nabla u$。项目采用非守恒对流形式，而

$$
\nabla\cdot(\boldsymbol b u)
=\boldsymbol b\cdot\nabla u+(\nabla\cdot\boldsymbol b)u.
$$

当速度恒定或 $\nabla\cdot\boldsymbol b=0$ 时，两种形式相同。

## 3. 先写收支形式，再写项目形式

标量 CDR 方程先写成局部收支：

$$
\boxed{
d u_t
=\nabla\cdot(a\nabla u)
-\boldsymbol b\cdot\nabla u
-cu+f.
}
$$

右边依次是扩散净流入、对流贡献、反应贡献和外部源。

把扩散、对流和反应项移到左边，得到项目形式：

$$
\boxed{
d u_t-\nabla\cdot(a\nabla u)
+\boldsymbol b\cdot\nabla u+cu=f.
}
$$

这种顺序可以避免把左边的正负号误解成物质只能流入或只能流出。

## 4. 初值、边界与稳态

初值说明开始时整个区域的状态：

$$
u(\boldsymbol x,0)=u_0(\boldsymbol x).
$$

Dirichlet 边界条件直接指定边界上的值：

$$
u|_{\partial\Omega}=g(\boldsymbol x,t).
$$

例如细杆满足

$$
u_t=u_{xx},\qquad u(0,t)=1,\quad u(1,t)=0.
$$

热量从左端向右扩散。经过足够长时间后，若 $u_t=0$，则 $u_{xx}=0$，结合边界条件得到

$$
u(x)=1-x.
$$

瞬态研究到达稳态之前的过程；稳态研究时间变化停止后的空间分布。

## 5. 为什么要变成弱形式

先看

$$
-u''=f,\qquad u(0)=u(1)=0.
$$

P1 有限元使用相连的小线段。它们在单元内部二阶导数为零，在节点处普通二阶导数又不存在，所以不能直接要求强形式逐点成立。

取边界为零的测试函数 $v$。乘 $v$、积分并分部积分：

$$
\int_0^1(-u'')v\,dx=\int_0^1fv\,dx,
$$

$$
\int_0^1u'v'\,dx-[u'v]_0^1=\int_0^1fv\,dx.
$$

因为 $v(0)=v(1)=0$，得到

$$
\boxed{
\int_0^1u'v'\,dx=\int_0^1fv\,dx.
}
$$

弱形式没有把折线变光滑，而是把对 $u$ 的要求从二阶导数降低到一阶导数。测试函数像一个带权探测器，检查它所关注区域内的加权平衡。

完整标量 CDR 的弱形式是

$$
\boxed{
\int_\Omega d u_t v
+\int_\Omega a\nabla u\cdot\nabla v
+\int_\Omega(\boldsymbol b\cdot\nabla u)v
+\int_\Omega cuv
=\int_\Omega fv.
}
$$

只有扩散项经过分部积分。

## 6. 网格、帽子函数与 FEALPy

把 $[0,1]$ 划分为两个单元：

```text
x₀=0  ─────────  x₁=0.5  ─────────  x₂=1
```

P1 空间在每个小区间上使用直线。每个节点对应一个帽子函数 $\phi_i$，满足

$$
\phi_i(x_j)=
\begin{cases}
1,&i=j,\\
0,&i\ne j.
\end{cases}
$$

近似解写成

$$
u_h(x,t)=\sum_jU_j(t)\phi_j(x).
$$

基函数 $\phi_j(x)$ 描述固定的空间形状；系数 $U_j(t)$ 是待求的节点值。固定网格下，随时间变化的是 $U_j(t)$，不是帽子函数。

这与 FEALPy 的对象直接对应：

```python
mesh = IntervalMesh.from_interval_domain([0, 1], nx=2)
space = LagrangeFESpace(mesh, p=1)
```

- `mesh` 保存节点和单元；
- `space` 保存有限元自由度与基函数规则；
- `space.basis(...)` 计算 $\phi_j$；
- `space.grad_basis(...)` 计算 $\nabla\phi_j$。

选择 $v=\phi_i$ 会得到矩阵的第 $i$ 行；展开 $u_h$ 后，每个 $U_j$ 对应一列。所有单元贡献在共享自由度处相加，这就是装配。

## 7. 标量 CDR 的四类矩阵

将 $u_h=\sum_jU_j\phi_j$ 代入弱形式，并依次取 $v=\phi_i$。

容量质量矩阵：

$$
M_{ij}=\int_\Omega d\phi_j\phi_i,
\qquad M\dot{\boldsymbol U}.
$$

扩散矩阵：

$$
K_{ij}=\int_\Omega a\nabla\phi_j\cdot\nabla\phi_i.
$$

对流矩阵：

$$
B_{ij}=\int_\Omega(\boldsymbol b\cdot\nabla\phi_j)\phi_i.
$$

反应矩阵与载荷：

$$
R_{ij}=\int_\Omega c\phi_j\phi_i,
\qquad
F_i=\int_\Omega f\phi_i.
$$

合起来：

$$
\boxed{
M\dot{\boldsymbol U}+(K+B+R)\boldsymbol U=\boldsymbol F.
}
$$

令 $A=K+B+R$，可写成

$$
M\dot{\boldsymbol U}+A\boldsymbol U=\boldsymbol F.
$$

FEALPy 中的对应关系是：

```python
ScalarMassIntegrator(coef=d)       # M
ScalarDiffusionIntegrator(coef=a)  # K
ScalarConvectionIntegrator(coef=b) # B
ScalarMassIntegrator(coef=c)       # R
ScalarSourceIntegrator(source=f)   # F
```

## 8. 从 $U^n$ 推进到 $U^{n+1}$

空间积分完成后，连续变量 $x$ 已经进入矩阵条目，剩下的是关于节点向量 $\boldsymbol U(t)$ 的常微分方程组。

后向欧拉使用

$$
\dot{\boldsymbol U}^{n+1}
\approx
\frac{\boldsymbol U^{n+1}-\boldsymbol U^n}{\Delta t},
$$

得到

$$
\boxed{
(M^{n+1}+\Delta t A^{n+1})\boldsymbol U^{n+1}
=M^{n+1}\boldsymbol U^n+\Delta t\boldsymbol F^{n+1}.
}
$$

初值给出 $\boldsymbol U^0$。每一步用已知的 $\boldsymbol U^n$ 组成右端，再求 $\boldsymbol U^{n+1}$。

BDF2 是 Backward Differentiation Formula of order 2：

$$
\dot{\boldsymbol U}^{n+1}
\approx
\frac{3\boldsymbol U^{n+1}-4\boldsymbol U^n+\boldsymbol U^{n-1}}{2\Delta t}.
$$

整理后：

$$
\boxed{
\left(\frac32M^{n+1}+\Delta tA^{n+1}\right)\boldsymbol U^{n+1}
=M^{n+1}\left(2\boldsymbol U^n-\frac12\boldsymbol U^{n-1}\right)
+\Delta t\boldsymbol F^{n+1}.
}
$$

第一次推进时只有 $\boldsymbol U^0$，所以项目先用一次后向欧拉得到 $\boldsymbol U^1$，之后再使用 BDF2。

## 9. 从一个物质扩展到多个物质

设有 $m$ 个未知量

$$
\boldsymbol u=(u_1,\ldots,u_m)^T,
$$

反应矩阵为 $C=(C_{rs})$。第 $r$ 条方程中的反应项是

$$
(C\boldsymbol u)_r
=\sum_{s=1}^mC_{rs}u_s.
$$

行 $r$ 表示正在写第 $r$ 条方程，列 $s$ 表示影响来自第 $s$ 个未知量。

每个分量展开为

$$
u_{s,h}=\sum_jU_{s,j}\phi_j.
$$

对第 $r$ 条方程取测试函数 $\phi_i$，定义反应块

$$
R^{rs}_{ij}
=\int_\Omega C_{rs}\phi_j\phi_i.
$$

于是第 $r$ 条离散方程为

$$
M\dot{\boldsymbol U}_r
+(K+B)\boldsymbol U_r
+\sum_{s=1}^mR^{rs}\boldsymbol U_s
=\boldsymbol F_r.
$$

两分量时，空间矩阵为

$$
\begin{pmatrix}
K+B+R^{11}&R^{12}\\
R^{21}&K+B+R^{22}
\end{pmatrix}.
$$

块 $R^{rs}$ 表示第 $s$ 个物质怎样影响第 $r$ 条方程。扩散和对流在本项目中分别作用于各分量，因此 $K+B$ 只出现在块对角线上。

例如

$$
C=
\begin{pmatrix}
1&-1\\
-1&1
\end{pmatrix}
$$

给出

$$
\partial_tu_1=-u_1+u_2,\qquad
\partial_tu_2=u_1-u_2.
$$

两个变化率相加为零，所以没有其他作用时，$u_1+u_2$ 保持不变。

### 9.1 为什么两个物质、一个 P1 线段会有四个未知数

一个 P1 线段单元有左右两个节点，因此有两个标量基函数

$$
\phi_1(x)=1-x,\qquad \phi_2(x)=x.
$$

“一个单元”不等于“一个基函数”。每个物质都要使用这两个基函数：

$$
u_{1,h}=U_{1,1}\phi_1+U_{1,2}\phi_2,
$$

$$
u_{2,h}=U_{2,1}\phi_1+U_{2,2}\phi_2.
$$

所以总未知向量是

$$
\boldsymbol U=
\begin{pmatrix}
U_{1,1}\\U_{1,2}\\U_{2,1}\\U_{2,2}
\end{pmatrix}
=
\begin{pmatrix}
\boldsymbol U_1\\\boldsymbol U_2
\end{pmatrix}.
$$

等价地，向量有限元空间有四个向量基函数：

$$
\begin{pmatrix}\phi_1\\0\end{pmatrix},\quad
\begin{pmatrix}\phi_2\\0\end{pmatrix},\quad
\begin{pmatrix}0\\\phi_1\end{pmatrix},\quad
\begin{pmatrix}0\\\phi_2\end{pmatrix}.
$$

一般地，$m$ 个分量、每个分量 $N$ 个空间自由度，一共有 $mN$ 个未知数。

### 9.2 为什么 $C$ 会扩展成有限元块矩阵

反应矩阵

$$
C=
\begin{pmatrix}
2&-1\\-1&2
\end{pmatrix}
$$

只描述同一空间位置上两个物质的关系。空间弱形式还会产生标量质量型矩阵

$$
S_{ij}=\int_\Omega\phi_j\phi_i.
$$

因此 $C$ 中的每个数会扩展成一个空间块：

$$
R=
\begin{pmatrix}
2S&-S\\-S&2S
\end{pmatrix}.
$$

第一块行是物质 1 的一组节点方程，第二块行是物质 2 的一组节点方程。非对角块不为零时，两组新时刻节点值必须联立求解：

$$
\begin{cases}
2S\boldsymbol U_1-S\boldsymbol U_2=\boldsymbol F_1,\\
-S\boldsymbol U_1+2S\boldsymbol U_2=\boldsymbol F_2.
\end{cases}
$$

FEALPy 的 `BlockForm` 负责拼接这些块，稀疏求解器一次返回所有分量的长向量，而不是先完整求出一个物质再求另一个物质。

## 10. 程序中的装配顺序

```text
PDE 数据 a、b、C、d、f
        ↓
网格 mesh
        ↓
有限元空间 space
        ↓
积分器计算单元矩阵
        ↓
BilinearForm / LinearForm 装配标量矩阵
        ↓
BlockForm 拼接多分量矩阵
        ↓
DirichletBC 处理边界
        ↓
spsolve 求解新时刻的节点向量
```

数学讲义只解释程序对象背后的数学责任。逐行源码、数组形状和 FEALPy 内部调用将在程序学习讲义中单独展开。

## 11. 适定性：有解、唯一并且不失控

“写出了方程”不自动意味着问题可以可靠求解。Hadamard 意义的适定性包含：

1. 存在性：至少有一个解；
2. 唯一性：相同数据不能对应两个不同的解；
3. 连续依赖性：数据小幅变化时，解的变化可以被控制。

对零 Dirichlet 边界的稳态标量问题，取解自身作为测试函数可得

$$
A(u,u)=\int_\Omega a|\nabla u|^2
+\int_\Omega\left(c-\frac12\nabla\cdot\boldsymbol b\right)u^2.
$$

真正进入能量的是有效反应

$$
q=c-\frac12\nabla\cdot\boldsymbol b.
$$

若 $a\ge a_0>0$，并且 $q$ 非负或其负部分没有超过扩散通过 Poincaré 不等式能够控制的程度，则双线性形式有界且强制。Lax–Milgram 定理随后给出唯一弱解和连续依赖性。

多分量问题中，能量看到反应矩阵的对称部分

$$
C_s=\frac{C+C^T}{2},
$$

以及

$$
Q=C_s-\frac12(\nabla\cdot\boldsymbol b)I.
$$

反对称部分不直接贡献 $\boldsymbol u^TC\boldsymbol u$，但仍会改变耦合运动。含时问题还要求 $d\ge d_0>0$；若 $d$ 随时间变化，能量推导中会出现 $d_t$。

这里需要区分三件事：连续问题是否适定、离散矩阵的条件数是否良好、离散解是否保持单调。这三者相关但不等价。大常速度的散度可以为零，不直接破坏连续能量条件，却会增大单元 Péclet 数并诱发离散振荡。

## 12. 强对流、边界层与 SUPG

考虑一维问题

$$
-a u''+b u'=0,\qquad u(0)=0,\quad u(1)=1,\quad b>0.
$$

在一个宽度为 $h$ 的单元上，扩散与对流的尺度分别约为 $a\Delta u/h^2$ 和 $b\Delta u/h$。定义

$$
Pe=\frac{|b|h}{2a}.
$$

当 $a$ 变小，扩散相对减弱、$Pe$ 增大，出口附近形成厚度约为 $a/b$ 的薄边界层。连续精确解仍然单调；非物理锯齿来自粗网格上的中心型 Galerkin 离散。

一维 P1 均匀网格的内部节点满足

$$
\left(-\frac ah-\frac b2\right)U_{i-1}
+\frac{2a}{h}U_i
+\left(-\frac ah+\frac b2\right)U_{i+1}=0.
$$

当 $Pe>1$ 时，右邻点系数变为正，相邻增量的比例

$$
\frac{U_{i+1}-U_i}{U_i-U_{i-1}}
=\frac{1+Pe}{1-Pe}
$$

变为负数，增量交替换号，产生锯齿。项目算例 $a=0.001,b=1,h=1/32$ 给出 $Pe=15.625$。

SUPG 是 Streamline Upwind/Petrov–Galerkin。它在每个单元增加

$$
\int_K\tau(\boldsymbol b\cdot\nabla v_h)R(u_h),
$$

其中 $R(u_h)$ 是完整强残差。简单一维 P1 情形中，新增项等价于沿流线把扩散从 $a$ 增加到

$$
a_{\mathrm{eff}}=a+\tau b^2,
$$

从而降低有效 Péclet 数和振荡风险。项目的完整 SUPG 还包括时间、变系数扩散、反应耦合和源项残差。SUPG 提高稳定性，但不会增加网格自由度；没有振荡不等于薄边界层已经被准确解析。

## 13. 怎样验证计算结果

### 13.1 代数残差只验证求解器

$$
r_{\mathrm{rel}}
=\frac{\|A\boldsymbol U-\boldsymbol F\|_2}
{\max(1,\|\boldsymbol F\|_2)}.
$$

小残差说明线性求解器准确解出了已经装配的系统。如果扩散符号写反，错误系统仍可能具有接近零的残差，所以还必须检查边界、精确解误差和收敛行为。

### 13.2 制造解是求解器的标准砝码

先选择光滑真解 $u$，再由强形式反推出 $f$，并使用相应边界和初值。程序只接收题目数据，计算结果再与真解比较。制造解不是为了替代真实任务，而是为了校准求解器；真实任务通常没有解析答案。

验证包含两类互补测试：

- 精确再现：真解属于离散空间时，误差应接近舍入误差；
- 收敛测试：真解不属于离散空间时，加密网格或缩小时间步，误差应按预期阶数下降。

### 13.3 误差与收敛阶

向量 $L^2$ 误差与梯度误差分别为

$$
E_0=\left[\int_\Omega\sum_r(u_r-u_{h,r})^2\right]^{1/2},
$$

$$
E_1=\left[\int_\Omega\sum_r|\nabla u_r-\nabla u_{h,r}|^2\right]^{1/2}.
$$

若 $E(h)\approx Ch^p$，网格宽度减半时观测阶为

$$
p\approx\log_2\frac{E(h)}{E(h/2)}.
$$

平滑问题的常见预期是 P1 的 $L^2$ 误差二阶、梯度误差一阶，P2 分别三阶和二阶。后向欧拉通常一阶，BDF2 通常二阶。真实问题没有精确解时，还需结合网格加密、边界误差、物理性质和不同方法对照。

## 14. 第一轮数学主线

本项目的数学流程可以压缩为：

```text
局部收支
  → CDR 强形式
  → 初值与边界
  → 弱形式降阶
  → 网格与基函数
  → M、K、B、R、F
  → 多分量块系统
  → BE / BDF2 时间推进
  → 强对流时的 SUPG
  → 残差、制造解与收敛验证
```

第一轮学习已经覆盖这条主线。详细证明、SUPG 参数推导和误差估计可作为后续专题深化，不影响进入程序讲义的学习。
