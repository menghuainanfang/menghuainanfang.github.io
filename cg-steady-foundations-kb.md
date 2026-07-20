---
title: 连续有限元（CG）稳态问题理论基础
status: foundation
scope: CG 有限元 / 稳态椭圆问题 / 弱形式与误差分析
prerequisites:
  - 多元微积分与线性代数基础
  - Sobolev 空间的初步概念（$H^1$, $H_0^1$, $L^2$）
  - 有限元离散的基本概念（单元、节点、基函数、自由度）
related:
  - ./cg-steady-foundations.md
  - ../dg-steady-foundations.md
---


**CG 有限元方法的理论核心是一条闭环：将 PDE 转化为弱形式 → 用 Lax–Milgram 保证适定性 → 在协调子空间离散 → 用 Céa 引理将求解误差转化为逼近误差 → 用插值理论给出收敛阶。**

整条链条的可信性建立在三个支柱上：

$$
\boxed{
\text{Sobolev 空间的完备性}
\quad\oplus\quad
\text{双线性型的连续性与强制性}
\quad\oplus\quad
\text{有限元空间的逼近能力}.
}
$$

其中，连续性 + 强制性 → 适定性（Lax–Milgram）；协调离散 + Galerkin 正交 → 拟最优性（Céa）；插值理论 → 具体收敛阶（Bramble–Hilbert）。三者缺一不可。

---

## 1. 模型问题：Poisson 方程

考虑经典 Poisson 型边值问题：

$$
\begin{cases}
-\nabla\cdot(\kappa\nabla u)=f, & \text{in }\Omega,\\
u=0, & \text{on }\partial\Omega.
\end{cases}
$$

其中：
- $\Omega\subset\mathbb{R}^d$（$d=2,3$）为有界 Lipschitz 区域；
- $\kappa(x)$ 满足一致椭圆性：$0<\kappa_0\le \kappa(x)\le \kappa_1$；
- $f\in H^{-1}(\Omega)$（常见情况下 $f\in L^2(\Omega)$）。

这个模型是很多复杂系统线性化后的核心子问题，也是 CG 理论最标准的入口。第 12 节会展示它如何嵌入更一般的二阶椭圆框架。

> **📌 贯穿全文的算例：单位正方形上的 Poisson 方程**
> 为便于理解后续各节的概念，建议读者在脑海中保持一个具体例子：
> $$
> -\Delta u = 2\pi^2\sin(\pi x)\sin(\pi y),\quad (x,y)\in[0,1]^2,\quad u|_{\partial\Omega}=0,
> $$
> 其精确解为 $u(x,y)=\sin(\pi x)\sin(\pi y)\in C^\infty(\overline\Omega)$。
> - **均匀网格 + P1 元**：$H^1$ 误差 $\propto h$，$L^2$ 误差 $\propto h^2$
> - **验证方式**：在 $8\times 8, 16\times 16, 32\times 32, 64\times 64$ 网格上计算误差，在 log-log 图上检验斜率
>
> 这个例子贯穿 CG 理论的始终——从弱形式的推导到收敛阶的数值验证。

---

## 2. 函数空间与核心工具

### 2.1 基本空间速查

CG 理论涉及的 Sobolev 空间及其在有限元分析中的角色：

| 空间 | 定义 | 在 CG 中的角色 |
|------|------|---------------|
| $L^2(\Omega)$ | 平方可积函数空间 | 基础空间；$L^2$-误差的度量空间 |
| $H^1(\Omega)$ | $\{v\in L^2: \partial_{x_i}v\in L^2\}$ | 允许非零边界迹，用于 Neumann/Robin 或 lifting |
| $H_0^1(\Omega)$ | $\{v\in H^1: \gamma v=0\}$ $=\overline{C_c^\infty(\Omega)}^{H^1}$ | 齐次 Dirichlet 的试探/测试空间；Poincaré 使 $\|\nabla v\|$ 成为范数 |
| $H^{-1}(\Omega)$ | $(H_0^1(\Omega))'$ | 右端源项的广义空间；允许分布型右端 |
| $H^{1/2}(\partial\Omega)$ | 迹空间 | 边界数据的自然空间 |

### 2.2 三个关键概念的直观理解

**迹算子与 $H_0^1$**：在 $H^1(\Omega)$ 中，函数只有"几乎处处"意义，不能逐点取边界值。迹算子 $\gamma: H^1(\Omega)\to H^{1/2}(\partial\Omega)$ 是连续线性算子，对光滑函数与经典边界限制一致。$H_0^1$ 定义为核 $\{v\in H^1: \gamma v=0\}$，等价于 $C_c^\infty$ 的 $H^1$ 闭包。**直观**：$H_0^1$ 把齐次 Dirichlet 边界条件"编码"进了函数空间本身，不必逐点强加。

**$H^{-1}$ 对偶空间**：$H^{-1}:=(H_0^1)'$，其中的元素是作用在 $H_0^1$ 上的连续线性泛函。若 $f\in L^2$，通过 $\langle f,v\rangle=\int fv$ 嵌入；更一般地，形如 $-\nabla\cdot g$（$g\in L^2(\Omega)^d$）的分布也属于 $H^{-1}$，作用为 $\langle -\nabla\cdot g,v\rangle:=\int g\cdot\nabla v$。**意义**：右端不必是经典函数——分布型源项可统一纳入同一框架。

**三个空间的选用逻辑**：$H_0^1$ 是能量泛函的自然定义域（边界项自动消失），$H^1$ 处理非齐次边界（通过 lifting 转化为齐次问题），$H^{-1}$ 让右端项的来源更灵活。对 Poisson 齐次 Dirichlet 问题，$V=H_0^1$ 是最自然的选择。

### 2.3 核心工具速查

以下不等式和定理在 CG 理论中频繁使用（设 $\Omega$ 为有界 Lipschitz 域）：

| 工具 | 公式 / 结论 | 典型用途 |
|------|------------|---------|
| **Cauchy–Schwarz** | $|\int_\Omega fg|\le \|f\|_{L^2}\|g\|_{L^2}$ | 乘积项上界估计 |
| **Poincaré** | $\|v\|_{L^2}\le C_P\|\nabla v\|_{L^2}$（$v\in H_0^1$） | $H_0^1$ 中 $\|\nabla v\|$ 等价于全范数 |
| **Sobolev 嵌入** | $H^1\hookrightarrow L^q$：$d=2,q<\infty$；$d=3,q\le6$ | 非线性项的 $L^q$ 控制 |
| **迹定理** | $\|\gamma v\|_{L^2(\partial\Omega)}\le C\|v\|_{H^1}$ | 边界积分的有界性 |
| **Young 不等式** | $ab\le\frac{\varepsilon}{2}a^2+\frac{1}{2\varepsilon}b^2,\;\varepsilon>0$ | 能量估计中吸项 |
| **Lax–Milgram** | $a$ 连续+强制 $\Rightarrow$ $\exists!\,u$ 使 $a(u,v)=\ell(v)$，且 $\|u\|\le\frac{1}{\alpha}\|\ell\|$ | 弱解存在唯一性的核心引擎 |

---

## 3. 弱形式与适定性

### 3.1 从强形式到弱形式

对任意测试函数 $v\in H_0^1(\Omega)$，乘以方程两端并在 $\Omega$ 上积分，利用 Green 公式（分部积分），边界项因 $v|_{\partial\Omega}=0$ 消失：

$$
\boxed{
a(u,v)=\ell(v),\quad \forall v\in H_0^1(\Omega),
}
$$

其中

$$
a(u,v)=\int_\Omega \kappa\nabla u\cdot\nabla v\,dx,\qquad
\ell(v)=\langle f,v\rangle_{H^{-1},H_0^1}.
$$



### 3.2 连续性与强制性

**连续性**：存在常数 $M>0$ 使得
$$
|a(w,v)|\le M\|w\|_{H^1}\|v\|_{H^1},\quad \forall w,v\in H_0^1.
$$

**强制性（椭圆性）**：存在 $\alpha>0$ 使得
$$
a(v,v)\ge \alpha\|v\|_{H^1}^2,\quad \forall v\in H_0^1.
$$
（由 $\kappa\ge\kappa_0>0$ 与 Poincaré 不等式联合保证。）

### 3.3 适定性结论

由 **Lax–Milgram 定理**：弱解 $u\in H_0^1(\Omega)$ 存在且唯一，并满足先验估计

$$
\|u\|_{H^1}\le \frac{1}{\alpha}\|f\|_{H^{-1}}.
$$

这意味着解连续依赖于数据——问题的"适定性"得到保证，离散分析才有意义。

---

## 4. CG 离散：协调有限元空间

### 4.1 离散问题

取形状规则（shape-regular）剖分 $\mathcal{T}_h$，定义连续分片多项式空间（$P_k$ 元）：

$$
V_h = \{v_h\in C^0(\overline\Omega): v_h|_K\in P_k(K),\;\forall K\in\mathcal{T}_h\}\cap H_0^1(\Omega).
$$

离散问题：求 $u_h\in V_h$，使

$$
\boxed{
a(u_h,v_h)=\ell(v_h),\quad \forall v_h\in V_h.
}
$$

关键性质：$V_h\subset H_0^1$（**协调元**），双线性型在子空间上仍连续且强制 → 离散问题自动有唯一解。

### 4.2 为什么协调性重要

协调离散（$V_h\subset V$）带来的直接收益：
1. 离散问题的连续性与强制性**免费继承**自连续问题；
2. Galerkin 正交性自然成立（第 5 节）；
3. Céa 引理直接适用（第 5.2 节）。

非协调离散（如 DG、CR 元）需要额外处理一致性误差（Strang 引理，见第 9.3 节）。



---

## 5. Galerkin 正交与 Céa 引理

### 5.1 Galerkin 正交性

设误差 $e=u-u_h$，由于 $V_h\subset V$，在弱形式中取 $v=v_h$ 并相减：

$$
\boxed{
a(e,v_h)=0,\quad \forall v_h\in V_h.
}
$$

**几何解释**：$u_h$ 是 $u$ 在 $V_h$ 上的 $a(\cdot,\cdot)$-正交投影。

### 5.2 Céa 引理（拟最优性）

利用 Galerkin 正交性，对任意 $v_h\in V_h$：

$$
\alpha\|u_h-v_h\|^2 \le a(u_h-v_h,u_h-v_h) = a(u-v_h,u_h-v_h) \le M\|u-v_h\|\|u_h-v_h\|,
$$

消去一个因子并用三角不等式得到：

$$
\boxed{
\|u-u_h\|_{H^1}
\le \frac{M}{\alpha}
\inf_{v_h\in V_h}\|u-v_h\|_{H^1}.
}
$$

**这条不等式是整个 CG 误差分析的核心**：它将"求解误差"转化为"最佳逼近误差"——CG 解在能量范数下近似最优，误差上限仅由有限元空间的逼近能力决定。

---

## 6. 插值误差与收敛阶

### 6.1 标准收敛结果

若真解 $u\in H^{k+1}(\Omega)$，对 $P_k$-CG 在准一致网格上有：

$$
\boxed{
\begin{aligned}
\|u-u_h\|_{H^1(\Omega)} &\le C h^k |u|_{H^{k+1}(\Omega)}, \\[4pt]
\|u-u_h\|_{L^2(\Omega)} &\le C h^{k+1}|u|_{H^{k+1}(\Omega)} \quad\text{(Aubin–Nitsche)}.
\end{aligned}
}
$$

常见结论速查：

| 单元阶数 | $H^1$-误差阶 | $L^2$-误差阶 | 对偶正则性要求 |
|---------|------------|------------|-------------|
| $P_1$（线性元） | $O(h)$ | $O(h^2)$ | $u\in H^2$，域为凸 |
| $P_2$（二次元） | $O(h^2)$ | $O(h^3)$ | $u\in H^3$ |
| $P_k$（$k$ 次元） | $O(h^k)$ | $O(h^{k+1})$ | $u\in H^{k+1}$ + 对偶正则性 |

### 6.2 直觉解释（类比泰勒展开）

在一维区间 $[0,h]$ 上，$u$ 的 $k$ 次泰勒多项式 $p_k$ 满足：

$$
|u(x)-p_k(x)|\le Ch^{k+1}\max|u^{(k+1)}|,\qquad
|(u-p_k)'(x)|\le Ch^{k}\max|u^{(k+1)}|.
$$

推广到多维：Bramble–Hilbert 定理保证 Lagrange 插值 $I_hu$ 与最佳多项式逼近同阶——

$$
\|u-I_hu\|_{H^m}\le Ch^{k+1-m}|u|_{H^{k+1}}\quad(m=0,1).
$$

取 $m=1$ 得 $H^1$ 误差 $O(h^k)$；取 $m=0$ 得 $L^2$ 误差 $O(h^{k+1})$。


---

## 7. 误差分析的完整框架

### 7.1 误差分解（核心思想）

设 $I_hu\in V_h$ 是某个插值/投影，误差可分解为两个正交成分：

$$
u-u_h = (u-I_hu) + (I_hu-u_h) =: \eta + \xi.
$$

| 成分 | 名称 | 控制机制 |
|------|------|---------|
| $\eta$ | **逼近误差** | 有限元空间的插值能力（Bramble–Hilbert） |
| $\xi$ | **离散误差** | Galerkin 方程 + 稳定性（通常 $\|\xi\|\lesssim\|\eta\|$） |

### 7.2 标准先验误差链

将 Céa 引理与插值估计串接，得到 CG 误差分析的标准三步：

$$
\boxed{
\|u-u_h\|_{H^1}
\;\le\;
\underbrace{C\inf_{v_h\in V_h}\|u-v_h\|_{H^1}}_{\text{Céa 拟最优性}}
\;\le\;
\underbrace{C\|u-I_hu\|_{H^1}}_{\text{取 Lagrange 插值}}
\;\le\;
\underbrace{C h^k |u|_{H^{k+1}}}_{\text{Bramble–Hilbert}}.
}
$$

### 7.3 $L^2$ 超收敛（Aubin–Nitsche 对偶技巧）

设对偶问题 $-\nabla\cdot(\kappa\nabla z)=u-u_h$（齐次 Dirichlet）有 $H^2$ 正则性，则：

$$
\|u-u_h\|_{L^2}^2 = a(z,u-u_h) = a(z-I_hz,u-u_h) \lesssim h\|z\|_{H^2}\|u-u_h\|_{H^1},
$$

约去一个 $\|u-u_h\|_{L^2}$ 因子即得高一阶估计。**注意**：这步对区域与系数正则性敏感——非凸域或奇异解下 $L^2$ 阶可能退化。

### 7.4 收敛讨论的前提检查清单

做理论/数值实验时，建议明确以下前提，避免"结论正确但条件省略"：

1. **区域正则性**：凸域 → $H^2$ 正则性；凹角 → 降阶；
2. **系数正则性**：$\kappa$ 分片常数 vs 光滑 → 影响常数与正则性；
3. **网格条件**：shape-regular（单元不畸变） vs quasi-uniform（全局尺寸一致）；
4. **边界条件类型**：Dirichlet（$H_0^1$） / Neumann（$H^1$） / Robin → 影响函数空间与强制性；
5. **解正则性不足**：角点奇异、界面低正则性 → 收敛阶可能显著低于理论预测。

---

## 8. 非齐次 Dirichlet 条件

### 8.1 Lifting 方法

考虑 $u=g$ on $\partial\Omega$（$g\in H^{1/2}(\partial\Omega)$）。解空间为仿射空间 $V_g:=\{v\in H^1:\gamma v=g\}$，不是线性空间。标准策略是"lifting + 齐次化"：

取 $w\in H^1$ 满足 $\gamma w=g$，令 $u=w+\tilde u$，则 $\tilde u\in H_0^1$ 满足

$$
\boxed{
a(\tilde u,v)=\ell(v)-a(w,v),\quad \forall v\in H_0^1.
}
$$

**核心**：非齐次边界 → 右端修正，主算子不变——整个 CG 理论框架无需修改。

### 8.2 离散与误差结构

离散时构造 $w_h$ 使 $\gamma w_h=g_h$，求 $\tilde u_h\in V_h$：$a(\tilde u_h,v_h)=\ell(v_h)-a(w_h,v_h)$，最终 $u_h=w_h+\tilde u_h$。

当 $g_h\neq g$（边界近似）时，误差需用 **Strang 分解**：

$$
|u-u_h|
\lesssim
\underbrace{\inf_{v_h\in V_h}|u-v_h|}_{\text{逼近误差}}
+
\underbrace{\sup_{v_h\in V_h}\frac{|a(u,v_h)-\ell_h(v_h)|}{\|v_h\|}}_{\text{一致性误差}}.
$$

一致性项主要由 $g_h-g$ 与 $w_h-w$ 驱动。**工程准则**：$g_h$ 的逼近阶必须与有限元阶匹配，否则一致性误差将主导并拉低整体收敛阶。



---

## 9. Variational Crimes（变分犯罪）

经典误差分析假设 $a_h=a$, $\ell_h=\ell$。实际实现中，数值积分和几何逼近会导致偏离——称为 **variational crimes**。

### 9.1 数值积分

实际刚度矩阵通过数值积分公式 $Q_K$ 计算：

$$
a_h(u_h,v_h)=\sum_{K}Q_K(\kappa\nabla u_h\cdot\nabla v_h),\quad\text{而非 }\sum_K\int_K\kappa\nabla u_h\cdot\nabla v_h.
$$

若积分公式对次数 $\le m$ 的多项式精确，$P_k$ 元需 $m\ge 2k-2$ 才能维持 $O(h^k)$。

| 单元 | 刚度矩阵被积多项式次数 | Gauss 积分最低阶 |
|------|---------------------|----------------|
| $P_1$ | $0$ | 1 点（$m=1$） |
| $P_2$ | $2$ | $2\times2$（$m=3$） |
| $P_k$ | $2k-2$ | $\ge 2k-1$（保险） |

### 9.2 曲边界逼近

当 $\Omega$ 有曲边界时，计算域 $\Omega_h$ 为多边形近似，产生几何误差。若 $\text{dist}(\partial\Omega_h,\partial\Omega)=O(h^2)$，$P_1$ 元保持 $O(h)$；高阶元需使用 **isoparametric 元** 使几何逼近阶与有限元阶一致。

### 9.3 Strang 引理

变分犯罪的统一分析工具是 Strang 引理：

$$
\boxed{
\|u-u_h\|_{H^1}
\le
C\left(
\inf_{v_h\in V_h}\|u-v_h\|_{H^1}
+
\sup_{w_h\in V_h}
\frac{|a(v_h,w_h)-a_h(v_h,w_h)|}{\|w_h\|_{H^1}}
\right).
}
$$

**第一项** → 逼近误差；**第二项** → 一致性误差。核心准则：只要第二项保持 $O(h^k)$，整体方法仍最优。

---

## 10. 自适应有限元（AFEM）

### 10.1 基本循环

AFEM 通过后验误差估计驱动局部加密：

$$
\boxed{
\text{SOLVE} \;\rightarrow\; \text{ESTIMATE} \;\rightarrow\; \text{MARK} \;\rightarrow\; \text{REFINE}.
}
$$

### 10.2 残差型误差指示子

最经典的残差型局部误差指示子：

$$
\eta_K^2
=
h_K^2\|f+\nabla\cdot(A\nabla u_h)\|_{L^2(K)}^2
+
\frac12
\sum_{e\subset\partial K}
h_e\|[A\nabla u_h\cdot n]\|_{L^2(e)}^2,
\qquad
\eta^2 = \sum_{K}\eta_K^2.
$$

- 第一项：单元内部残差（$u_h$ 不满足强形式 PDE）；
- 第二项：跨单元法向通量跳跃（$u_h$ 在单元边界不可导）。

### 10.3 可靠性与效率

| 性质 | 不等式 | 含义 |
|------|--------|------|
| **可靠性** | $\|u-u_h\|_{H^1}\le C_{\text{rel}}\,\eta$ | $\eta$ 是误差的安全上界（不会漏报） |
| **效率** | $\eta_K\le C_{\text{eff}}\|u-u_h\|_{H^1(\omega_K)}$ | $\eta_K$ 不过分高估（不会浪费自由度） |

### 10.4 Dörfler 标记策略

给定 $\theta\in(0,1)$，选最小集合 $\mathcal M\subset\mathcal T_h$ 使

$$
\sum_{K\in\mathcal M}\eta_K^2 \ge \theta \sum_{K\in\mathcal T_h}\eta_K^2.
$$

该策略保证每次迭代至少消除固定比例的误差，是 AFEM 收敛理论的基石。

### 10.5 收敛与最优复杂度

存在 $\rho\in(0,1)$ 使得 quasi-error 收缩：

$$
\|u-u_{h_{k+1}}\|_{H^1}^2 + \gamma\eta_{k+1}^2
\le
\rho\left(\|u-u_{h_k}\|_{H^1}^2 + \gamma\eta_k^2\right).
$$

若解属于逼近类 $\mathcal A^s$，则 $\|u-u_h\|_{H^1}\le C N^{-s}$——AFEM 可达到理论最优收敛率。

---

## 11. 偏流主导问题与稳定化

### 11.1 问题与现象

当对流项主导扩散项时（如对流扩散方程 $-\varepsilon\Delta u + \mathbf{b}\cdot\nabla u = f$，$\varepsilon\ll|\mathbf{b}|h$），定义单元 Péclet 数：

$$
\text{Pe}_K = \frac{|\mathbf{b}|h_K}{2\varepsilon}.
$$

$\text{Pe}_K\gg 1$ 时，标准 CG（等价于中心差分）缺乏数值耗散，产生非物理振荡。



### 11.2 SUPG 稳定化

**SUPG（Streamline Upwind Petrov–Galerkin）** 是最常用的稳定化方法，在测试函数方向加入沿流线的扰动：

$$
a(u_h,v_h)
+
\sum_{K}
\tau_K
(\mathbf{b}\cdot\nabla u_h,\;\mathbf{b}\cdot\nabla v_h)_K
=
(f,v_h)
+
\sum_{K}
\tau_K(f,\;\mathbf{b}\cdot\nabla v_h)_K,
\qquad
\tau_K \approx \frac{h_K}{2|\mathbf{b}|}.
$$

### 11.3 稳定化方法速查

| 方法 | 机制 | 适用场景 |
|------|------|---------|
| **SUPG** | 流线方向扰动测试函数 | 对流扩散，中等 Péclet |
| **GLS** | 残差最小二乘稳定化 | 更一般的非对称算子 |
| **CIP** | 梯度跳跃惩罚：$\sum_e\gamma h_e\int_e[\nabla u_h][\nabla v_h]$ | 对流主导 + 希望保持连续性 |
| **DG（对比基准）** | 上风数值通量天然稳定 | 高 Péclet 问题 |

**CG vs DG 的核心差异**：DG 通过数值通量机制天然具备上风稳定性 → 无需额外稳定化项；CG 是中心型离散 → 必须显式引入稳定化。这是 DG 在对流主导问题中的关键优势。

---

## 12. 一般二阶椭圆问题框架

Poisson 方程只是最简单的情形。更一般地，CG 方法的理论基础建立在如下统一模型上：

$$
\begin{cases}
-\nabla\cdot(A(x)\nabla u)+b(x)\cdot\nabla u+c(x)u=f, & x\in\Omega,\\
u=g, & x\in\partial\Omega.
\end{cases}
$$

其中：
- $A(x)\in\mathbb{R}^{d\times d}$ 对称正定，满足一致椭圆性 $\alpha|\xi|^2\le\xi^T A\xi\le\beta|\xi|^2$；
- $b(x)\in\mathbb{R}^d$ 对流项，$c(x)\ge0$ 反应项。

弱形式：求 $u\in V=H_0^1$，使 $a(u,v)=(f,v)$，其中

$$
a(u,v)=\int_\Omega A\nabla u\cdot\nabla v+(b\cdot\nabla u)v+cuv\,dx.
$$

> **📌 Poisson 方程是特例**：当 $A=I$, $b=0$, $c=0$ 时，上述框架退化为第 1 节的 Poisson 问题。因此，前文以 Poisson 方程为载体建立的所有定理（弱形式 → Lax-Milgram → Céa → 误差阶）均可直接迁移到一般框架——只需验证 $a(\cdot,\cdot)$ 的连续性与强制性条件。

---

## 13. 定理链条与导航

### 13.1 核心定理依赖图

CG 理论的所有结果按照以下依赖关系构成一条完整的逻辑链：

```
Sobolev 空间 (H¹, H₀¹, H⁻¹)
  │
  └─→ 弱形式 a(u,v)=ℓ(v)
        │
        ├─ 连续性: |a(w,v)| ≤ M‖w‖‖v‖
        ├─ 强制性: a(v,v) ≥ α‖v‖²
        │
        └─→ Lax–Milgram: ∃! 弱解
              │
              └─→ Galerkin 离散 (V_h ⊂ V)
                    │
                    ├─ 离散可解 (继承连续+强制)
                    ├─ Galerkin 正交: a(e,v_h)=0
                    │
                    └─→ Céa 拟最优性: ‖u-u_h‖ ≤ (M/α) inf‖u-v_h‖
                          │
                          └─→ 插值估计 (Bramble–Hilbert): ‖u-I_h u‖ ≤ Ch^k|u|_{H^{k+1}}
                                │
                                └─→ 先验误差阶: ‖u-u_h‖_{H¹}=O(h^k), ‖u-u_h‖_{L²}=O(h^{k+1})
```

### 13.2 扩展路径

| 扩展方向 | 所需工具 | 产出 |
|---------|---------|------|
| $H^1$ → $L^2$ 超收敛 | Aubin–Nitsche 对偶技巧 + 正则性假设 | $\|u-u_h\|_{L^2}=O(h^{k+1})$ |
| 先验 → 后验 | 残量估计 + 局部单元指标 $\eta_K$ | AFEM 自适应闭环 |
| 协调 → 非协调 | Strang 引理 + 一致性误差估计 | DG / CR 元等非协调方法 |
| 扩散 → 对流主导 | SUPG / GLS / CIP 稳定化 | 高 Péclet 数问题的稳定离散 |
| 线性 → 非线性 | 单调性 / 紧性 / 不动点 / 离散 Brouwer | 非线性椭圆问题的适定与误差分析 |

### 13.3 知识库导航（CG → 多物理耦合）

```text
本文件 (CG 稳态理论基础)
  │
  ├─→ 第 1–7 节：CG 核心骨架（弱形式 → 误差阶）
  │      └─ 映射到 CH/AC 椭圆子步（化学势、扩散算子）
  │
  ├─→ 第 8 节：非齐次边界
  │      └─ 映射到 NS 的 inflow/outflow 边界处理
  │
  ├─→ 第 10 节：AFEM 自适应
  │      └─ 映射到相界面附近局部加密策略
  │
  ├─→ 第 11 节：稳定化
  │      └─ 映射到对流主导的温度/组分传输
  │
  └─→ 第 12 节：一般椭圆框架
         └─ 映射到粘弹性本构的扩散-反应子问题
```
---

## 14. 常见误解

### 误解 1：CG 的 $H^1$ 误差 $O(h^k)$ 和 $L^2$ 误差 $O(h^{k+1})$ 总是成立

不完全正确。$H^1$ 误差阶要求 $u\in H^{k+1}$；$L^2$ 超收敛额外要求对偶问题的 $H^2$ 正则性。在非凸域、角点奇异或系数跳跃时，正则性可能不足，观测到的阶会显著低于理论预测。

### 误解 2：有限元空间阶数越高越好

不正确。高阶元对解的光滑性要求更高。在激波、界面、角点等低正则性区域，高阶元不仅不会提高精度，还可能产生 Gibbs 振荡。$hp$ 混合策略才是应对复杂解结构的正确思路。

### 误解 3：网格越密，误差就越小

原则上正确但不完整。当一致性误差主导时（数值积分精度不足、几何逼近阶不够等），加密网格可能不降误差，甚至因条件数恶化而增大舍入误差。Strang 引理（第 9.3 节）给出了完整的误差结构。

### 误解 4：稳定化只是"加点人工粘性"就行

不准确。稳定化的关键是**一致性**——稳定化项必须随 $h\to0$ 以正确速率消失，否则会污染解的精度。SUPG 的一致性来自沿流线方向的残量加权（而非简单地加 $\varepsilon_{\text{art}}\Delta u_h$）。

### 误解 5：Lax–Milgram 用完就可以忘了

不正确。Lax–Milgram 是 CG 理论的发动机——它不光给出存在唯一性，还给出了稳定性估计 $\|u\|\le\frac{1}{\alpha}\|\ell\|$。这条估计在非线性问题线性化后的每一步、在时间推进的每个子步都会被反复调用。理解它的两个条件（连续性 + 强制性）各自失效的后果，是调试数值崩溃的基本功。

### 误解 6：弱形式只是为了推导方便

不正确。弱形式是解的存在唯一性证明的数学基础（Sobolev 空间框架下的自然形式），也是有限元离散的出发点（Galerkin 方法直接建立在弱形式上）。它不是推导技巧，而是整个 CG 大厦的地基。

---

## 15. 参考文献与阅读顺序

### 15.1 CG 有限元的主教材

1. S. Brenner and R. Scott, *The Mathematical Theory of Finite Element Methods*, 3rd ed., Springer, 2008.
   - CG 数学理论的**首选参考**。建议依次阅读：第 0 章（基础）、第 1–2 章（Sobolev 空间与变分问题）、第 4–5 章（有限元空间与插值）、第 0 章末尾的 Céa 引理。

2. A. Ern and J.-L. Guermond, *Theory and Practice of Finite Elements*, Applied Mathematical Sciences 159, Springer, 2004.
   - 覆盖 CG 与 DG，理论严谨且兼顾实现。第 1–3 章适合打基础；第 4–5 章覆盖混合元与稳定化。

3. T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*, Dover, 2000.
   - **工科入门首选**。物理直觉好，公式推导清晰。适合快速建立有限元的整体图景。

### 15.2 自适应有限元与后验误差

4. M. Ainsworth and J. T. Oden, *A Posteriori Error Estimation in Finite Element Analysis*, Wiley, 2000.
   - 后验误差估计的权威参考。残差型、恢复型估计均有系统处理。

5. W. Bangerth and R. Rannacher, *Adaptive Finite Element Methods for Differential Equations*, Birkhäuser, 2003.
   - 以 DWR（对偶加权残量）目标导向自适应为主线，适合需要控制特定输出量的场景。

### 15.3 Sobolev 空间与 PDE 理论

6. L. C. Evans, *Partial Differential Equations*, 2nd ed., AMS, 2010.
   - 第 5 章（Sobolev 空间）和第 6 章（二阶椭圆方程）是 CG 理论的数学基础。适合需要严谨证明时查阅。

### 15.4 稳定化与对流扩散

7. H.-G. Roos, M. Stynes, and L. Tobiska, *Robust Numerical Methods for Singularly Perturbed Differential Equations*, 2nd ed., Springer, 2008.
   - 对流扩散问题的稳定化方法大全。SUPG、GLS、CIP 及其它方法的理论分析。

---

## 16. 可复用结论

对于后续所有基于 CG 的多物理耦合研究，建议始终先回答以下五个问题：

1. **函数空间是什么？** 试探/测试空间与边界条件是否匹配？$H_0^1$ 还是 $H^1$？需要 lifting 吗？
2. **双线性型是否连续且强制？** 若否——是系数问题、边界条件问题、还是算子结构问题？可否通过稳定化修复？
3. **离散是协调的吗？** 若是 → Céa 引理直接适用；若否 → 需要 Strang 引理，一致性误差是多少？
4. **真解正则性足够吗？** 若正则性不足 → 预期收敛阶降阶多少？可否通过局部 $h$ 加密补偿？
5. **观测到的收敛阶与理论是否一致？** 不一致时——是代码 bug、正则性不足、还是一致性误差主导？

只有这五个问题同时闭合，才能将"一个能跑的程序"提升为"一个具有清楚数学保证和可验证数值性质的 CG 求解器"。

---

## 17. 预备知识学习资源

### 17.1 有限元方法入门

| 教材 | 难度 | 最适合 |
|------|------|--------|
| T. J. R. Hughes, *The Finite Element Method*, Dover | 入门→中级 | **首选推荐**：物理直觉好，公式推导清晰，适合工科背景 |
| K. J. Bathe, *Finite Element Procedures*, 第 2 版 | 入门→高级 | 偏重固体力学和结构，代码细节丰富 |
| S. Brenner & R. Scott, *The Mathematical Theory of FEM*, Springer | 中级→高级 | 偏数学理论，Sobolev 空间和误差分析的处理严谨 |
| Elman, Silvester & Wathen, *Finite Elements and Fast Iterative Solvers*, Oxford | 中级 | 偏流体和迭代求解器，对理解"为什么自适应影响求解器"有用 |

**在线资源**：
- [FEniCSx Tutorial](https://fenicsproject.org/tutorial/) — 用 Python 交互式学习 FEM，每个教程 10–30 分钟
- [deal.II Tutorial](https://www.dealii.org/developer/doxygen/tutorial.html) — C++ 库，教程覆盖从 Poisson 到 Navier-Stokes 的自适应计算
- [NGSolve Tutorials](https://docu.ngsolve.org/latest/) — Python 接口，支持高阶元、混合元与 hybrid 离散

### 17.2 Sobolev 空间与 PDE 理论

| 教材 | 最适合 |
|------|--------|
| L. C. Evans, *Partial Differential Equations*, AMS | 第 5 章：Sobolev 空间的系统介绍；第 6 章：椭圆方程 |
| Brenner & Scott 教材第 0–1 章 | 面向有限元的 Sobolev 空间最小子集 |
| R. A. Adams & J. J. F. Fournier, *Sobolev Spaces*, Academic Press | 纯数学的完整处理，适合深入研究 |

### 17.3 需要复习的数学基础

| 主题 | 推荐快速查阅的资源 |
|------|-------------------|
| 多元微积分（梯度、散度、分部积分/Green 公式） | 任意本科微积分教材的多元部分 |
| 线性代数（正定矩阵、特征值、范数等价） | Strang, *Linear Algebra and Its Applications*, 第 6 章 |
| 泛函分析入门（Hilbert 空间、对偶、Riesz 表示） | Kreyszig, *Introductory Functional Analysis with Applications*, 第 3–4 章 |
| 变分法基础（Euler–Lagrange 方程、最小位能原理） | Gelfand & Fomin, *Calculus of Variations*, 第 1 章 |

---

## 18. 实践工具与代码

### 18.1 有限元与 PDE 求解框架

| 工具 | 语言 | 特点 |
|------|------|------|
| [FEniCSx](https://fenicsproject.org/) | Python/C++ | **快速原型首选**。Python 接口友好，UFL 语言使弱形式直接可写 |
| [deal.II](https://www.dealii.org/) | C++ | 最完善的自适应 FEM 库，tutorial 覆盖 Poisson → NS |
| [NGSolve](https://ngsolve.org/) | Python/C++ | 高阶元、混合元支持好，Python 中直接使用 |
| [MFEM](https://mfem.org/) | C++ | LLNL 开发的并行自适应 FEM 库，支持 $h$、$p$、$hp$ |
| [FreeFEM++](https://freefem.org/) | 自有 DSL | 变分形式直接编写，适合教学和快速探索 |
| [libMesh](https://libmesh.github.io/) | C++ | UT Austin 开发的并行自适应 FEM 库 |

### 18.2 建议的起手实践（CG 稳态问题）

1. **第 1 步**（1–2 小时）：用 FEniCSx 求解贯穿算例（$\sin\pi x\sin\pi y$ 的 Poisson 问题），在 $8\times8$ 到 $64\times64$ 的均匀网格上验证 $P_1$ 元的 $H^1$ 一阶、$L^2$ 二阶收敛。
2. **第 2 步**（2–3 小时）：用 deal.II step-3/step-4 跑通 $P_1$-CG 的完整流程——弱形式 → 组装 → 求解 → 误差计算 → 收敛率检验。
3. **第 3 步**（3–5 小时）：在 L 形域上求解 $-\Delta u=1$（角点奇异），观察 $H^1$ 和 $L^2$ 收敛阶的退化，对比均匀网格与自适应网格（deal.II step-6）的效果。
4. **第 4 步**（5–8 小时）：实现一维对流扩散方程 $-\varepsilon u''+u'=0$（$\varepsilon=10^{-4}$）的 CG 求解，观察 $\text{Pe}\gg1$ 时的振荡，加入 SUPG 稳定化后对比。

---

## 19. 迁移到多物理耦合的模板

将 CG 模板迁移到 CH/AC-NS-热-本构耦合时，每个子模块都应回答以下问题：

| 检查项 | 具体问题 | 对应 CG 理论节点 |
|--------|---------|----------------|
| **函数空间** | 试探/测试空间是什么？与边界条件是否匹配？ | 第 2 节（Sobolev 空间） |
| **适定性** | 双线性/非线性形式是否连续？是否（广义）强制/单调？ | 第 3 节（弱形式与适定性） |
| **离散稳定性** | 离散后能否得到类似 Céa 的拟最优性（或稳定界）？ | 第 5 节（Céa 引理） |
| **误差来源** | 主导项来自逼近、线性化、时间离散、还是一致性？ | 第 7 节（误差分析框架） |
| **收敛验证** | 观测阶是否与理论一致？不一致是否可由正则性不足解释？ | 第 6 节（收敛阶）+ 第 7.4 节（检查清单） |

> **参考实践路线**：在稳态 CG 上完成完整定理链条的数值验证 → 将每个定理映射到 CH/AC/NS/热/本构的对应子模块 → 对比 CG 与 DG 在稳定性、局部守恒与实现复杂度上的差异 → 根据具体物理问题的需求选择或组合离散方案。
