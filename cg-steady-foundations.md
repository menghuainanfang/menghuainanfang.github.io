# 连续有限元（CG）稳态问题理论基础笔记

> 面向后续多物理耦合（CH/AC-NS-能量-本构）之前的理论打底：先从典型稳态椭圆问题出发，梳理 CG 的函数分析框架、适定性、离散稳定性与收敛性。

---

## 1. 从模型问题开始

我们先考虑经典 Poisson 型边值问题：

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

这个模型是很多复杂系统线性化后的核心子问题，也是 CG 理论最标准的入口。

---

## 2. Sobolev 空间与弱导数（理论准备）

### 2.1 基本空间

- $L^2(\Omega)$：平方可积函数空间；
- $H^1(\Omega)=\{v\in L^2: \partial_{x_i}v\in L^2\}$；
- $H_0^1(\Omega)$：迹为 0 的 $H^1$ 子空间；
- 对偶空间 $H^{-1}(\Omega)=(H_0^1(\Omega))'$。

常用范数：
$$
\|v\|_{H^1(\Omega)}^2=\|v\|_{L^2(\Omega)}^2+\|\nabla v\|_{L^2(\Omega)}^2.
$$
在 $H_0^1$ 中，由 Poincaré 不等式可用半范数等价：
$$
\|v\|_{H^1(\Omega)}\simeq \|\nabla v\|_{L^2(\Omega)}.
$$

### 2.1.1 什么是“迹（trace）”？为什么 $H_0^1$ 要用“迹为 0”定义？

对光滑函数而言，“边界值”就是把 $v(x)$ 限制到 $\partial\Omega$。但在 $H^1(\Omega)$ 中，函数一般只有“几乎处处”意义，不保证逐点可取值，因此需要用函数分析方式定义边界值。

在有界 Lipschitz 区域上，存在连续线性算子
$$
\gamma: H^1(\Omega)\to H^{1/2}(\partial\Omega),
$$
称为**迹算子**，它满足：
1. 对任意光滑函数，$\gamma v$ 与经典边界限制 $v|_{\partial\Omega}$ 一致；
2. 连续性估计
$$
\|\gamma v\|_{H^{1/2}(\partial\Omega)}\le C\|v\|_{H^1(\Omega)}.
$$

其中 $H^{1/2}(\partial\Omega)$ 可以理解为“边界上的半阶 Sobolev 空间”。一个常见定义是用 Gagliardo 型半范数：
$$
|w|_{H^{1/2}(\partial\Omega)}^2
:=\int_{\partial\Omega}\int_{\partial\Omega}
\frac{|w(x)-w(y)|^2}{|x-y|^{d}}\,dS_xdS_y,
$$
并配上
$$
\|w\|_{H^{1/2}(\partial\Omega)}^2
:=\|w\|_{L^2(\partial\Omega)}^2+|w|_{H^{1/2}(\partial\Omega)}^2.
$$
（边界维数是 $d-1$，故分母指数写成 $d=(d-1)+1$。）

**为什么迹算子是线性的？**
因为“取边界值”本质上是限制映射，满足
$$
\gamma(\alpha u+\beta v)=\alpha\gamma u+\beta\gamma v.
$$
在 $H^1$ 中它由光滑函数上的限制映射通过稠密性与连续延拓得到，所以线性被保留。

于是可定义
$$
H_0^1(\Omega)=\{v\in H^1(\Omega):\gamma v=0\text{ on }\partial\Omega\}.
$$
并且有经典等价关系
$$
H_0^1(\Omega)=\overline{C_c^\infty(\Omega)}^{\|\cdot\|_{H^1}}.
$$
其直观理由是：
- $C_c^\infty(\Omega)$ 函数显然在边界为零，因此其闭包包含于“零迹空间”；
- 反过来，若 $\gamma v=0$，可通过延拓、平滑与截断构造一列内部紧支光滑函数在 $H^1$ 逼近 $v$。

**直观理解**：$H_0^1$ 就是“满足齐次 Dirichlet 边界条件”的自然能量空间。这样在弱形式中把边界条件“编码”进试探/测试空间，不需要逐点强加边界值。

### 2.1.2 为什么要区分 $H^1$ 与 $H_0^1$？

- $H^1(\Omega)$：允许非零边界迹，适用于 Neumann/Robin 或非齐次 Dirichlet（配合 lifting）。
- $H_0^1(\Omega)$：边界迹为零，适用于齐次 Dirichlet；可直接用 Poincaré 不等式得到范数等价。

但构造 $H_0^1$ **不只是为了“应付 Poisson 齐次边界”**，更深层原因是：

1. **能量与变分结构的自然定义域**
   很多椭圆算子的能量泛函（如 $\int|\nabla v|^2$）在 $H_0^1$ 上最自然；该空间对一阶导数稳定，且边界项自动消失。

2. **算子论中的核心空间**
   把 $A=-\Delta$ 看作从 $H_0^1$ 到 $H^{-1}$ 的映射时，$H_0^1$ 提供了自然的定义域与闭性框架，是谱理论/半群理论/弱解理论的基石之一。

3. **“零边界”提供规范化（gauge fixing）**
   在纯 Neumann 问题中，解只到常数不唯一；而 $H_0^1$ 通过边界锚定去掉这种自由度，使强制性与唯一性更直接。

4. **离散与连续框架一致**
   有限元空间 $V_h\subset H_0^1$ 时，离散边界条件、能量估计、Céa 引理同属一个 Hilbert 结构，理论与算法强一致。

对于 Poisson 齐次 Dirichlet 问题，选 $V=H_0^1(\Omega)$ 的直接收益是：
1. 双线性型 $a(u,v)=\int\kappa\nabla u\cdot\nabla v$ 在该空间上天然闭合；
2. 借助 Poincaré，不需要额外 $\|u\|_{L^2}$ 控制即可得到强制性；
3. 有限元离散时直接取 $V_h\subset H_0^1$，边界处理与理论一致。

### 2.1.3 对偶空间 $H^{-1}$ 如何理解？

定义
$$
H^{-1}(\Omega):=(H_0^1(\Omega))'.
$$
其中元素是作用在 $H_0^1$ 上的连续线性泛函。其范数写作
$$
\|F\|_{H^{-1}}:=\sup_{0\ne v\in H_0^1}\frac{|\langle F,v\rangle|}{\|v\|_{H^1}}.
$$

这里“右端是普通函数”通常指 $f\in L^2(\Omega)$（甚至更光滑如 $L^p$ 或连续函数），即可以逐点/积分意义直接写 $\int fv$。而在弱框架里，右端可放宽为分布型对象，只要它对 $H_0^1$ 的作用连续即可。

若 $f\in L^2(\Omega)$，可通过
$$
\langle f,v\rangle=\int_\Omega fv\,dx
$$
看成 $H^{-1}$ 元素（因为 $H_0^1\hookrightarrow L^2$ 连续）。更一般地，像 $-\nabla\cdot g$（$g\in L^2(\Omega)^d$）也属于 $H^{-1}$，作用方式为
$$
\langle -\nabla\cdot g,v\rangle:=\int_\Omega g\cdot\nabla v\,dx.
$$

你提到“离散非线性里冻结一些项移到右端”，这通常是**线性化/迭代策略**（Picard、牛顿的某种变体）导致的代数重排；是否能解释弱解存在性，关键在于该右端在每一步是否属于当前测试空间的对偶并满足有界性。也就是说：
- “能放进 $H^{-1}$”是每步弱问题可解的函数分析条件；
- 但整体非线性问题的存在性还需额外工具（单调性、紧性、不动点、离散 Brouwer/Schauder 等）。

**意义**：引入 $H^{-1}$ 后，右端项不必是经典函数，弱导数/分布意义下的源项可统一进入同一框架。

### 2.2 若干关键工具（后续会频繁使用）

- **Cauchy–Schwarz 不等式**；
- **Poincaré 不等式**（控制 $L^2$ 范数）；
- **Sobolev 嵌入定理**（正则性/非线性估计常用）；
- **迹定理**（处理边界条件）；
- **Lax–Milgram 定理**（弱解存在唯一性的核心）。

下面给出常用的“可直接调用版本”（设 $\Omega$ 为有界 Lipschitz 域）：

1. **Cauchy–Schwarz**
$$
\left|\int_\Omega fg\,dx\right|\le \|f\|_{L^2(\Omega)}\|g\|_{L^2(\Omega)}.
$$

2. **Poincaré 不等式**（$v\in H_0^1(\Omega)$）
$$
\|v\|_{L^2(\Omega)}\le C_P\|\nabla v\|_{L^2(\Omega)}.
$$
同时也可写成较弱但常用形式
$$
\|v\|_{L^2(\Omega)}\le C\|v\|_{H^1(\Omega)}.
$$
并可推出
$$
\|v\|_{H^1(\Omega)}\le C\|\nabla v\|_{L^2(\Omega)}.
$$

3. **Sobolev 嵌入（示例）**
- $d=2,3$ 时：$H^1(\Omega)\hookrightarrow L^q(\Omega)$，
  $$
  \begin{cases}
  1\le q<\infty,& d=2,\\
  1\le q\le 6,& d=3,
  \end{cases}
  $$
并有
$$
\|v\|_{L^q}\le C\|v\|_{H^1}.
$$
“嵌入连续”意味着：在 $H^1$ 里有界的函数，其 $L^q$ 范数也统一受控，不会在较弱范数下失控。

直观上，$H^1$ 控制了函数振荡（梯度）与整体大小（$L^2$），所以函数不能“又尖又窄且振幅无限大”；因此可进入某些更强的可积性空间。

4. **迹定理**
$$
\|\gamma v\|_{H^{1/2}(\partial\Omega)}\le C\|v\|_{H^1(\Omega)},\qquad
\|\gamma v\|_{L^2(\partial\Omega)}\le C\|v\|_{H^1(\Omega)}.
$$

5. **Lax–Milgram 定理（简式）**

设 $V$ 为 Hilbert 空间，$a(\cdot,\cdot)$ 在 $V$ 上连续且强制：
$$
|a(w,v)|\le M\|w\|_V\|v\|_V,\qquad a(v,v)\ge \alpha\|v\|_V^2\ (\alpha>0),
$$
且 $\ell\in V'$ 连续，则存在唯一 $u\in V$ 使
$$
a(u,v)=\ell(v),\ \forall v\in V,
$$
并满足
$$
\|u\|_V\le \frac{1}{\alpha}\|\ell\|_{V'}.
$$

6. **Young 不等式（常与 Cauchy 配合）**
$$
ab\le \frac{\varepsilon}{2}a^2+\frac{1}{2\varepsilon}b^2,\quad \varepsilon>0.
$$
常用于把乘积项吸收到能量估计左侧。

---

## 3. 弱形式与适定性

对任意测试函数 $v\in H_0^1(\Omega)$，分部积分得弱形式：

$$
a(u,v)=\ell(v),\quad \forall v\in H_0^1(\Omega),
$$
其中
$$
a(u,v)=\int_\Omega \kappa\nabla u\cdot\nabla v\,dx,
\qquad
\ell(v)=\langle f,v\rangle_{H^{-1},H_0^1}.
$$

### 3.1 连续性
存在常数 $M>0$ 使得
$$
|a(w,v)|\le M\|w\|_{H^1}\|v\|_{H^1},\quad \forall w,v\in H_0^1.
$$

### 3.2 强制性（椭圆性）
存在 $\alpha>0$ 使得
$$
a(v,v)\ge \alpha\|v\|_{H^1}^2,\quad \forall v\in H_0^1.
$$

于是由 **Lax–Milgram**：弱解 $u\in H_0^1(\Omega)$ 存在且唯一，并满足先验估计
$$
\|u\|_{H^1}\le C\|f\|_{H^{-1}}.
$$

---

## 4. CG 离散：有限元空间与离散问题

取形状规则（shape-regular）剖分 $\mathcal{T}_h$，定义连续分片多项式空间（以 $P_k$ 为例）
$$
V_h\subset H_0^1(\Omega),
$$
即标准 CG 空间。

离散问题：求 $u_h\in V_h$，使
$$
a(u_h,v_h)=\ell(v_h),\quad \forall v_h\in V_h.
$$

因为 $V_h\subset H_0^1$，且双线性型在子空间上仍连续且强制，故离散问题也有唯一解。

---

## 5. Galerkin 正交与 Céa 引理

### 5.1 Galerkin 正交性
设误差 $e=u-u_h$，则
$$
a(e,v_h)=0,\quad \forall v_h\in V_h.
$$

### 5.2 Céa 引理（拟最优性）
$u_h$  是 Galerkin 变分问题的离散唯一解 ，满足上述正交性，那么
$$
a(u-u_h , v_h)  = 0 ,\quad \forall v_h \in V_h
$$
可以写做 $u - u_h = u-v_h + v_h - u_h$ 的分裂，代入：
$$
a(v_h - u_h , v_h - u_h)  = a(u- v_h , v_h - u_h)
$$
左侧 $\ge \alpha \|v_h - u_h\|^2$ , 右侧 $\le M\|u-v_h\|\|v_h - u_h\|$  消去和转化得到：
$$
\|u-u_h\|_{H^1}
\le
\frac{M}{\alpha}
\inf_{v_h\in V_h}\|u-v_h\|_{H^1}.
$$

解释：CG 解在能量范数意义下“近似最优”，误差本质上由“空间逼近能力”决定。

---

## 6. 插值误差与收敛阶

若 $u\in H^{k+1}(\Omega)$，对 $P_k$-CG 有（准一致网格下）
$$
\|u-u_h\|_{H^1(\Omega)}\le C h^k |u|_{H^{k+1}(\Omega)}.
$$

结合 Aubin–Nitsche 对偶技巧（需额外正则性）可得
$$
\|u-u_h\|_{L^2(\Omega)}\le C h^{k+1}|u|_{H^{k+1}(\Omega)}.
$$

常见结论（$k=1$, 即 P1 元）：
- $H^1$-误差 $O(h)$；
- $L^2$-误差 $O(h^2)$（在对偶正则性满足时）。

这边我们做一个通俗点的解释，我们可以参考泰勒展开（虽然这不严谨），设 $u \in H^{k+1}(0,h)$  在区间 $[0,h]$  上（先从一维开始考虑。
我们在 $x = 0$  处做泰勒展开：
$$
u(x) = \sum_{j=  0}^k\frac{u^{(j)}(0)}{j!} + \frac{u^{(k+1)}(\xi_x)}{(k+1)!}x^{k+1}
$$
那我们知道前面一部分是 $k$  次多项式，设为 $p_k(x)$ , 则
$$
u(x) - p_k(x) = \frac{u^{(k+1)}(\xi_x)}{(k+1)!}x^{k+1}
$$
显然我们有
$$
|u(x) - p_k(x)| \le Ch^{k+1}
$$
这我们可以类比 $L^2$  误差，那 $H^1$  误差就是求个导数， 
$$
(u - p_k)'(x) = \frac{u^{(k+1)}(\xi_x)}{k!}x^{k}
$$
因此
$$
\|u - p_k\|_{H^1} ∼ h^k
$$
推广到二维，我们设单元直径为 $h$ , 在某一点 $x_0$ 进行泰勒展开：
$$
u(x) = \sum_{|\alpha|\le k}\frac{D^{\alpha}u(x_0)}{\alpha !}(x-x_0)^\alpha + R_{k+1}(x)
$$
其中余项满足
$$
|R_{k+1}(x)|\le Ch^{k+1}\max_{|\beta| = k+1}|D^{\beta} u|
$$
进而也同样有
$$
\|u - p_k\|_{L^2} \le Ch^{k+1}|u|_{H^{k+1}}, \quad \|u - p_k\|_{H^1} \le Ch^{k}|u|_{H^{k+1}}
$$
那我们为什么不使用泰勒多项式而是使用插值来解释？因为 Lagrange 插值 $I_h u$  的局部插值和多项式逼近是同阶误差。这得益于 Bramble-Hibert 定理的解释：
$$
\| u - I_h u\|_{H^m} \le Ch^{k+1 - m}|u|_{H^{k+1}}
$$
这个定理的证明中就使用了泰勒多项式作为中间量联系插值和原函数的。

---

## 7. 收敛性讨论中常见“隐含条件”

做理论时建议明确写出以下前提，避免“结论正确但条件省略”：

1. **区域正则性**：凸域往往更容易获得 $H^2$ 正则性；
2. **系数正则性**：$\kappa$ 是否分片常数/光滑会影响常数与正则性；
3. **网格条件**：shape-regular 与 quasi-uniform 的区别；
4. **边界条件类型**：Dirichlet / Neumann / Robin 对函数空间不同；
5. **解正则性不足时**：阶数可能降阶（例如角点奇异）。

---

## 8. 从稳态到后续问题的迁移

该框架可直接迁移到后续更复杂系统的“单步子问题”：

- CH/AC 中的椭圆子步（化学势、扩散算子）；
- NS 的 Stokes/Poisson 子问题；
- 温度方程隐式步对应椭圆型离散系统；
- 粘弹性本构方程线性化后的扩散-反应型子问题。

建议实践路线：

1. 先在该模型上整理**完整定理链条**（弱形式→Lax-Milgram→Céa→误差阶）；
2. 再将每个定理映射到你后续多物理模型中的对应模块；
3. 最后比较 CG 与 DG 在稳定性、局部守恒与实现复杂度上的差异。

---
## 9. 可继续补充的专题

---

## 10. 为什么必须建立这套理论？（研究动机）

从“会实现算法”到“能判断算法是否可靠”，中间缺的就是这条理论链。对你当前研究方向（多物理耦合、强非线性、可能存在界面与奇异）尤其重要：

1. **弱解框架是复杂模型可解性的最低门槛**
   很多真实问题未必有经典解，但可在 Sobolev 空间中建立弱解；没有这一步，后续离散就失去数学对象。

2. **适定性决定了数值模拟是否“问题本身稳定”**
   若连续问题不适定，离散再精细也无法得到可信结果；Lax–Milgram 提供“存在、唯一、连续依赖数据”的基线保证。

3. **Céa 引理把“求解误差”转化为“逼近误差”**
   它告诉我们：在标准椭圆问题里，CG 的误差上限主要受有限元空间逼近能力支配。这是选择元素阶数、网格策略的理论依据。

4. **误差估计连接理论与实验**
   你在数值实验里观察到的收敛阶（例如 P1 的 $H^1$ 一阶、$L^2$ 二阶）并非经验现象，而是可由定理预测、并用于反向检查代码正确性。

5. **为后续 CH/NS/本构耦合提供模块化分析模板**
   多物理系统虽复杂，但常可拆成若干椭圆/抛物子步，每个子步都能复用“弱形式—稳定性—误差”这套骨架。

---

## 11. 这些定理之间的逻辑关系（建议记忆为一条流水线）

可以把整套 CG 理论记成如下依赖图：

$$
\text{Sobolev 空间}
\Rightarrow \text{弱形式}
\Rightarrow \text{(连续性+强制性)}
\Rightarrow \text{Lax–Milgram 适定性}
\Rightarrow \text{Galerkin 离散可解}
\Rightarrow \text{Galerkin 正交}
\Rightarrow \text{Céa 拟最优性}
\Rightarrow \text{插值估计}
\Rightarrow \text{先验误差阶}.
$$

进一步：

- 若再加 **对偶问题正则性**，就能从 $H^1$-误差推进到 $L^2$-误差（Aubin–Nitsche）；
- 若加入 **残量估计**，可得到后验误差并驱动自适应网格（AFEM）；
- 若放松“完全一致离散”（数值积分、几何近似），需引入 Strang 引理评估一致性误差。

这说明：你现在写的稳态 CG 理论，不是孤立知识点，而是后续高级主题的“母结构”。

---

## 12. 误差分析补充：分解、估计与实践含义

### 12.1 误差分解（核心思想）

设 $I_h u\in V_h$ 是某个插值/投影，误差可分解为
$$
u-u_h=(u-I_hu)+(I_hu-u_h)=:\eta+\xi.
$$

- $\eta$：逼近误差（由空间能力决定）；
- $\xi$：离散误差（由 Galerkin 方程与稳定性控制）。

在对称强制问题中，借助 Galerkin 正交性通常可得
$$
\|\xi\|_{H^1}\lesssim \|\eta\|_{H^1},
$$
因此总体误差与最佳逼近同阶。

### 12.2 Céa + 插值 = 标准先验误差

把 Céa 引理与插值估计组合：
$$
\|u-u_h\|_{H^1}
\le C\inf_{v_h\in V_h}\|u-v_h\|_{H^1}
\le C\|u-I_hu\|_{H^1}
\le C h^k |u|_{H^{k+1}}.
$$

这三步分别对应：
1. 拟最优性（离散稳定）；
2. 选择具体近似函数（构造）；
3. 使用插值理论（逼近能力）。

### 12.3 $L^2$ 误差为何会“高一阶”

设对偶问题
$$
-\nabla\cdot(\kappa\nabla z)=u-u_h
$$
有足够正则性（如 $z\in H^2$），则可把
$$
\|u-u_h\|_{L^2}^2
$$
转写为双线性型并利用 Galerkin 正交性消去低阶项，最终得到高一阶估计：
$$
\|u-u_h\|_{L^2}\lesssim h\|u-u_h\|_{H^1}.
$$

这一步对区域与系数正则性较敏感，因此在非凸域/奇异解下常出现“理论二阶、数值降阶”的现象。

### 12.4 一致性误差：Strang 引理（扩展视角）

若实际离散是
$a_h(\cdot,\cdot),\ell_h(\cdot)$
（例如数值积分、几何逼近），则误差上界常写成两部分：

$$
\text{总误差} \lesssim \underbrace{\text{逼近误差}}_{\text{approximation}} + \underbrace{\text{一致性误差}}_{\text{consistency}}.
$$

这解释了为什么工程实现中“网格很细但误差不降”：可能不是空间阶数问题，而是一致性误差主导。

### 12.5 面向你研究任务的误差分析清单

后续将该模板迁移到 CH/AC-NS-热-本构耦合时，建议每个子模块都回答：

- 采用的函数空间是什么？是否与边界条件匹配？
- 双线性/非线性形式是否满足连续性与（广义）强制/单调性？
- 离散后是否能得到类似 Céa 的拟最优性（或稳定界）？
- 误差主导项来自逼近、线性化、时间离散还是一致性？
- 观测到的收敛阶是否与理论预期一致，若不一致是否可由正则性不足解释？

---

# 13. 非齐次 Dirichlet 条件的处理
考虑
$$
\begin{cases}
-\nabla\cdot(\kappa\nabla u)=f, & x\in\Omega,\\
u=g, & x\in\partial\Omega,
\end{cases}
\qquad g\in H^{1/2}(\partial\Omega).
$$
此时解空间不是 $H_0^1(\Omega)$，而是仿射空间
$$
V_g:=\{v\in H^1(\Omega):\gamma v=g\}.
$$

## 13.1 lifting 与齐次化
取 $w\in H^1(\Omega)$ 满足 $\gamma w=g$，记
$$
u=w+\tilde u,\qquad \tilde u\in H_0^1(\Omega).
$$
则问题等价为：求 $\tilde u\in H_0^1(\Omega)$，使
$$
a(\tilde u,v)=\ell(v)-a(w,v),\qquad \forall v\in H_0^1(\Omega),
$$
其中
$$
a(p,v):=\int_\Omega \kappa\nabla p\cdot\nabla v\,dx,\qquad
\ell(v):=\langle f,v\rangle.
$$
核心点：非齐次边界被转化为右端修正，主算子与齐次情形一致。

## 13.2 能量表述
若
$$
E(v)=\frac12\int_\Omega \kappa|\nabla v|^2\,dx-\int_\Omega fv\,dx,
$$
则原问题是
$$
u=\arg\min_{v\in V_g}E(v),
$$
等价于
$$
\tilde u=\arg\min_{\tilde v\in H_0^1(\Omega)}E(w+\tilde v).
$$
因此 lifting 将“仿射约束最小化”转换为“线性空间最小化”。

## 13.3 有限元离散
取 $V_h\subset H_0^1(\Omega)$，构造离散 lifting $w_h$（边界满足 $\gamma w_h=g_h$），求 $\tilde u_h\in V_h$：
$$
a(\tilde u_h,v_h)=\ell(v_h)-a(w_h,v_h),\qquad \forall v_h\in V_h.
$$
最终
$$
u_h=w_h+\tilde u_h.
$$

lifting 常用两类：插值 lifting（实现简单）与 harmonic lifting（能量更小、通常更稳）。

## 13.4 非齐次边界下的误差结构（Céa/Strang 视角）
写作
$$
u=w+\tilde u,\qquad u_h=w_h+\tilde u_h.
$$
当 $g_h=g$（边界精确表示）时，$\tilde u_h$ 满足标准 Galerkin 正交，直接得到 Céa 型估计：
$$
|u-u_h|_{H^1}
\lesssim
|w-w_h|_{H^1}
+
\inf_{v_h\in V_h}|\tilde u-v_h|_{H^1}.
$$
若 $w_h$ 与有限元阶匹配，整体阶数与齐次问题一致。
当 $g_h\neq g$（边界近似）时出现一致性误差，需用 Strang 分解：
$$
|u-u_h|
\lesssim
\underbrace{\inf_{v_h\in V_h}|u-v_h|}_{\text{逼近误差}}
+
\underbrace{\sup_{v_h\in V_h}\frac{|a(u,v_h)-\ell_h(v_h)|}{\|v_h\|}}_{\text{一致性误差}}.
$$
一致性项主要由 $g_h-g$ 与 $w_h-w$ 驱动。

## 13.5 结论与实用准则
非齐次 Dirichlet 不改变主理论链条：
$$
\text{Sobolev 空间}\rightarrow \text{弱形式}\rightarrow \text{Lax--Milgram}
\rightarrow \text{Galerkin}\rightarrow \text{误差估计}.
$$
其额外成本仅在边界数据处理。工程上应保证 $g_h$ 与有限元阶匹配，否则会引入主导的一致性误差并拉低收敛表现。

# 14. Variational Crimes 对收敛性的影响
经典有限元误差分析假设离散问题与连续变分形式完全一致，但实际计算中常出现 **variational crimes**，即离散双线性形式 $a_h(\cdot,\cdot)$ 与连续形式 $a(\cdot,\cdot)$ 不完全一致。常见来源包括数值积分误差与几何逼近误差，这些误差通过 **一致性误差** 影响最终收敛阶。
## 14.1 数值积分误差
有限元离散需要计算  
$$a(u_h,v_h)=\sum_{K\in\mathcal T_h}\int_K A\nabla u_h\cdot\nabla v_h$$
实际实现中通常采用数值积分公式 $Q_K$ 近似，从而得到离散形式  
$$a_h(u_h,v_h)=\sum_{K\in\mathcal T_h}Q_K(A\nabla u_h\cdot\nabla v_h)$$
若积分公式对次数 $\le m$ 的多项式精确，当 $m\ge 2k-2$ 时，$P_k$ 有限元仍保持最优收敛阶  
$\|u-u_h\|_{H^1(\Omega)}=O(h^k)$。因此实际计算通常选取 **Gauss 积分阶 $\ge 2k-1$** 以避免一致性误差主导。

## 14.2 曲边界逼近
当 $\Omega$ 具有曲边界时，计算区域通常采用多边形 $\Omega_h$ 近似。此时离散双线性形式变为  
$a_h(u_h,v_h)=\int_{\Omega_h}A\nabla u_h\cdot\nabla v_h$，误差主要来自边界位置与法向量的几何近似。
若边界逼近满足 $\text{dist}(\partial\Omega_h,\partial\Omega)=O(h^2)$，则 $P_1$ 元仍保持 $\|u-u_h\|_{H^1}=O(h)$。对于高阶有限元，一般需要使用 **isoparametric 元** 以保证几何逼近阶与有限元阶一致，从而维持 $O(h^k)$ 收敛率。

## 14.3 Strang 误差分解
变分犯罪对误差的影响可通过 Strang 引理描述：
$$
\|u-u_h\|_{H^1}
\le
C\left(
\inf_{v_h\in V_h}\|u-v_h\|_{H^1}
+
\sup_{w_h\in V_h}
\frac{|a(v_h,w_h)-a_h(v_h,w_h)|}{\|w_h\|_{H^1}}
\right).
$$

第一项为有限元空间的逼近误差，第二项为一致性误差。只要后一项保持 $O(h^k)$，整体方法仍具有最优收敛阶。

## 14.4 与 lifting 的关系

在非齐次 Dirichlet 条件的处理过程中，常采用分解 $u=w+G$。实际计算中使用离散 lifting $G_h$，即 $u_h=w_h+G_h$。此时离散问题中既包含 $G_h$ 的逼近误差，也包含由数值积分或几何近似产生的变分犯罪项。

整体误差可理解为三部分的叠加：有限元逼近误差、lifting 逼近误差以及一致性误差。只要 $G_h$ 与离散双线性形式的误差均保持 $O(h^k)$ 量级，则整体解仍满足 $\|u-u_h\|_{H^1}=O(h^k)$。
**总结**  
Variational crimes 通过 Strang 引理中的一致性误差影响有限元收敛性。只要数值积分精度、几何逼近阶以及 lifting 逼近阶与有限元空间阶数匹配，离散方法仍保持最优收敛阶。

# 15. 自适应有限元方法（AFEM）：后验误差估计与收敛理论
自适应有限元方法（Adaptive Finite Element Method, AFEM）通过局部网格加密提高计算效率，其核心思想是根据 **后验误差估计（a posteriori error estimation）** 自动调整网格分辨率，从而在给定计算成本下获得更高精度。
典型 AFEM 算法遵循循环结构
$$
\text{SOLVE} \;\rightarrow\; \text{ESTIMATE} \;\rightarrow\; \text{MARK} \;\rightarrow\; \text{REFINE}.
$$
其中
- **SOLVE**：求解当前网格上的离散问题  
- **ESTIMATE**：计算后验误差指示子  
- **MARK**：选择需要加密的单元  
- **REFINE**：局部细化网格

## 15.1 后验误差估计
考虑椭圆问题
$$
-\nabla\cdot(A\nabla u)=f \quad \text{in } \Omega,
$$
其有限元解为 $u_h\in V_h$。后验误差估计通过计算 **残差型误差指示子** 来衡量局部误差。
常见形式为
$$
\eta_K^2
=
h_K^2\|f+\nabla\cdot(A\nabla u_h)\|_{L^2(K)}^2
+
\frac12
\sum_{e\subset\partial K}
h_e\|[A\nabla u_h\cdot n]\|_{L^2(e)}^2 .
$$
其中
- $h_K$ 为单元尺度  
- $h_e$ 为边长  
- $[\cdot]$ 表示跨单元通量跳跃
全局误差估计量为
$$
\eta^2 = \sum_{K\in\mathcal T_h} \eta_K^2 .
$$

## 15.2 可靠性与效率
理想的后验误差估计需要满足 **可靠性（reliability）** 与 **效率（efficiency）**。
可靠性表示误差估计控制真实误差
$$
\|u-u_h\|_{H^1(\Omega)}
\le
C_{\text{rel}}\eta .
$$
效率表示估计量不严重高估误差
$$
\eta_K
\le
C_{\text{eff}}
\|u-u_h\|_{H^1(\omega_K)} .
$$
其中 $\omega_K$ 为单元 $K$ 的邻域。可靠性保证误差估计的安全性，而效率保证加密策略的有效性。

## 15.3 标记策略（Marking）
误差估计得到局部指标 $\eta_K$ 后，需要确定加密区域。最常用的是 **Dörfler 标记策略**。
给定参数 $\theta\in(0,1)$，寻找最小集合 $\mathcal M\subset\mathcal T_h$ 使
$$
\sum_{K\in\mathcal M}\eta_K^2
\ge
\theta
\sum_{K\in\mathcal T_h}\eta_K^2 .
$$
该策略保证至少消除固定比例的误差，从而推动误差下降。

## 15.4 网格细化
局部加密通常采用保持网格形状规则性的细化策略，例如
- newest vertex bisection (NVB)
- longest edge bisection
- red–green refinement
这些算法保证
- 网格形状正则性  
- 单元尺度满足 $h_{k+1}\le h_k$  
- 加密仅局部发生
从而避免全局自由度爆炸。

## 15.5 AFEM 收敛理论
AFEM 的收敛理论通常建立在 **误差收缩性质** 之上。对于标准残差估计与 Dörfler 标记，可证明存在 $0<\rho<1$ 使
$$
\|u-u_{h_{k+1}}\|_{H^1}^2
+
\gamma\eta_{k+1}^2
\le
\rho
\left(
\|u-u_{h_k}\|_{H^1}^2
+
\gamma\eta_k^2
\right).
$$
该不等式称为 **quasi-error contraction**，其中
$$
\text{quasi-error}
=
\|u-u_h\|_{H^1}^2 + \gamma\eta^2 .
$$
由此可得到
- AFEM 全局收敛  
- 误差单调下降

## 15.6 最优复杂度
进一步可证明 AFEM 具有 **最优逼近复杂度**。若解属于逼近类 $\mathcal A^s$，则
$$
\|u-u_h\|_{H^1}
\le
C N^{-s},
$$
其中 $N$ 为自由度数量。这意味着 AFEM 能在复杂解结构（如边界层、奇异点）附近自动加密，并达到理论最优收敛率。

## 15.7 与移动网格方法的关系
AFEM 与移动网格方法均属于 **自适应网格技术**，但机制不同：
- AFEM：通过 **拓扑改变（refinement）** 调整网格  
- 移动网格：通过 **节点位置变化（mesh redistribution）** 调整网格
两者均依赖误差指标或控制度量，但 AFEM 更侧重离散误差理论，而移动网格通常基于连续度量或能量泛函。

**总结**
AFEM 通过后验误差估计驱动局部网格加密，其理论基础包括：
- 残差型误差估计  
- 可靠性与效率  
- Dörfler 标记策略  
- quasi-error 收缩  
- 最优复杂度
这些结果保证了 AFEM 在复杂解结构问题中具有稳定的收敛性与高效性。


# 16. 稳定化技术：偏流主导问题中的 CG 局限
对于扩散主导问题，标准连续有限元方法（CG）具有良好的稳定性与收敛性。然而当问题 **偏流主导（advection-dominated）** 时，数值解往往出现非物理振荡，这反映了离散系统稳定性的退化。
考虑典型对流扩散方程  
$$-\varepsilon\Delta u + \mathbf{b}\cdot\nabla u = f \quad  in \Omega$$
当扩散系数 $\varepsilon \ll |\mathbf{b}|h$ 时，对流项主导，问题表现出明显的边界层或内部层结构。
此时可引入 **单元 Péclet 数**
$$
\text{Pe}_K = \frac{|\mathbf{b}|h_K}{2\varepsilon}.
$$
当 $\text{Pe}_K \gg 1$ 时，标准 CG 方法通常产生振荡解，即使网格已经相当细。

## 16.1 振荡产生的原因
标准 Galerkin 离散形式为  
$$a(u_h,v_h) = (\varepsilon\nabla u_h,\nabla v_h) + (\mathbf{b}\cdot\nabla u_h,v_h)$$
当扩散项很小时，离散系统主要由对流项控制，而 CG 方法本质上是 **中心差分型离散**，缺乏必要的数值耗散，从而在层结构附近产生振荡。该现象与有限差分中中心差分在高 Péclet 数下失稳的机制类似。
因此需要在离散系统中引入 **稳定化项（stabilization term）**。

## 16.2 SUPG 稳定化
最常见的稳定化方法是 **SUPG（Streamline Upwind Petrov–Galerkin）**。其思想是在测试函数方向加入沿流线方向的扰动。
离散形式写为
$$
a(u_h,v_h)
+
\sum_{K\in\mathcal T_h}
\tau_K
(\mathbf{b}\cdot\nabla u_h,\mathbf{b}\cdot\nabla v_h)_K
=
(f,v_h)
+
\sum_{K\in\mathcal T_h}
\tau_K(f,\mathbf{b}\cdot\nabla v_h)_K .
$$
其中 $\tau_K$ 为稳定化参数，常取
$$
\tau_K \approx \frac{h_K}{2|\mathbf{b}|}.
$$
该方法在流线方向引入适度数值耗散，从而抑制振荡，同时保持较高精度。

## 16.3 其他稳定化方法
除了 SUPG 外，常见稳定化方法还包括：

- **GLS（Galerkin Least Squares）**：通过残差最小二乘引入稳定化项  
- **CIP（Continuous Interior Penalty）**：在梯度跳跃上加入惩罚项  
- **Edge stabilization**：通过边界跳跃增强稳定性
例如 CIP 方法在边界 $e$ 上加入
$$
\sum_{e\in\mathcal E_h}
\gamma h_e
\int_e
[\nabla u_h]\,[\nabla v_h] .
$$
这些方法在保持 CG 空间连续性的同时提供额外数值耗散。

## 16.4 与 DG 方法的关系
相比 CG，DG 方法在单元界面允许解不连续，并通过 **数值通量（numerical flux）** 自然引入上风稳定化。例如在对流问题中使用上风通量即可获得稳定离散。
因此在偏流主导问题中：
- CG 需要额外稳定化项（SUPG、CIP 等）
- DG 则通过数值通量机制天然具备稳定性
这也是 DG 方法在高 Péclet 数问题中具有优势的重要原因。

**总结**
当问题偏流主导时，标准 CG 方法由于缺乏数值耗散会产生非物理振荡。稳定化技术通过在离散系统中引入沿流线方向或界面方向的附加项恢复稳定性，其中 SUPG、GLS 与 CIP 是最常用的方法。

---

# 17. 一般二阶椭圆问题框架
Poisson 方程只是最简单的情形。更一般地，CG 方法通常从如下二阶椭圆边值问题出发：
$$
\begin{cases}
-\nabla\cdot(A(x)\nabla u)+b(x)\cdot\nabla u+c(x)u=f, & x\in\Omega,\\
u=g, & x\in\partial\Omega.
\end{cases}
$$
其中
- $\Omega\subset\mathbb{R}^d$ 为有界 Lipschitz 区域；
- $A(x)\in\mathbb{R}^{d\times d}$ 为对称正定矩阵；
- $b(x)\in\mathbb{R}^d$ 为对流项；
- $c(x)\ge0$ 为反应项；
- $f$ 为给定源项。
### 椭圆性假设
为保证问题适定，需要假设扩散矩阵满足一致椭圆性：
$$
\alpha |\xi|^2
\le
\xi^T A(x)\xi
\le
\beta |\xi|^2,
\qquad
\forall \xi\in\mathbb{R}^d,
$$
其中 $0<\alpha\le\beta$ 为常数。
该条件保证扩散算子在能量意义下是强正定的。
### 弱形式
令试探空间$V=H_0^1(\Omega),$对任意测试函数 $v\in V$，乘以 $v$ 并积分分部得到
$$
a(u,v)=(f,v),
\qquad \forall v\in V,
$$
其中双线性型
$$
a(u,v)
=
\int_\Omega
A\nabla u\cdot\nabla v
+
(b\cdot\nabla u)v
+
cuv
\,dx,
$$
右端泛函
$$
(f,v)=\int_\Omega f v\,dx .
$$
因此弱问题为：
> 求 $u\in V$ 使得  
>$$
a(u,v)=(f,v),\qquad \forall v\in V.$$
### 连续问题的适定性
如果满足：
1. **连续性**
$$
|a(u,v)|\le C\|u\|_{H^1}\|v\|_{H^1},
$$
2. **强制性（coercivity）**
$$
a(v,v)\ge \alpha \|v\|_{H^1}^2 ,
$$
则根据 **Lax–Milgram 定理**，弱问题存在唯一解$u\in H_0^1(\Omega).$
### Poisson 方程作为特例
当
$$
A=I,\quad b=0,\quad c=0
$$
时，上述一般模型退化为
$$
-\Delta u=f,
$$
即前述 Poisson 问题。
因此 Poisson 方程可以看作 **一般二阶椭圆问题的最基本代表**，而 CG 方法的大部分理论均在该统一框架下建立。

