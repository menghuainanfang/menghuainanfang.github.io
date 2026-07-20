# 连续有限元（CG）稳态问题 —— 入门学习计划

> **背景**：大二升大三，信息与计算科学专业，刚加入导师课题组。目标是用 4 周时间从现有基础达到能阅读课题组代码、理解 CG-FEM 核心理论、并动手实现简单求解器的水平。
>
> **本计划基于学长笔记 [`cg-steady-foundations.md`](cg-steady-foundations.md) 的理论体系构建，学习路径与笔记结构一一对应。**
>
> **配套文件**：能力测试题与学习建议见 [`考题与学习建议.md`](考题与学习建议.md)


---

## 第一部分：知识体系全景图

学长笔记的逻辑链是：**Sobolev 空间 → 弱形式 → Lax-Milgram 适定性 → Galerkin 离散 → Céa 引理 → 插值误差 → 收敛阶**。然后向外扩展：非齐次边界、变分犯罪、自适应、稳定化、一般椭圆框架。

### 五层知识结构

```
┌──────────────────────────────────────────────────────────────────┐
│              CG 稳态有限元 —— 知识体系（对应学长笔记）              │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 0：先修基础（补漏）                                    ║ │
│  ║  • 线性代数：矩阵运算、正定性、稀疏矩阵、CSR、条件数          ║ │
│  ║  • 多元微积分：分部积分、散度定理、重积分变量替换              ║ │
│  ║  • ODE：两点边值问题、Green 函数、差分法                      ║ │
│  ║  • 数值方法：高斯求积、Lagrange 插值、分片多项式              ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                          ↓                                        │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 1：函数空间与理论基础 —— 笔记 §2                       ║ │
│  ║  • L², H¹, H₀¹, H⁻¹ 定义与关系                               ║ │
│  ║  • 弱导数、迹算子、Sobolev 嵌入                               ║ │
│  ║  • Poincaré 不等式、Cauchy-Schwarz、Young 不等式              ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                          ↓                                        │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 2：弱形式与适定性 —— 笔记 §3, §17                      ║ │
│  ║  • 从强形式到弱形式的标准推导                                  ║ │
│  ║  • 连续性 + 强制性 → Lax-Milgram → 存在唯一性                 ║ │
│  ║  • 一般二阶椭圆方程的统一框架                                  ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                          ↓                                        │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 3：CG 离散与误差分析 —— 笔记 §4, §5, §6, §12           ║ │
│  ║  • 有限元空间 V_h ⊂ H₀¹（分片多项式）                        ║ │
│  ║  • 离散问题的唯一可解性                                        ║ │
│  ║  • Galerkin 正交性 → Céa 引理（拟最优性）                     ║ │
│  ║  • 插值估计 → H¹ 误差 O(h^k), L² 误差 O(h^{k+1})             ║ │
│  ║  • Aubin-Nitsche 对偶论证                                     ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                          ↓                                        │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 4：扩展专题 —— 笔记 §13–§17                            ║ │
│  ║  • 非齐次 Dirichlet BC（lifting 技术）                        ║ │
│  ║  • Variational Crimes（数值积分、曲边界、Strang 引理）        ║ │
│  ║  • 自适应 FEM（AFEM）：后验估计 → 标记 → 加密 → 收敛          ║ │
│  ║  • 对流主导问题的稳定化（SUPG, GLS, CIP）                     ║ │
│  ║  • 一般椭圆框架 → 多物理耦合的可复用模板                       ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                          ↓                                        │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  Layer 5：编程实现（贯穿全过程）                               ║ │
│  ║  • 1D FEM 从零实现（Python）                                   ║ │
│  ║  • 2D 三角形元 Poisson 求解器                                 ║ │
│  ║  • FEniCSx 入门 → 复现笔记中的典型算例                        ║ │
│  ║  • 收敛阶数值验证（用精确解对比）                              ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

### 前置知识详细清单

#### 线性代数（最优先补漏）

| 知识点 | 在 FEM 中的用途 | 重要程度 |
|--------|----------------|:---:|
| 矩阵乘法、转置、逆 | 刚度矩阵的组装与操作 | ★★★ |
| 对称正定矩阵 | 刚度矩阵的核心性质，保证唯一可解 | ★★★ |
| 稀疏矩阵与 CSR 存储 | 实际 FEM 计算中 99% 以上的元素为零 | ★★★ |
| 条件数 $\kappa(A)$ | 影响求解精度与迭代收敛速度 | ★★☆ |
| 特征值 | 理解条件数、共轭梯度收敛速度 | ★★☆ |
| LU 分解、前代/回代 | 直接法求解 $Ku = f$ 的基础 | ★★☆ |
| 最小二乘法与正规方程 | Galerkin 投影的代数本质 | ★★☆ |

#### 微积分与数学分析

| 知识点 | 在 FEM 中的用途 | 重要程度 |
|--------|----------------|:---:|
| 分部积分 | 弱形式推导的核心工具（把导数从解转移到测试函数）| ★★★ |
| 散度定理 / Green 恒等式 | 多维弱形式推导的必需品 | ★★★ |
| 泰勒展开 | 理解插值误差阶的直接类比 | ★★★ |
| 重积分与变量替换 / Jacobian | 参考单元 → 物理单元的坐标变换 | ★★☆ |
| 函数光滑性（$C^0, C^1$, 弱导数）| 理解为什么 hat function 可以作为 FEM 基函数 | ★★☆ |

#### 常微分方程

| 知识点 | 在 FEM 中的用途 | 重要程度 |
|--------|----------------|:---:|
| 两点边值问题（BVP）| 一维 FEM 的模型问题 $-u''=f$ | ★★★ |
| Green 函数 | Aubin-Nitsche 对偶论证的核心工具 | ★★☆ |
| 有限差分法 | 与 FEM 对比理解离散化的不同思路 | ★★☆ |

#### 数值计算基础

| 知识点 | 在 FEM 中的用途 | 重要程度 |
|--------|----------------|:---:|
| Lagrange 插值 | 单元形函数的构造（一维 P1/P2 元）| ★★★ |
| 分片多项式 | 有限元空间的定义 | ★★★ |
| 高斯求积 | 刚度矩阵与载荷向量的数值积分 | ★★★ |
| 条件数与数值稳定性 | 理解网格很细时求解为何变难 | ★★☆ |
| 多项式逼近误差 | 理解 $H^1$ 误差 $O(h^k)$ 的来源 | ★★☆ |

### 核心理论链（笔记 §11 的流水线，建议背下来）

```
Sobolev 空间 → 弱形式 → (连续性+强制性) → Lax-Milgram 适定性
    → Galerkin 离散可解 → Galerkin 正交 → Céa 拟最优性
    → 插值估计 → 先验误差阶
    → [+对偶正则性] → Aubin-Nitsche → L² 高一阶
    → [+残量估计] → 后验估计 → AFEM
    → [+数值积分/几何近似] → Strang 引理 → 一致性误差
```

---

## 第二部分：四周学习计划

### 总体策略

| 原则 | 说明 |
|------|------|
| **以笔记为纲** | 学长笔记是理论主线，每天的学习对应笔记的具体章节 |
| **先会用、再理解** | 1D FEM 第一周就动手写代码，不要等"学完理论再编程" |
| **理论-代码双向验证** | 每学完一个理论点，立即用代码验证（收敛阶、误差分布等） |
| **线性代数融入全程** | 刚度矩阵组装 = 线性代数；求解 = 线性方程组；误差 = 特征值/条件数 |

### 每周目标

| 周次 | 主题 | 对应笔记 | 核心产出 |
|:---:|------|:---:|---------|
| 第 1 周 | 基础修补 + 弱形式入门 + 1D FEM | §1–§4 | 线性代数过关 + 推导弱形式 + **一维 FEM Python 求解器** |
| 第 2 周 | CG 离散 + 误差分析 | §4–§6, §12 | 理解 Céa 引理 + 数值验证收敛阶 |
| 第 3 周 | 2D FEM + 编程实战 | §4–§6（二维化）| 2D Poisson 求解器 + FEniCSx 入门 |
| 第 4 周 | 扩展专题 + 知识库输出 | §13–§17 | 完整知识库文档 + 课题组代码阅读 |

---

## 第三部分：第 1 周详细日计划

> **目标**：(1) 线性代数与先修知识达标；(2) 理解 Sobolev 空间基本概念；(3) 掌握从 PDE 到弱形式的推导；(4) **从零写出 1D FEM 求解器**。
>
> 配套测试题见 [`考题与学习建议.md`](考题与学习建议.md)，每天下午的练习对应其中的题目编号。

---

### 第 1 天 — 线性代数速通：矩阵与正定性

| 项目 | 内容 |
|------|------|
| **上午 3h** | 复习：矩阵乘法、转置、逆、对称矩阵、正定矩阵的定义与判定（顺序主子式 / $x^T A x > 0$）。阅读 Strang §1.1–§1.6 或任意线性代数教材前 2 章 |
| **下午 2h** | 考题文件：完成线性代数题（正定矩阵判定、局部→全局组装）|
| **晚间编程 1.5h** | **代码 1**：用 NumPy 创建刚度矩阵，验证对称正定性，求解 $Ku = f$ |

```python
# 代码 1：一维 FEM 刚度矩阵初体验
import numpy as np
import matplotlib.pyplot as plt

# 创建 N=5 的一维刚度矩阵
N = 5
h = 1.0 / (N - 1)
K = np.zeros((N, N))
for i in range(N - 1):
    Ke = np.array([[1, -1], [-1, 1]]) / h
    K[i:i+2, i:i+2] += Ke

print("刚度矩阵 K:\n", K)
print("对称?", np.allclose(K, K.T))
eigvals = np.linalg.eigvalsh(K)
print("特征值:", eigvals)
print("全部正定?", np.all(eigvals > 0))
```

---

### 第 2 天 — 线性代数深化：稀疏性与条件数

| 项目 | 内容 |
|------|------|
| **上午 3h** | 学习：稀疏矩阵存储（CSR 格式）、条件数定义 $\kappa(A) = \|A\|\|A^{-1}\|$、共轭梯度法思想 |
| **下午 2h** | 考题文件：完成稀疏矩阵与条件数相关题目 |
| **晚间编程 1.5h** | **代码 2**：比较稠密 vs 稀疏存储，观察不同 $N$ 下条件数的增长 |

```python
# 代码 2：稀疏存储与条件数
import numpy as np
from scipy.sparse import diags

# 用稀疏格式创建 N=100 的刚度矩阵
N = 100
h = 1.0 / (N - 1)
main_diag = 2 * np.ones(N) / h
off_diag = -np.ones(N - 1) / h
K_sparse = diags([off_diag, main_diag, off_diag], [-1, 0, 1], format='csr')

print("非零元素数:", K_sparse.nnz)
print("稠密存储元素数:", N * N)
print("稀疏比:", K_sparse.nnz / (N * N) * 100, "%")

# 观察条件数随 N 增长
for N_test in [10, 50, 100, 500]:
    h_t = 1.0 / (N_test - 1)
    K_t = diags([-np.ones(N_test-1)/h_t, 2*np.ones(N_test)/h_t, -np.ones(N_test-1)/h_t],
                [-1, 0, 1], format='array').toarray()
    K_bc = K_t[1:-1, 1:-1]  # 施加 Dirichlet BC 后
    cond = np.linalg.cond(K_bc)
    print(f"N={N_test:4d}, cond(K)≈{cond:.2e}, N²={N_test**2}")
```

---

### 第 3 天 — 微积分回顾 + Sobolev 空间初识

| 项目 | 内容 |
|------|------|
| **上午 3h** | 阅读笔记 **§2.1–§2.2**：$L^2$, $H^1$, $H_0^1$, 弱导数定义、迹算子的直观理解、Poincaré 不等式。重点理解：**$H_0^1$ = "平方可积且一阶弱导数平方可积 + 边界为零"** |
| **下午 2h** | 考题文件：完成分部积分、散度定理、函数光滑性相关题目 |
| **晚间编程 1.5h** | **代码 3**：数值验证 Poincaré 不等式，直观感受 $H^1$ 函数 vs 非 $H^1$ 函数 |

```python
# 代码 3：Sobolev 空间的数值直觉
import numpy as np
import matplotlib.pyplot as plt

def poincare_ratio(f, df, N=1000):
    """数值计算 ||f||_L² / ||f'||_L²"""
    x = np.linspace(0, 1, N)
    f_vals = f(x)
    df_vals = df(x)
    norm_L2 = np.sqrt(np.trapz(f_vals**2, x))
    norm_H1_semi = np.sqrt(np.trapz(df_vals**2, x))
    return norm_L2 / norm_H1_semi

functions = [
    ("sin(πx)", lambda x: np.sin(np.pi*x), lambda x: np.pi*np.cos(np.pi*x)),
    ("x(1-x)", lambda x: x*(1-x), lambda x: 1-2*x),
    ("sin(2πx)", lambda x: np.sin(2*np.pi*x), lambda x: 2*np.pi*np.cos(2*np.pi*x)),
]

for name, f, df in functions:
    r = poincare_ratio(f, df)
    print(f"函数 {name:12s}: ||f||_L² / ||f'||_L² = {r:.4f} (均 ≤ 1/π ≈ 0.3183)")

# 画出 H¹ 函数 vs 非 H¹ 函数
x = np.linspace(-1, 1, 2000)
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, np.abs(x), label='|x| ∈ H¹')
plt.legend(); plt.title('H¹ 函数')
plt.subplot(1, 3, 2)
plt.plot(x, np.sign(x), label='sign(x) ∉ H¹')
plt.legend(); plt.title('非 H¹ 函数')
plt.subplot(1, 3, 3)
h = 0.001
df_approx = (np.abs(x + h) - np.abs(x - h)) / (2*h)
plt.plot(x, df_approx, label="|x| 的弱导数 ≈ sign(x)")
plt.legend(); plt.title('弱导数（差分近似）')
plt.tight_layout(); plt.savefig('sobolev_intuition.png', dpi=150)
plt.show()
```

---

### 第 4 天 — 弱形式推导 + Lax-Milgram（核心日！）

| 项目 | 内容 |
|------|------|
| **上午 3h** | 阅读笔记 **§1, §3, §17**：Poisson 方程的弱形式推导全过程。**重点：分部积分 → 边界项消失 → 降阶**。理解 Lax-Milgram 定理的条件（连续性 + 强制性）与结论 |
| **下午 2h** | 手推 5 遍：从 $-\nabla \cdot (\kappa \nabla u) = f$ 到 $a(u,v)=\ell(v)$ 的完整推导（含 Green 恒等式）。完成考题文件中的 BVP、能量泛函、Lax-Milgram 相关题目 |
| **晚间编程 1.5h** | **代码 4**：Ritz 法 vs FEM 对比实现，理解不同基函数下的逼近效果 |

```python
# 代码 4：Ritz 法求解 -u'' = 1
import numpy as np
import matplotlib.pyplot as plt

def exact(x):
    return 0.5 * x * (1 - x)

# --- 方法1：全局多项式基 (Ritz) ---
def ritz_global(N_terms):
    """用 sin(nπx) 基函数做 Ritz 近似"""
    K = np.zeros((N_terms, N_terms))
    f = np.zeros(N_terms)
    for m in range(N_terms):
        n = m + 1
        K[m, m] = (n * np.pi)**2 / 2
        f[m] = (1 - (-1)**n) / (n * np.pi)
    c = np.linalg.solve(K, f)
    return c

# --- 方法2：有限元基 (hat functions) ---
def fem_1d_linear(N_elements):
    """一维线性 FEM"""
    N_nodes = N_elements + 1
    h = 1.0 / N_elements
    K = np.zeros((N_nodes, N_nodes))
    F = np.zeros(N_nodes)
    for e in range(N_elements):
        Ke = np.array([[1, -1], [-1, 1]]) / h
        Fe = np.array([h/2, h/2])
        K[e:e+2, e:e+2] += Ke
        F[e:e+2] += Fe
    K_bc = K[1:-1, 1:-1]
    F_bc = F[1:-1]
    u_interior = np.linalg.solve(K_bc, F_bc)
    u = np.zeros(N_nodes)
    u[1:-1] = u_interior
    return u

# 绘图比较
x_fine = np.linspace(0, 1, 200)
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x_fine, exact(x_fine), 'k-', linewidth=2, label='精确解')
for N in [1, 2, 3, 5]:
    c = ritz_global(N)
    u_approx = sum(c[i] * np.sin((i+1)*np.pi*x_fine) for i in range(N))
    plt.plot(x_fine, u_approx, '--', label=f'Ritz N={N}')
plt.legend(); plt.title('Ritz 法（全局基）')
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.plot(x_fine, exact(x_fine), 'k-', linewidth=2, label='精确解')
for Ne in [3, 5, 10]:
    x_nodes = np.linspace(0, 1, Ne+1)
    u_fem = fem_1d_linear(Ne)
    plt.plot(x_nodes, u_fem, 'o-', label=f'FEM N_e={Ne}')
plt.legend(); plt.title('FEM（分片线性基）')
plt.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig('ritz_vs_fem.png', dpi=150)
plt.show()
```

---

### 第 5 天 — 一维 FEM 全流程手算 + 编程（最重要的一天！）

| 项目 | 内容 |
|------|------|
| **上午 3h** | 阅读笔记 **§4**。手算：$[0,1]$ 上 4 等分，求解 $-u''=1$，$u(0)=u(1)=0$。**完整走通：剖分 → 局部刚度矩阵 → 全局组装 → 施加 BC → 求解 → 与精确解比较** |
| **下午 2h** | 考题文件：完成差分法、Lagrange 插值、hat function 相关题目 |
| **晚间编程 1.5h** | **代码 5**：从零实现通用一维 FEM 求解器（这是核心代码，后续每天都用） |

```python
# 代码 5：从零实现一维线性 FEM 求解器（通用版）
import numpy as np
import matplotlib.pyplot as plt

def fem_1d_solver(f_func, kappa_func, a, b, N_elements, u_left, u_right):
    """
    一维 FEM 求解器：-(κ(x)u')' = f(x),  x ∈ (a,b)
                      u(a)=u_left, u(b)=u_right
    """
    N_nodes = N_elements + 1
    h = (b - a) / N_elements
    x_nodes = np.linspace(a, b, N_nodes)

    K = np.zeros((N_nodes, N_nodes))
    F = np.zeros(N_nodes)

    for e in range(N_elements):
        x_left = x_nodes[e]
        x_right = x_nodes[e+1]
        h_e = x_right - x_left

        x_mid = (x_left + x_right) / 2
        kappa_mid = kappa_func(x_mid)

        # 局部刚度矩阵
        Ke = kappa_mid * np.array([[1, -1], [-1, 1]]) / h_e

        # 局部载荷（梯形积分）
        f_left = f_func(x_left)
        f_right = f_func(x_right)
        Fe = np.array([f_left, f_right]) * h_e / 2

        # 组装到全局
        K[e:e+2, e:e+2] += Ke
        F[e:e+2] += Fe

    # 施加 Dirichlet BC（划行划列法）
    K_bc = K[1:-1, 1:-1].copy()
    F_bc = F[1:-1].copy()
    F_bc[0]  -= K[1, 0] * u_left
    F_bc[-1] -= K[-2, -1] * u_right

    u_interior = np.linalg.solve(K_bc, F_bc)
    u = np.zeros(N_nodes)
    u[0] = u_left
    u[-1] = u_right
    u[1:-1] = u_interior

    return x_nodes, u, K

# --- 测试与验证 ---
# 例题: -u'' = 1, u(0)=u(1)=0, 精确解 u = x(1-x)/2
x_fem, u_fem, K_fem = fem_1d_solver(
    f_func=lambda x: 1.0,
    kappa_func=lambda x: 1.0,
    a=0, b=1, N_elements=10,
    u_left=0, u_right=0
)

x_exact = np.linspace(0, 1, 200)
u_exact = 0.5 * x_exact * (1 - x_exact)

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(x_exact, u_exact, 'k-', lw=2, label='精确解')
plt.plot(x_fem, u_fem, 'ro-', ms=6, label='FEM (10单元)')
plt.legend(); plt.grid(True, alpha=0.3)
plt.xlabel('x'); plt.ylabel('u'); plt.title("一维 FEM 求解 -u''=1")

plt.subplot(1, 2, 2)
error = np.abs(u_fem - 0.5 * x_fem * (1 - x_fem))
plt.semilogy(x_fem, error, 'bo-')
plt.xlabel('x'); plt.ylabel('|误差|'); plt.title('节点误差分布')
plt.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig('fem1d_solver.png', dpi=150)
plt.show()

print("最大节点误差:", np.max(error))
```

---

### 第 6 天 — 误差分析：Céa 引理与收敛阶

| 项目 | 内容 |
|------|------|
| **上午 3h** | 阅读笔记 **§5, §6, §12**：Galerkin 正交性、Céa 引理推导、先验误差估计。**核心定理链：Galerkin 正交 → Céa 拟最优 → 插值估计 → $O(h^k)$** |
| **下午 2h** | 考题文件：完成泰勒展开→误差阶、Green 函数→正则性相关题目 |
| **晚间编程 1.5h** | **代码 6**：数值验证收敛阶（log-log 图——理解 FEM 最关键的数字实验） |

```python
# 代码 6：数值验证收敛阶
import numpy as np
import matplotlib.pyplot as plt

def compute_convergence_rates():
    """用不同网格密度求解，验证 H¹ 和 L² 收敛阶"""
    N_list = [4, 8, 16, 32, 64, 128]
    errors_H1 = []
    errors_L2 = []
    h_list = []

    for Ne in N_list:
        x_nodes, u_fem, _ = fem_1d_solver(
            f_func=lambda x: 1.0,
            kappa_func=lambda x: 1.0,
            a=0, b=1, N_elements=Ne,
            u_left=0, u_right=0
        )
        h = 1.0 / Ne
        h_list.append(h)

        u_exact_nodes = 0.5 * x_nodes * (1 - x_nodes)

        # L² 误差（梯形积分）
        error_L2 = np.sqrt(np.trapz((u_fem - u_exact_nodes)**2, x_nodes))
        errors_L2.append(error_L2)

        # H¹ 半范数误差：逐单元比较导数
        error_H1 = 0
        for e in range(Ne):
            he = x_nodes[e+1] - x_nodes[e]
            du_fem = (u_fem[e+1] - u_fem[e]) / he
            xm = (x_nodes[e] + x_nodes[e+1]) / 2
            du_ex = 0.5 - xm
            error_H1 += (du_fem - du_ex)**2 * he
        error_H1 = np.sqrt(error_H1)
        errors_H1.append(error_H1)

    # 拟合收敛阶
    h_arr = np.array(h_list)
    err_H1 = np.array(errors_H1)
    err_L2 = np.array(errors_L2)

    coeff_H1 = np.polyfit(np.log(h_arr), np.log(err_H1), 1)
    coeff_L2 = np.polyfit(np.log(h_arr), np.log(err_L2), 1)

    # 绘图
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.loglog(h_arr, err_H1, 'bo-', label=f'H¹ 误差 (斜率≈{coeff_H1[0]:.2f})')
    plt.loglog(h_arr, err_L2, 'rs-', label=f'L² 误差 (斜率≈{coeff_L2[0]:.2f})')
    plt.loglog(h_arr, h_arr * err_H1[0]/h_arr[0], 'b--', alpha=0.3, label='O(h)')
    plt.loglog(h_arr, h_arr**2 * err_L2[0]/h_arr[0]**2, 'r--', alpha=0.3, label='O(h²)')
    plt.xlabel('h'); plt.ylabel('误差'); plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    plt.title('收敛阶验证（P1 线性元）')

    plt.subplot(1, 2, 2)
    for Ne in [4, 8, 16]:
        xn, uf, _ = fem_1d_solver(lambda x: 1.0, lambda x: 1.0, 0, 1, Ne, 0, 0)
        plt.plot(xn, uf, 'o-', label=f'Ne={Ne}')
    xf = np.linspace(0, 1, 200)
    plt.plot(xf, 0.5*xf*(1-xf), 'k-', lw=2, label='精确解')
    plt.legend(); plt.title('网格加密 → 收敛'); plt.grid(True, alpha=0.3)

    plt.tight_layout(); plt.savefig('convergence_study.png', dpi=150)
    plt.show()

    print(f"P1 线性元收敛阶:")
    print(f"  H¹ 误差: O(h^{coeff_H1[0]:.2f}) (理论: O(h))")
    print(f"  L² 误差: O(h^{coeff_L2[0]:.2f}) (理论: O(h²))")

compute_convergence_rates()
```

---

### 第 7 天 — 周总结 + 知识库框架

| 时间 | 任务 |
|------|------|
| **上午 2h** | 回顾本周笔记与代码；整理未理解透彻的概念清单 |
| **下午 2h** | 搭建个人有限元知识库的 Markdown 框架（按笔记 §1→§17 建目录）；归档本周代码 |
| **晚间** | 休息。若本周目标有遗留，查漏补缺 |

**第 1 周检查清单：**

- [ ] 能独立手推 Poisson 方程的弱形式（5 分钟以内）
- [ ] 能手算 4 单元一维 FEM（从剖分到求解）
- [ ] 能用 Python 从零实现一维 FEM 求解器（不看参考答案）
- [ ] 理解刚度矩阵为何对称正定、稀疏
- [ ] 理解 Lax-Milgram 的条件（连续性 + 强制性）
- [ ] 数值验证了 P1 元的 O(h) 和 O(h²) 收敛阶
- [ ] 能解释 C⁰ 连续性、CG 的含义、弱形式为何降低求导阶数
- [ ] 能解释 Céa 引理（FEM 解是"能量范数下的最佳逼近"）

---

## 第四部分：第 2–4 周概要计划

### 第 2 周：误差分析深化 + 二维 FEM 理论

| 天 | 主题 | 笔记对应 | 核心内容 | 编程任务 |
|:---:|------|:---:|---------|---------|
| 1 | Galerkin 正交 + Céa 引理 | §5 | 手推正交性→拟最优性的完整证明 | 数值验证 Galerkin 正交（$a(u-u_h, v_h)=0$）|
| 2 | 插值误差理论 | §6, §12 | Bramble-Hilbert 引理思想、泰勒展开类比 | 数值观察不同光滑解下的收敛阶差异 |
| 3 | Aubin-Nitsche 对偶 | §12.3 | $L^2$ 为何高一阶？对偶论证的完整推导 | 数值验证角点奇异 → 收敛阶退化 |
| 4 | 二维 Poisson 弱形式 | §3(二维化) | Green 恒等式在二维的推导；三角形剖分 | 用 FEniCS 求解第一个二维 Poisson 问题 |
| 5 | 二维形函数 | §4(二维化) | 面积坐标、线性三角形元形函数、参考单元 | 手动推导三角形元的形函数与局部刚度矩阵 |
| 6 | 二维组装 + 数值积分 | §14.1 | 参考单元上的高斯求积；等参变换；Jacobian | Python 实现三角形元刚度矩阵计算 |
| 7 | 周总结 | — | 二维求解器雏形；整理理论笔记 | **里程碑：运行自己的 2D Poisson 代码** |

### 第 3 周：扩展专题 + 工具实操

| 天 | 主题 | 笔记对应 | 核心内容 |
|:---:|------|:---:|---------|
| 1 | 非齐次 Dirichlet BC | §13 | lifting 技术、插值 lifting vs harmonic lifting |
| 2 | Variational Crimes | §14 | 数值积分误差、曲边界逼近、Strang 引理 |
| 3 | FEniCSx 实战 I | — | 安装 + Tutorial 1-3：Poisson、混合边界 |
| 4 | FEniCSx 实战 II | — | 参数扫描、收敛阶验证、后处理 |
| 5 | AFEM 概念 | §15 | SOLVE→ESTIMATE→MARK→REFINE 循环、Dörfler 标记 |
| 6 | 稳定化方法 | §16 | 对流扩散、Péclet 数、SUPG 直观理解 |
| 7 | 周总结 | — | FEniCS 算例集 + 理论-代码映射表 |

### 第 4 周：知识库输出 + 课题衔接

| 天 | 主题 | 核心内容 |
|:---:|------|---------|
| 1–2 | 撰写个人 FEM 知识库 | 按笔记 §2→§17 的结构，用自己的语言重新阐述每个概念；每个理论点配一段验证代码 |
| 3 | 知识库代码化 | 将前 3 周写的代码整理为可复用的模块（`fem1d.py`, `fem2d.py`, `convergence.py`）|
| 4 | 阅读课题组代码 | 在导师/学长指导下接触课题组 FEM 代码，定位你已理解的模块 |
| 5 | 准备汇报 | 制作 15 分钟 PPT：弱形式→离散→误差→数值验证 |
| 6 | 向课题组汇报 | 获取反馈，调整后续方向 |
| 7 | 前瞻规划 | 从稳态 CG → 时间依赖问题 → 多物理耦合的路线图 |

---

## 第五部分：推荐资源

### 教材

| 优先级 | 书名 | 说明 |
|:---:|------|------|
| ⭐⭐⭐ | **Larson & Bengzon — The Finite Element Method** | 理论与 MATLAB 代码并重，最适合自学 |
| ⭐⭐⭐ | **Johnson — Numerical Solution of PDE by FEM** | 经典薄本，适合快速建立理论框架 |
| ⭐⭐ | **曾攀《有限元分析及应用》** | 中文入门，工程视角 |
| ⭐⭐ | **Strang — Linear Algebra and Its Applications** | 线性代数速查 |
| ⭐ | **Brenner & Scott — Mathematical Theory of FEM** | 进阶理论参考 |

### 在线资源

- **Wolfgang Bangerth (Colorado State) FEM 视频课** (YouTube)：从零讲起，非常清晰
- **FEniCSx 官方教程**：https://fenicsproject.org/tutorial/
- **MIT 18.085/18.086 (Strang)**：计算科学与工程
- **Deal.II 教程**：https://www.dealii.org/developer/doxygen/tutorial.html（虽用 C++，但理论解释极好）

### 编程工具箱

| 工具 | 用途 |
|------|------|
| Python + NumPy + SciPy | 实现 1D/2D FEM 原型 |
| Matplotlib | 可视化（网格、解、误差、收敛图）|
| FEniCSx | 科研级 FEM（Python 接口）|
| Jupyter Notebook | 边写边记录，产出可复现的学习笔记 |

### 环境搭建

```bash
pip install numpy scipy matplotlib
```

验证安装：
```bash
python -c "import numpy; import scipy; import matplotlib; print('环境OK')"
```

---

## 附录：快速参考卡片

### Lax-Milgram 定理速记

> **条件**：$a(\cdot,\cdot)$ 连续（$|a|\le M\|u\|\|v\|$）+ 强制（$a(v,v)\ge \alpha\|v\|^2$），$\ell$ 连续
> **结论**：$\exists! u: a(u,v)=\ell(v)$，且 $\|u\|\le \frac{1}{\alpha}\|\ell\|_{V'}$

### Céa 引理速记

> **核心不等式**：$\|u-u_h\|_{H^1} \le \frac{M}{\alpha} \inf_{v_h\in V_h} \|u-v_h\|_{H^1}$
> **含义**：FEM 解的误差 ≤ 常数 × 最佳逼近误差

### 误差阶速记

| 单元阶数 | $H^1$ 误差 | $L^2$ 误差 |
|:---:|:---:|:---:|
| P1（线性元）| $O(h)$ | $O(h^2)$ |
| P2（二次元）| $O(h^2)$ | $O(h^3)$ |
| Pk | $O(h^k)$ | $O(h^{k+1})$ |

### 一维 FEM 局部刚度矩阵

$$K^e = \frac{1}{h}\begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

### 一步安装

```bash
pip install numpy scipy matplotlib
```

---

> **下一步行动**：
> 1. 安装环境：`pip install numpy scipy matplotlib`
> 2. 做配套考题 [`考题与学习建议.md`](考题与学习建议.md)，定位薄弱环节
> 3. 按第 1 周日计划开始，编程任务**当天完成**——动手写代码是理解 FEM 最快的路径
> 4. 导师或学长有新笔记，随时发给我调整计划
