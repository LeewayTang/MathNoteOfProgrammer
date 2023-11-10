**目录**
 [TOC]


# 1 时间序列分析

**随机时间序列**
按照时间顺序排列的一组随机变量
$$
X_1, X_2, ..., X_t,...
$$
表示一个随机事件的时间序列，记为$\{X_t, t\in T\}$


**观察值序列**
用$x_1, x_2, ..., x_n$表示随机序列的n个有序观察值，称为长度为n的观察值序列。
观察值序列是随机时间序列的一个实现，一次抽样过程。

> 研究目的：揭示随机时间序列的性质
> 研究方法：通过观察值序列进行推断

**时间序列分析方法**
- 描述性时序分析：直观比较、绘图
- 统计时序分析：数理统计学分析
  - 频率分析（频谱分析、谱分析）：假设任何一种无趋势的时间序列都可以分解成若干不同频率的周期波动。早期借助傅里叶分析，然后是傅里叶变换，最后引入了最大熵谱估计理论。是非常有效的动态数据分析方法。
  - 时域分析：从序列自相关的角度揭示时间序列的发展规律，用数学模型进行拟合。

## 1-1 时间序列的预处理
### 1-1-1 平稳性检验

**时间序列的概率分布**

m是正整数，$t_1,t_2,...,t_m\in T$，则m维随机向量$(X_{t_1},X_{t_2},...,X_{t_m})$的联合概率分布是
$$
F_{t_1,t_2,...,t_m}(x_1,x_2,...,x_m)=P(X_{t_1}\le x_1,X_{t_2}\le x_2,...,X_{t_m}\le x_m)
$$
这些有限维的分布函数的全体称为时间序列的概率分布族。
但是分布函数根本得不到……只能研究**序列的低阶矩（特征统计量）**。

**时间序列的均值**
对于时间序列$\{X_t, t\in T\}$，任意时刻的序列值$X_t$是一个随机变量。记分布函数为$F_t(x)$，若
$$
\int^{\infin}_{-\infin}x\mathrm{d}F_t(x)<\infin
$$
则一定存在常数$\mu_t$，它满足
$$
\mu_t=E(X_t)=\int^{\infin}_{-\infin}x\mathrm{d}F_t(x)
$$
称$\mu_t$是序列在t时刻的均值函数。取遍所有的观察时刻，就构成均值序列。

**时间序列的方差**
当$\int^{\infin}_{-\infin}x^2\mathrm{d}F_t(x)<\infin$，定义时间序列的方差函数
$$
\sigma^2_t=D(X_t)=E(X_t-\mu_t)^2=\int^{\infin}_{-\infin}(x-\mu_t)^2\mathrm{d}F_t(x)
$$
取遍所有的观察时刻，就构成方差序列。

**时间序列的自协方差**
衡量时间序列在不同时间点的变化程度（波动性）。
对于时间序列$\{X_t, t\in T\}$，任取$t,s\in T$，自协方差函数定义为：
$$
r(t,s)=E(X_t-\mu_t)(X_s-\mu_s)
$$

**时间序列的自相关系数**
衡量时间序列不同时间点的线性相关程度（随机性）。
对于时间序列$\{X_t, t\in T\}$，任取$t,s\in T$，自相关系数定义为：
$$
\rho(t,s)=\frac{r(t,s)}{\sqrt{D(X_t)\cdot D(X_s)}}
$$

**严平稳**
序列的所有统计性质都不会随着时间变化。

**宽平稳**
序列的低阶矩（二阶）不随时间变化。即
1. $\forall t\in T,E(X^2_t)<\infin$
2. $\forall t\in T,E(X_t)=\mu$
3. $\forall t,s\in T,r(t,s)=r(t+\tau,s+\tau)\Rightarrow D(X_t)=r(t,t)=r(0)$

由性质3定义宽平稳时间序列的延迟k自协方差，即**k阶自协方差函数**$\gamma(k)=\gamma(t,t+k),\forall k\in\Nu$，延迟k自相关系数，即**k阶自相关系数**$\rho_k=\frac{\gamma(t,t+k)}{\sqrt{D(X_t)\cdot D(X_{t+k})}}=\frac{\gamma(k)}{\gamma(0)}$
### 1-1-2 纯随机性检验

**纯随机序列**
序列值之间不相互影响。数学定义是：如果时间序列$\{X_t, t\in T\}$满足性质：
$$
E(X_t)=\mu,\forall t\in T\\
\gamma(t,s)=\begin{cases}
    \sigma^2,&t=s \\
    0,&t\ne s
\end{cases},\forall t,s\in T
$$
则称其为白噪声序列（纯随机序列），记为$X_t\sim WN(\mu,\sigma^2)$
理论上，白噪声序列的自协方差和自相关系数都是0，但是观察值序列有限的时候，可能只是接近0
残差序列具有纯随机性则表明序列信息提取充分。

**方差齐性**
序列的每个变量的方差都相等。
根据Markov定理，只有方差齐性的序列，用最小二乘法得到的未知参数估计值才是准确有效的。

**Barlett定理**
一个纯随机时间序列的n期观察序列$\{x_t,t=1,...,n\}$的延迟k（k不为0）期的样本自相关系数将近似服从均值为0，方差为$\frac 1n$的正态分布。
根据Barlett定理检验纯随机性：
**纯随机性检验**
- 原假设：延迟期数小于等于m的序列值之间相互独立（$H_0:\rho_1=\rho_2=L=\rho_m=0,\forall m\ge 1$）
- 备择假设：延迟期数小于等于m的序列值之间不相互独立（$H_1:\exist \rho_k\ne 0,\forall m\ge 1,k\le m$）
- 检验统计量：由Barlett定理得到Q统计量$Q=n\sum\limits^m_{k=1}\hat{\rho_k}^2\sim\chi^2(m)$，在大样本下更精确的统计量是**QLB统计量**$LB=n(n+2)\sum\limits^m_{k=1}(\frac{\hat{\rho_k}^2}{n-k})\sim\chi^2(m)$
- 判别原则：当Q大于$\chi_{1-\alpha}^2(m)$分位点，或者P(Q)<α，可以以1-α的置信度拒绝原假设，即该序列不是白噪声序列。

## 1-2 分析平稳时间序列
目的：提取平稳序列中的信息。
方法：用平稳时间序列模型拟合目标平稳序列。
### 1-2-1 工具
**p阶差分**
差分就是序列中相邻值作差，反应序列的变化率。
$$
\nabla^px_t=\nabla^{p-1}(\nabla x_t)=\nabla^{p-1}x_t-\nabla^{p-1}x_{t-1}
$$
**k步差分**
$$
\nabla_kx_t=x_t-x_{t-k}
$$
**延迟算子**
相当于将当前序列值的时间回溯一个时刻，记为B。
$$
x_{t-1}=Bx_t \\
... \\
x_{t-p}=B^px_t
$$
**线性差分方程**
可以表征一个序列所包含的特征。
1. 常用于分析自协方差函数和自相关函数；
2. 特征根用于判别稳定性。
$$
z_t+a_1z_{t-1}+...+a_pz_{t-p}=h(t)
$$
当h(t)=0时为齐次。

求解的方式可以参考线性微分方程。如对于齐次线性差分方程求解：
1. 求对应特征方程的特征根。特征方程$\lambda^p+a_1\lambda^{p-1}+...+a_p=0$有p个根，即特征根$\lambda_1...\lambda_p$
2. 求通解，分不同情形：不相等的实根、相等的实根、复数根。

非齐次线性差分方程还需要一个特解。

### 1-2-2 随机时间序列模型
使用一个时间序列的历史值和随机扰动项建立的数学模型$X_t=F(X_{t-1},X_{t-2},...,\mu_t)$，重点：
1. 具体形式
2. 时序变量的滞后期
3. 随机扰动项的结构

平稳时间序列模型有：
**AR模型**
即p阶自回归（AutoRegression）模型，AR(p)的t时刻序列值受到过去p时间段内历史值的影响，以及一个和历史序列值无关的白噪声残差。
$$
\begin{cases}
    x_t=\varphi_0+\varphi_1x_{t-1}+\varphi_2x_{t-2}+...+\varphi_px_{t-p}+\epsilon_t \\
    \varphi_p\ne 0 \\
    \epsilon_t 是白噪声 \\
    E(x_s\epsilon_t)=0,\forall s<t
\end{cases}
$$

$\varphi_0=0$时成为中心化AR(p)模型。
对于任意的AR序列模型，可以通过平移将其转为中心化的AR序列模型。构造方式是：
令$\mu=\frac{\varphi_0}{1-\sum\limits_{i=1}^p\varphi_i}$，则$x_t-\mu=\varphi(x_{t-1}-\mu)+...+\varphi_p(x_{t-p}-\mu)+\epsilon_t$
- 均值：$E(x_t)=\frac{\varphi_0}{1-\sum\limits_{i=1}^p\varphi_i}$
- （k阶）自相关系数（k阶自协方差除以方差）的通解形式：$\rho_k=\frac{\gamma_k}{\gamma_0}=\sum\limits^p_{i=1}c_i\lambda_i^k,\mid\lambda_i\mid<1\land c_i不恒为0$
- AR模型自相关系数的特点：
  - 指数衰减，证明AR序列的相关性往往出现在短期。
  - 拖尾性：不会在后期变为0，始终有非0的取值。


**MA模型**
即q阶移动平均（MovingAverage）模型，MA(q)的t时刻序列值受到一个白噪声序列过去q时间段内历史值的影响。
$$
\begin{cases}
    x_t=\mu+\epsilon_t+\theta_1\epsilon_{t-1}+...+\theta_q\epsilon_{t-q} \\
    \theta_q\ne 0 \\
    \{\epsilon_t\}是零均值白噪声序列 
\end{cases}
$$

$\mu=0$时称为中心化MA(q)模型。
- 均值是常数：$E(x_t)=E(\mu+\epsilon_t-\theta_1\epsilon_{t-1}-...-\theta_q\epsilon_{t-q})=\mu$
- 方差是常数：$D(x_t)=D(\mu+\epsilon_t-\theta_1\epsilon_{t-1}-...-\theta_q\epsilon_{t-q})=(1+\theta_1^2+...+\theta_q^2)\sigma^2_\epsilon$

**ARMA模型**

即自回归移动平均模型ARMA(p,q)
$$
\begin{cases}
  x_t=\phi_0+\phi_1x_{t-1}+...+\phi_px_{t-p}+\epsilon_t-\theta_1\epsilon_{t-1}-\theta_q\epsilon_{t-q} \\
  \{\epsilon_t\}是零均值白噪声序列 \\
  E(x_s\epsilon_t)=0,\forall s<t
\end{cases}
$$

$\phi_0=0$时称为中心化ARMA(p,q)模型。

**模型的平稳性**
即该模型生成的时间序列是平稳的，互为充要条件。

