
# 运行方法

运行 main.m 文件即可运行程序。

# 算例

- tb = $0.5/\pi$
- $x\in [0,2\pi]$
- $u_t + (\frac{1}{2}u^2)_x=0$
- $u0 = \sin(x)$



周期边界条件

# 参数设置
- 基函数的次数为 2
- dt = 0.01h
- 不使用限制器
- 数值通量选择 Local Lax-Friedrichs
- 时间离散为三次龙格库塔法

# 误差阶

| n   | 误差1范数 | 误差阶   | 误差无穷范数 | 误差阶   |
|-----|-----------|----------|--------------|----------|
| 8   | 0.0036    |          | 6.84E-04     |          |
| 16  | 0.0018    | 1        | 1.93E-04     | 1.83 |
| 32  | 5.91E-04  | 1.606913 | 3.53E-05     | 2.45 |
| 64  | 1.59E-04  | 1.898529 | 5.59E-06     | 2.66 |
| 128 | 4.16E-05  | 1.928129 | 7.89E-07     | 2.83 |
| 256 | 1.07E-05  | 1.955691 | 1.06E-07     | 2.90 |
| 512 | 2.73E-06  | 1.9734   | 1.40E-08     | 2.92 |

# 参考文献 
[1] Cockburn, Bernardo, and Chi-Wang Shu. "TVB Runge-Kutta local projection discontinuous Galerkin finite element method for conservation laws. II. General framework." Mathematics of computation 52.186 (1989): 411-435.