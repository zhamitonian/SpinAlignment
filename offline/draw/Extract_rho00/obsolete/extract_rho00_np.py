import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# 定义角分布拟合函数（带幅度参数）
def angular_distribution(cos_theta, A, rho00):
    """
    W(theta) = A * 0.75 * (1 - rho00 + (3*rho00 - 1) * cos^2(theta))
    A: 幅度参数（归一化因子）
    rho00: 自旋密度矩阵元
    """
    return A * 0.75 * (1 - rho00 + (3 * rho00 - 1) * cos_theta**2)

# 读取数据
df = pd.read_csv("./images/nsig_results.txt", sep=r"\s+",
                 names=["xp_center", "xp_width", "helicity_center", 
                        "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                 skiprows=1)

# 过滤掉全零的无效数据
df_valid = df[(df['nsig'] > 0) | (df['xp_center'] > 0)]

# 按 xp_center 分组
grouped = df_valid.groupby('xp_center')

# 获取唯一的 xp_center 值
unique_xp = sorted(df_valid['xp_center'].unique())
# 过滤掉 0.0 值（无效数据）
unique_xp = [x for x in unique_xp if x > 0]

# 计算子图布局
n_groups = len(unique_xp)
n_cols = 4  # 每行4个图
n_rows = int(np.ceil(n_groups / n_cols))

# 创建子图
fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
axes = axes.flatten() if n_groups > 1 else [axes]

# 存储拟合结果
xp_values = []
rho00_values = []
rho00_errors = []

# 对每个 xp_center 值绘图和拟合
for idx, xp_value in enumerate(unique_xp):
    ax = axes[idx]
    
    # 获取该组数据
    group = grouped.get_group(xp_value)
    
    # 提取数据
    cos_theta = group['helicity_center'].values  # 已经是 cos(theta) 值
    nsig = group['nsig'].values
    nsig_err = group['nsig_err'].values
    
    # 执行拟合（带幅度参数）
    try:
        # 初始猜测：A 为数据平均值，rho00 为 0.3
        A_init = np.mean(nsig)
        # 使用加权最小二乘拟合（考虑误差）
        popt, pcov = curve_fit(angular_distribution, cos_theta, nsig, 
                               sigma=nsig_err, 
                               p0=[A_init, 0.3], 
                               bounds=([0, -1], [np.inf, 1]))  # A > 0, rho00 ∈ [-1, 1]
        A_fit = popt[0]
        rho00_fit = popt[1]
        rho00_err = np.sqrt(pcov[1, 1])
        
        # 保存结果
        xp_values.append(xp_value)
        rho00_values.append(rho00_fit)
        rho00_errors.append(rho00_err)
        
        # 生成拟合曲线
        cos_theta_fine = np.linspace(-1, 1, 100)
        nsig_fit = angular_distribution(cos_theta_fine, A_fit, rho00_fit)
        
        fit_success = True
        fit_text = f'$\\rho_{{00}}$ = {rho00_fit:.3f} ± {rho00_err:.3f}'
    except Exception as e:
        print(f"拟合失败 (xp={xp_value}): {e}")
        fit_success = False
        fit_text = 'Fit failed'
    
    # 绘制带误差棒的数据点
    ax.errorbar(cos_theta, nsig, yerr=nsig_err, 
                fmt='o', markersize=5, capsize=3, capthick=1,
                label='Data', color='blue', alpha=0.7)
    
    # 绘制拟合曲线
    if fit_success:
        ax.plot(cos_theta_fine, nsig_fit, 'r-', linewidth=2, 
                label=fit_text)
    
    # 设置标题和标签
    ax.set_title(f'$x_p$ = {xp_value:.3f}', fontsize=12)
    ax.set_xlabel('$\\cos(\\theta_{helicity})$', fontsize=10)
    ax.set_ylabel('$N_{sig}$', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='best')

# 隐藏多余的子图
for idx in range(n_groups, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.savefig('nsig_vs_helicity_fits.png', dpi=150, bbox_inches='tight')
print(f"拟合图已保存为 nsig_vs_helicity_fits.png")
plt.show()

# ===== 绘制 rho00 vs xp 图 =====
if len(xp_values) > 0:
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    # 绘制 rho00 vs xp 带误差棒
    ax2.errorbar(xp_values, rho00_values, yerr=rho00_errors,
                fmt='o-', markersize=8, capsize=5, capthick=2,
                linewidth=2, color='darkblue', ecolor='red',
                label='$\\rho_{00}$ from fits')
    
    # 添加参考线 rho00 = 1/3 (无极化)
    ax2.axhline(y=1/3, color='gray', linestyle='--', linewidth=1.5,
                label='$\\rho_{00} = 1/3$ (unpolarized)')
    
    # 设置标签和标题
    ax2.set_xlabel('$x_p$ (Scaled Momentum)', fontsize=14)
    ax2.set_ylabel('$\\rho_{00}$ (Spin Density Matrix)', fontsize=14)
    ax2.set_title('Spin Alignment Parameter $\\rho_{00}$ vs $x_p$', fontsize=16)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=12)
    ax2.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig('rho00_vs_xp.png', dpi=150, bbox_inches='tight')
    print(f"rho00 结果图已保存为 rho00_vs_xp.png")
    plt.show()
    
    # 保存拟合结果到文本文件
    with open('rho00_results.txt', 'w') as f:
        f.write('# xp_center  rho00  rho00_error\n')
        for xp, rho00, rho00_err in zip(xp_values, rho00_values, rho00_errors):
            f.write(f'{xp:.4f}  {rho00:.6f}  {rho00_err:.6f}\n')
    print(f"拟合结果已保存到 rho00_results.txt")
    
    # 打印结果摘要
    print("\n===== 拟合结果摘要 =====")
    print(f"{'xp':>8s}  {'rho00':>10s}  {'error':>10s}")
    print("-" * 32)
    for xp, rho00, rho00_err in zip(xp_values, rho00_values, rho00_errors):
        print(f'{xp:8.4f}  {rho00:10.6f}  {rho00_err:10.6f}')
else:
    print("警告: 没有成功的拟合结果！")
