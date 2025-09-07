import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib

# 设置中文字体支持
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

def sinh_function(i, a, b):
    """ a * sinh(b * i)"""
    return a * np.sinh(b * i)

def fit_sinh_function(x_data, y_data, initial_guess=None):
    """拟合双曲正弦函数"""
    if initial_guess is None:
        # 自动估计初始参数
        y_max = np.max(y_data)
        x_max = x_data[np.argmax(y_data)]
        a_guess = y_max / np.sinh(0.1 * x_max) if x_max > 0 else 1.0
        b_guess = 0.1
        initial_guess = [a_guess, b_guess]
    
    try:
        # 进行拟合
        popt, pcov = curve_fit(sinh_function, x_data, y_data, p0=initial_guess)
        
        # 计算拟合误差
        perr = np.sqrt(np.diag(pcov))
        
        # 计算R²
        y_pred = sinh_function(x_data, *popt)
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        return popt, perr, r_squared
    except Exception as e:
        print(f"拟合失败: {e}")
        return None, None, None

def plot_fitting_results(x_data, y_data, popt, perr, r_squared):
    """绘制拟合结果"""
    plt.figure(figsize=(10,6))
    
    # 原始数据点
    plt.scatter(x_data, y_data, color='red', alpha=0.6, s=50, label='original data', zorder=3)
    
    # 拟合曲线
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 200)
    y_fit = sinh_function(x_fit, *popt)
    plt.plot(x_fit, y_fit, 'b-', linewidth=2, label=f'拟合曲线: y = {popt[0]:.3f} × sinh({popt[1]:.3f} × i)', zorder=2)
    # 设置图形属性
    plt.xlabel('i', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title('sinh function fitting results', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # 添加参数信息文本框
    textstr = f'parameter:\n'
    textstr += f'a = {popt[0]:.4f} ± {perr[0]:.4f}\n'
    textstr += f'b = {popt[1]:.4f} ± {perr[1]:.4f}\n'
    textstr += f'R² = {r_squared:.4f}'
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    return plt.gcf()

def fit_and_plot(x_data, y_data, file_name='sinh_fitting'):
    """拟合并绘制数据"""
    popt, perr, r_squared = fit_sinh_function(x_data, y_data)
    if popt is not None:
        plot_fitting_results(x_data, y_data, popt, perr, r_squared)
        plt.savefig(f'{file_name}.png', dpi=300, bbox_inches='tight')
