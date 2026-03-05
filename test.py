"""
硅片面型仿真程序
使用 Zernike 多项式模拟硅片面型误差，并提供多种可视化方式。
"""

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非交互式后端，服务器环境无需显示器
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# 设置中文字体（Windows 环境）
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False


# ────────────────────────────────────────────────
# Zernike 多项式（极坐标定义）
# ────────────────────────────────────────────────

def zernike_radial(n, m, rho):
    """计算 Zernike 径向多项式 R_n^m(rho)"""
    R = np.zeros_like(rho, dtype=float)
    for s in range((n - abs(m)) // 2 + 1):
        coef = ((-1) ** s * math.factorial(n - s) /
                (math.factorial(s) *
                 math.factorial((n + abs(m)) // 2 - s) *
                 math.factorial((n - abs(m)) // 2 - s)))
        R += coef * rho ** (n - 2 * s)
    return R


def zernike(n, m, rho, theta):
    """
    计算 Zernike 多项式 Z_n^m(rho, theta)
    n: 径向阶次, m: 角向频率
    """
    R = zernike_radial(n, abs(m), rho)
    if m > 0:
        return R * np.cos(m * theta)
    elif m < 0:
        return R * np.sin(-m * theta)
    else:
        return R


# ────────────────────────────────────────────────
# 硅片面型仿真
# ────────────────────────────────────────────────

class WaferSurface:
    """硅片面型仿真类"""

    def __init__(self, diameter_mm=300, grid_points=256):
        """
        参数:
            diameter_mm  : 硅片直径，单位 mm（默认 300mm，即 12 英寸）
            grid_points  : 网格点数
        """
        self.diameter = diameter_mm
        self.radius = diameter_mm / 2
        self.N = grid_points

        # 构建笛卡尔网格
        x = np.linspace(-self.radius, self.radius, self.N)
        y = np.linspace(-self.radius, self.radius, self.N)
        self.X, self.Y = np.meshgrid(x, y)

        # 极坐标（归一化半径）
        self.rho = np.sqrt(self.X**2 + self.Y**2) / self.radius
        self.theta = np.arctan2(self.Y, self.X)

        # 圆形孔径掩模
        self.mask = self.rho <= 1.0

        # 面型数据（nm）
        self.surface = np.full((self.N, self.N), np.nan)

    def generate(self, zernike_coeffs=None, noise_nm=2.0, seed=42):
        """
        生成仿真面型。

        参数:
            zernike_coeffs : dict，格式 {(n, m): 系数_nm}
                             若为 None 则使用预设的典型误差
            noise_nm       : 高频随机噪声幅值（nm）
            seed           : 随机种子
        """
        np.random.seed(seed)

        if zernike_coeffs is None:
            # 典型硅片面型误差：倾斜、离焦、像散、彗差、球差
            zernike_coeffs = {
                (1, -1): 5.0,    # Y 倾斜       5 nm
                (1,  1): 3.0,    # X 倾斜       3 nm
                (2,  0): 20.0,   # 离焦         20 nm
                (2, -2): 8.0,    # 45° 像散     8 nm
                (2,  2): -6.0,   # 0° 像散     -6 nm
                (3, -1): 4.0,    # Y 彗差       4 nm
                (3,  1): -3.0,   # X 彗差      -3 nm
                (4,  0): 10.0,   # 球差         10 nm
                (4, -4): 2.0,    # 四叶像差     2 nm
            }

        surface = np.zeros((self.N, self.N))

        for (n, m), coef in zernike_coeffs.items():
            Z = zernike(n, m, self.rho, self.theta)
            surface += coef * Z

        # 叠加高频随机噪声
        noise = noise_nm * np.random.randn(self.N, self.N)
        # 对噪声做简单平滑，模拟中频误差
        from scipy.ndimage import gaussian_filter
        noise = gaussian_filter(noise, sigma=3)
        surface += noise

        # 仅保留圆形孔径内的数据
        self.surface = np.where(self.mask, surface, np.nan)
        return self

    def pv(self):
        """Peak-to-Valley（nm）"""
        valid = self.surface[self.mask]
        return valid.max() - valid.min()

    def rms(self):
        """均方根误差（nm）"""
        valid = self.surface[self.mask]
        return np.sqrt(np.mean((valid - valid.mean())**2))

    def stats(self):
        print("=" * 40)
        print(f"  硅片直径   : {self.diameter} mm")
        print(f"  PV 面型误差: {self.pv():.2f} nm")
        print(f"  RMS 面型误差: {self.rms():.2f} nm")
        print("=" * 40)


# ────────────────────────────────────────────────
# 可视化函数
# ────────────────────────────────────────────────

def plot_all(wafer: WaferSurface):
    """一次性绘制 2D 色图、3D 曲面、X/Y 截面四幅图"""

    surface = wafer.surface.copy()
    X, Y = wafer.X, wafer.Y
    N = wafer.N
    r = wafer.radius

    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(f'硅片面型仿真  |  直径 {wafer.diameter} mm  |  '
                 f'PV={wafer.pv():.1f} nm  RMS={wafer.rms():.1f} nm',
                 fontsize=14, fontweight='bold')

    # ── 1. 2D 伪彩色图 ──────────────────────────────
    ax1 = fig.add_subplot(2, 2, 1)
    im = ax1.pcolormesh(X, Y, surface, cmap='RdYlBu_r', shading='auto')
    fig.colorbar(im, ax=ax1, label='面型高度 (nm)')
    # 绘制晶圆边缘圆
    theta_c = np.linspace(0, 2 * np.pi, 360)
    ax1.plot(r * np.cos(theta_c), r * np.sin(theta_c), 'k-', lw=1.5)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax1.set_title('面型误差 2D 分布图')

    # ── 2. 3D 曲面图 ────────────────────────────────
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    # 用 masked 数组让孔径外透明
    surf_masked = np.ma.array(surface, mask=~wafer.mask)
    ax2.plot_surface(X, Y, surf_masked, cmap='RdYlBu_r',
                     rstride=4, cstride=4, alpha=0.9, linewidth=0)
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    ax2.set_zlabel('高度 (nm)')
    ax2.set_title('面型误差 3D 曲面图')
    ax2.view_init(elev=30, azim=-60)

    # ── 3. X 方向截面 ───────────────────────────────
    ax3 = fig.add_subplot(2, 2, 3)
    mid = N // 2
    x_line = X[mid, :]
    y_line = surface[mid, :]
    valid = ~np.isnan(y_line)
    ax3.plot(x_line[valid], y_line[valid], 'b-', lw=1.5, label='Y=0 截面')
    ax3.axhline(0, color='k', lw=0.8, ls='--')
    ax3.fill_between(x_line[valid], y_line[valid], alpha=0.15)
    ax3.set_xlabel('X (mm)')
    ax3.set_ylabel('高度 (nm)')
    ax3.set_title('X 方向截面轮廓（Y=0）')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # ── 4. Y 方向截面 ───────────────────────────────
    ax4 = fig.add_subplot(2, 2, 4)
    y_col = Y[:, mid]
    x_col = surface[:, mid]
    valid2 = ~np.isnan(x_col)
    ax4.plot(y_col[valid2], x_col[valid2], 'r-', lw=1.5, label='X=0 截面')
    ax4.axhline(0, color='k', lw=0.8, ls='--')
    ax4.fill_between(y_col[valid2], x_col[valid2], alpha=0.15, color='r')
    ax4.set_xlabel('Y (mm)')
    ax4.set_ylabel('高度 (nm)')
    ax4.set_title('Y 方向截面轮廓（X=0）')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(r'D:\ai_coding_python\wafer_surface.png', dpi=150, bbox_inches='tight')
    print("图像已保存至 D:/ai_coding_python/wafer_surface.png")
    plt.close()


def plot_zernike_decomposition(wafer: WaferSurface, coeffs: dict):
    """绘制各 Zernike 项的单独贡献"""
    items = list(coeffs.items())
    n_items = len(items)
    cols = 3
    rows = (n_items + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(14, rows * 4))
    axes = axes.flatten()

    zernike_names = {
        (1, -1): 'Y 倾斜', (1, 1): 'X 倾斜',
        (2, 0): '离焦', (2, -2): '45°像散', (2, 2): '0°像散',
        (3, -1): 'Y 彗差', (3, 1): 'X 彗差',
        (4, 0): '球差', (4, -4): '四叶像差',
    }

    for idx, ((n, m), coef) in enumerate(items):
        Z = zernike(n, m, wafer.rho, wafer.theta)
        component = np.where(wafer.mask, coef * Z, np.nan)

        im = axes[idx].pcolormesh(wafer.X, wafer.Y, component,
                                   cmap='RdYlBu_r', shading='auto')
        fig.colorbar(im, ax=axes[idx], label='nm')
        theta_c = np.linspace(0, 2 * np.pi, 360)
        axes[idx].plot(wafer.radius * np.cos(theta_c),
                       wafer.radius * np.sin(theta_c), 'k-', lw=1)
        axes[idx].set_aspect('equal')
        name = zernike_names.get((n, m), f'Z({n},{m})')
        axes[idx].set_title(f'{name}\n系数={coef} nm')
        axes[idx].set_xlabel('X (mm)')
        axes[idx].set_ylabel('Y (mm)')

    for idx in range(n_items, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Zernike 各项面型分量', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(r'D:\ai_coding_python\wafer_zernike.png', dpi=150, bbox_inches='tight')
    print("Zernike 分解图已保存至 D:/ai_coding_python/wafer_zernike.png")
    plt.close()


# ────────────────────────────────────────────────
# 主程序
# ────────────────────────────────────────────────

if __name__ == '__main__':
    # 自定义 Zernike 系数（单位：nm），可根据需要修改
    coeffs = {
        (1, -1): 5.0,    # Y 倾斜
        (1,  1): 3.0,    # X 倾斜
        (2,  0): 20.0,   # 离焦
        (2, -2): 8.0,    # 45° 像散
        (2,  2): -6.0,   # 0° 像散
        (3, -1): 4.0,    # Y 彗差
        (3,  1): -3.0,   # X 彗差
        (4,  0): 10.0,   # 球差
        (4, -4): 2.0,    # 四叶像差
    }

    # 创建 300mm 硅片，256×256 网格
    wafer = WaferSurface(diameter_mm=300, grid_points=256)
    wafer.generate(zernike_coeffs=coeffs, noise_nm=2.0, seed=42)

    # 打印统计信息
    wafer.stats()

    # 绘制综合图（2D + 3D + 截面）
    plot_all(wafer)

    # 绘制 Zernike 各项分解图
    plot_zernike_decomposition(wafer, coeffs)
