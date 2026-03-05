import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def simulate_wafer():
    # 晶圆参数
    diameter = 300  # mm
    radius = diameter / 2
    N = 256
    
    # 创建网格
    x = np.linspace(-radius, radius, N)
    y = np.linspace(-radius, radius, N)
    X, Y = np.meshgrid(x, y)
    
    # 极坐标与掩膜
    rho = np.sqrt(X**2 + Y**2) / radius
    mask = rho <= 1.0
    
    # 生成基础形貌 (模拟离焦、倾斜和像散等低频误差)
    np.random.seed(100)
    surface = 15 * X/radius + 8 * (2*(rho**2)-1) + 12 * (X**2 - Y**2)/(radius**2)
    
    # 加入高频噪声
    noise = np.random.randn(N, N) * 5.0
    noise = gaussian_filter(noise, sigma=2)
    surface += noise
    
    # 应用掩膜，把晶圆外部的数据设为 NaN
    surface[~mask] = np.nan
    
    # 绘制结果
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, surface, cmap='plasma', shading='auto')
    plt.colorbar(label='Height (nm)')
    
    # 绘制边缘边界
    theta = np.linspace(0, 2*np.pi, 200)
    plt.plot(radius*np.cos(theta), radius*np.sin(theta), 'k-', linewidth=2)
    
    plt.title('Wafer Surface Simulation (test1.py)')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.axis('equal')
    plt.tight_layout()
    
    # 保存图片
    out_path = r'D:\ai_coding_python\test1_wafer_sim.png'
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f'Simulation successful. Output saved to: {out_path}')

if __name__ == '__main__':
    simulate_wafer()
