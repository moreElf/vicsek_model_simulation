import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# --- Fortranシミュレーションと合わせたパラメータ ---
# Fortranコードで設定した値と一致させてください
N = 100           # 粒子数
v0 = 0.1          # 粒子の速さ
total_steps = 4000  # Fortranでの総ステップ数
save_interval = 10  # Fortranでのデータ保存間隔

# --- 比較したい2つのシミュレーションの情報 ---
# ディレクトリ名と、そのシミュレーションで使ったetaの値を設定
run1_info = {'directory': 'traj', 'eta': 0.2}
run2_info = {'directory': 'traj2', 'eta': 0.4}


def calculate_phi_from_file(filename, N_particles, v0_speed):
    """
    単一のデータファイルから秩序パラメータを計算する。
    """
    try:
        # ファイルからvx, vyデータを読み込む (3列目と4列目)
        data = np.loadtxt(filename, usecols=(2, 3))
        if data.shape[0] != N_particles:
            print(f"Warning: {filename} の粒子数が {data.shape[0]} であり、期待値 {N_particles} と異なります。")
            return None
            
        vx = data[:, 0]
        vy = data[:, 1]
    except (IOError, IndexError) as e:
        print(f"Error reading or parsing file {filename}: {e}")
        return None

    # 秩序パラメータを計算
    # phi = (1 / (N*v0)) * |sum(v)|
    sum_vx = np.sum(vx)
    sum_vy = np.sum(vy)
    phi = (1 / (N_particles * v0_speed)) * np.sqrt(sum_vx**2 + sum_vy**2)
    return phi

def analyze_simulation_run(directory, N_particles, v0_speed, steps, interval):
    """
    指定されたディレクトリの全データファイルから秩序パラメータの時系列を計算する。
    """
    phi_vs_time = []
    time_steps = []

    # 保存されたステップ数を計算
    num_saved_files = steps // interval
    
    print(f"'{directory}' ディレクトリのデータを解析します...")
    for i in range(1, num_saved_files + 1):
        filename = os.path.join(directory, f'pos{i:04d}.d')
        if os.path.exists(filename):
            phi = calculate_phi_from_file(filename, N_particles, v0_speed)
            if phi is not None:
                phi_vs_time.append(phi)
                time_steps.append(i * interval)

    if not time_steps:
        print(f"エラー: '{directory}' にデータファイルが見つかりませんでした。")
        return None, None

    return np.array(time_steps), np.array(phi_vs_time)


# --- メインの描画処理 ---
if __name__ == '__main__':
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Vicsek Model Simulation Comparison', fontsize=16)

    # --- グラフ1: 秩序パラメータの時間変化の比較 (P17 左図) ---
    ax1.set_title('Order Parameter vs. Time')
    ax1.set_xlabel('t (step)')
    ax1.set_ylabel('Order parameter $\\varphi(t)$')
    ax1.set_xscale('log')
    ax1.set_ylim(0, 1.0)
    ax1.grid(True, which="both", ls="--")

    time1, phi1 = analyze_simulation_run(run1_info['directory'], N, v0, total_steps, save_interval)
    if time1 is not None:
        ax1.plot(time1, phi1, label=f"$\\eta={run1_info['eta']}$")

    time2, phi2 = analyze_simulation_run(run2_info['directory'], N, v0, total_steps, save_interval)
    if time2 is not None:
        ax1.plot(time2, phi2, label=f"$\\eta={run2_info['eta']}$")

    ax1.legend()


    # --- グラフ2: 秩序パラメータ vs ノイズ強度 (P17 右図) ---
    ax2.set_title('Order Parameter vs. Noise Intensity')
    ax2.set_xlabel('Noise intensity $\\eta$')
    ax2.set_ylabel('Order parameter $\\langle\\varphi\\rangle$')
    ax2.set_xlim(0, 1.0)
    ax2.set_ylim(0, 1.0)
    ax2.grid(True, ls="--")

    # ax2.text(0.5, 0.5, 
    #          'このグラフを完成させるには、\n'
    #          '様々なetaでシミュレーションを繰り返し、\n'
    #          '各結果の平衡状態での平均値を\n'
    #          'プロットする必要があります。\n\n'
    #          '（以下のサンプルデータを参照）',
    #          ha='center', va='center', wrap=True, bbox=dict(boxstyle="round,pad=0.5", fc="wheat", alpha=0.5))

    sample_eta =    [0.0,  0.1,  0.2,  0.3,  0.4,  0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sample_avg_phi =[0.99, 0.91, 0.68, 0.37, 0.15, 0.05, 0.02, 0.01, 0.01, 0.0, 0.0]
    
    # ★★★ 修正点: 'o-' を 'o' に変更し、線なしで点のみプロット ★★★
    ax2.plot(sample_eta, sample_avg_phi, 'o', color='purple', label='Sample Data')
    ax2.legend()
    

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('vicsek_comparison.png')
    plt.show()
    print("\n比較グラフを 'vicsek_comparison.png' として保存しました。")
