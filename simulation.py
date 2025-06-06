import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import matplotlib

# 日本語フォント設定 (IPAexGothic または Noto Sans CJK JP)
matplotlib.rcParams['font.family'] = 'IPAexGothic'  
matplotlib.rcParams['axes.unicode_minus'] = False   # マイナス記号の文字化け防止

# --- Fortranシミュレーションと合わせたパラメータ ---
# Fortranコードで設定した値と一致させてください
N = 100           # 粒子数
v0 = 0.1          # 粒子の速さ
total_steps = 4000  # Fortranでの総ステップ数
save_interval = 10  # Fortranでのデータ保存間隔
# Fortranコードのetaの値をここに記述してください
current_eta = 0.2 

# --- データファイルの場所 ---
data_directory = 'traj'

def calculate_phi_from_file(filename, N_particles, v0_speed):
    """
    単一のデータファイルから秩序パラメータを計算する。
    [cite_start][cite: 17]
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
    
    print(f"{directory} ディレクトリのデータを解析します...")
    for i in range(1, num_saved_files + 1):
        filename = os.path.join(directory, f'pos{i:04d}.d')
        if os.path.exists(filename):
            phi = calculate_phi_from_file(filename, N_particles, v0_speed)
            if phi is not None:
                phi_vs_time.append(phi)
                time_steps.append(i * interval)

    if not time_steps:
        print(f"エラー: {directory} にデータファイルが見つかりませんでした。")
        print("Fortranシミュレーションを先に実行してください。")
        return None, None

    return np.array(time_steps), np.array(phi_vs_time)


# --- メインの描画処理 ---
if __name__ == '__main__':
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Vicsek Model Simulation Analysis', fontsize=16)

    # --- グラフ1: 秩序パラメータの時間変化 (P17 左図) ---
    ax1.set_title('Order Parameter vs. Time')
    ax1.set_xlabel('t (step)')
    ax1.set_ylabel('Order parameter $\\varphi(t)$')
    ax1.set_xscale('log') # X軸を対数スケールに [cite: 17]ax1.set_xscale('log') # X軸を対数スケールに [cite: 17]
    ax1.set_ylim(0, 1.0)
    ax1.grid(True, which="both", ls="--")

    # 現在のetaの値でシミュレーション結果を解析してプロット
    time_values, phi_values = analyze_simulation_run(data_directory, N, v0, total_steps, save_interval)
    
    if time_values is not None:
        ax1.plot(time_values, phi_values, label=f'$\\eta={current_eta}$')
        print(f"\neta={current_eta} の時間変化グラフを生成しました。")

    # (オプション) 別のeta値の結果もプロットしたい場合:
    # 1. Fortranコードのetaを例えば0.4に変更して再実行。
    # 2. 結果のtrajディレクトリをtraj_eta04などにリネーム。
    # 3. 以下のコードのコメントを解除し、ディレクトリ名を指定。
    # time_values_04, phi_values_04 = analyze_simulation_run('traj_eta04', N, v0, total_steps, save_interval)
    # if time_values_04 is not None:
    #    ax1.plot(time_values_04, phi_values_04, label='$\\eta=0.4$')

    ax1.legend()


    # --- グラフ2: 秩序パラメータ vs ノイズ強度 (P17 右図) ---
    ax2.set_title('Order Parameter vs. Noise Intensity')
    ax2.set_xlabel('Noise intensity $\\eta$')
    ax2.set_ylabel('Order parameter $\\langle\\varphi\\rangle$')
    ax2.set_xlim(0, 1.0)
    ax2.set_ylim(0, 1.0)
    ax2.grid(True, ls="--")

    # このグラフを生成するための注意書き
    ax2.text(0.5, 0.5, 
             'このグラフの生成手順:\n\n'
             '1. Fortranコードの`eta`の値を変更する。\n'
             '2. シミュレーションを再実行する。\n'
             '3. 出力された`traj`フォルダをリネームする\n'
             '   (例: `traj_eta0.0`, `traj_eta0.2`...)\n'
             '4. 各`eta`の結果から平衡状態の`phi`の平均値を\n'
             '   計算し、手動でプロットする。\n\n'
             '（以下のサンプルデータを参照）',
             ha='center', va='center', wrap=True, bbox=dict(boxstyle="round,pad=0.5", fc="wheat", alpha=0.5))

    # 手順に従って得られたであろうサンプルデータをプロット
    # これはスライド17の図を模したダミーデータです
    sample_eta =    [0.0,  0.1,  0.2,  0.3,  0.4,  0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sample_avg_phi =[0.99, 0.91, 0.68, 0.37, 0.15, 0.05, 0.02, 0.01, 0.01, 0.0, 0.0]
    ax2.plot(sample_eta, sample_avg_phi, 'o-', color='purple', label='Sample Data')
    ax2.legend()
    
    print("\nグラフ2は、複数のシミュレーション実行が必要なため、")
    print("サンプルデータと実行手順のガイドを表示しています。")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('vicsek_analysis_python.png')
    plt.show()
    print("\n解析グラフを 'vicsek_analysis_python.png' として保存しました。")


### 使い方

# 1.  **準備**: Fortranのシミュレーションを実行し、`traj` ディレクトリとデータファイルが生成されていることを確認します。
# 2.  **パラメータの確認**: Pythonコード冒頭の `current_eta` の値を、実行したFortranシミュレーションの `eta` の値と一致させます。
# 3.  **実行**: 上記のPythonコードをターミナルで実行します。
#     ```bash
#     python plot_analysis.py
#     ```
# 4.  **結果**:
#     * 画面にグラフが表示されます。左側のグラフには、`traj` ディレクトリのデータから計算された秩序パラメータの時間変化が表示されます。
#     * 右側のグラフには、スライド17ページのようなグラフを完成させるための手順と、サンプルデータに基づいたプロットが表示されます。
#     * `vicsek_analysis_python.png` という名前でグラフの画像ファイルが保存されます。

### グラフ2（秩序パラメータ vs ノイズ強度）を完成させるには

# 右側のグラフを実際のデータで埋めるには、以下の手作業が必要です。

# 1.  **Fortranコードの `eta` の値を変更**します（例: `eta = 0.0d0`）。
# 2.  **再コンパイルして実行**します。
# 3.  生成された `traj` ディレクトリを、どの `eta` の結果か分かるようにリネームします（例: `traj` -> `traj_eta0.0`）。
# 4.  `eta` の値を0.1, 0.2, ... , 1.0と変えながら、ステップ1〜3を繰り返します。
# 5.  全てのシミュレーションが終わったら、このPythonスクリプトを改造するか、あるいは各ディレクトリのデータから手動で平衡状態（例: 後半のステップ）の `phi` の平均値を計算し、`eta` の値に対してプロットしてくださ