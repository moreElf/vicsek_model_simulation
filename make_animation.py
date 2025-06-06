import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob
from tqdm import tqdm

# --- 設定項目 ---
# アニメーション化したいデータディレクトリを指定
DATA_DIRECTORY = 'traj2'
# Fortranシミュレーションのパラメータ
L = 10.0  # ボックスサイズ
N = 100   # 粒子数
v0 = 0.1  # 速度

# --- データの読み込み ---
file_pattern = os.path.join(DATA_DIRECTORY, 'pos*.d')
file_list = sorted(glob.glob(file_pattern))

if not file_list:
    print(f"エラー: ディレクトリ '{DATA_DIRECTORY}' にデータファイルが見つかりません。")
    print("Fortranシミュレーションを先に実行してください。")
    exit()

all_data = []
print(f"'{DATA_DIRECTORY}' からデータを読み込んでいます...")
for filename in tqdm(file_list, desc="Loading data"):
    try:
        data = np.loadtxt(filename)
        if data.shape[0] == N and data.shape[1] == 4: # 列数もチェック
            all_data.append(data)
        else:
            print(f"Warning: ファイル {filename} のデータ形式が不正です (期待する形状: ({N}, 4))。スキップします。")
    except Exception as e:
        print(f"Warning: ファイル {filename} の読み込みに失敗しました: {e}。スキップします。")

if not all_data:
    print(f"エラー: ディレクトリ '{DATA_DIRECTORY}' から読み込める有効なデータフレームがありませんでした。")
    exit()

all_data = np.array(all_data)
num_frames = len(all_data)

if num_frames == 0:
    print(f"エラー: ディレクトリ '{DATA_DIRECTORY}' から読み込める有効なデータフレームがありませんでした。")
    exit()

# --- アニメーションのセットアップ ---
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')

# 最初のフレームのデータを取得してプロットオブジェクトを初期化
first_frame_data = all_data[0]
initial_positions_x = first_frame_data[:, 0]
initial_positions_y = first_frame_data[:, 1]
initial_velocities_x = first_frame_data[:, 2]
initial_velocities_y = first_frame_data[:, 3]

# 粒子（点）の初期プロット
particles_plot, = ax.plot(initial_positions_x, initial_positions_y, 'o', color='teal', markersize=5)

# ベクトル（矢印）の初期プロット
# ★★★ 修正点: quiverを最初のフレームのデータで初期化 ★★★
quiver_plot = ax.quiver(initial_positions_x, initial_positions_y,
                           initial_velocities_x, initial_velocities_y,
                           scale=1/v0, # スケールは1/v0のままか、あるいは調整が必要かもしれません
                           angles='xy', scale_units='xy', # units='xy'やangles='xy'も挙動に影響します
                           color='teal', alpha=0.8, width=0.005)

time_text = ax.set_title(f'Frame 0 / {num_frames -1}') # 初期フレームのタイトル

# --- アニメーションの更新関数 ---
def update(frame):
    data = all_data[frame]
    positions = data[:, 0:2]
    velocities = data[:, 2:4]

    particles_plot.set_data(positions[:, 0], positions[:, 1])
    
    # quiverの更新
    quiver_plot.set_offsets(positions)
    quiver_plot.set_UVC(velocities[:, 0], velocities[:, 1])

    time_text.set_text(f'Frame {frame} / {num_frames - 1}')
    return particles_plot, quiver_plot, time_text

# --- アニメーションの初期化関数 (blit=Trueの場合に推奨) ---
def init():
    # 最初のフレームの状態を再度設定（または初期状態に戻す）
    particles_plot.set_data(initial_positions_x, initial_positions_y)
    quiver_plot.set_offsets(np.column_stack([initial_positions_x, initial_positions_y]))
    quiver_plot.set_UVC(initial_velocities_x, initial_velocities_y)
    time_text.set_text(f'Frame 0 / {num_frames - 1}')
    return particles_plot, quiver_plot, time_text


# --- アニメーションの作成と保存 ---
# ★★★ 修正点: init_func を追加 ★★★
ani = animation.FuncAnimation(fig, update, frames=num_frames,
                              init_func=init, blit=True, interval=50)

output_filename_mp4 = "vicsek_animation_2.mp4"
output_filename_gif = "vicsek_animation_2.gif"

try:
    print(f"\nアニメーションを '{output_filename_mp4}' として保存しています... (数分かかる場合があります)")
    ani.save(output_filename_mp4, writer='ffmpeg', fps=20, dpi=100)
    print(f"'{output_filename_mp4}' の保存が完了しました。")
except Exception as e:
    print(f"\nmp4の保存に失敗しました: {e}")
    print("ffmpegがインストールされていない可能性があります。")
    try:
        print(f"代わりに '{output_filename_gif}' として保存を試みます...")
        ani.save(output_filename_gif, writer='imagemagick', fps=20)
        print(f"'{output_filename_gif}' の保存が完了しました。")
    except Exception as e2:
        print(f"gifの保存にも失敗しました: {e2}")
        print("ImageMagickがインストールされていない可能性があります。")
        print("アニメーションの保存ができませんでした。")
