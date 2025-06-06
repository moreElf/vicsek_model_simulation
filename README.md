# Vicsekモデル シミュレーション

自己駆動粒子の集団運動を記述するVicsekモデルのシミュレーションです。
このリポジトリには、Fortranによるシミュレーションコードと、Pythonによる結果の解析・可視化コードが含まれています。

## シミュレーションパラメータ (Simulation Parameters)

本シミュレーションは、以下のパラメータ設定に基づいて行われました。

### 固定パラメータ

以下のパラメータは、すべてのシミュレーションで共通の値に設定されています。

| パラメータ | 記号 | 値 |
| :--- | :---: | :--- |
| 粒子数 | `N` | `100` |
| ボックスサイズ | `L` | `10.0` |
| 相互作用半径 | `R₀` | `1.0` |
| 粒子の速さ | `v₀` | `0.1` |
| 時間ステップ | `dt` | `1.0` |

### 変更したパラメータ

本シミュレーションでは、系の秩序形成に与える影響を調査するため、**ノイズ強度 (η)** の値を変化させました。実行したシミュレーションと、対応するパラメータおよび結果の格納ディレクトリは以下の通りです。

| 相当ファイル名 | ノイズ強度 (η) | 結果格納ディレクトリ |
| :--- | :---: | :--- |
| vicsek_animation.mp4/ vicsek_animation.gif | `0.2` | `traj/` | 
| vicsek_animation2.mp4/ vicsek_animation2.gif | `0.4` | `traj2/` |

これらの異なるノイズ強度におけるシミュレーション結果を比較することで、ノイズが群れの秩序形成に与える影響を考察します。

## 実行方法

1. **コンパイル**:

ηの値を変化させたいため、二回のコンパイルを行うが、その際に、vicsek_simulation.f90の10行目でetaの値を変更し、55行目のコンパイル先のフォルダ名を変更してそれぞれ違うコンパイルをする。

# グラフの生成
```bash
python3 plot_comparison.py
```

# アニメーションの生成
```bash
python3 make_animation.py
```
