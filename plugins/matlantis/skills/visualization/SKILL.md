---
name: mt-visualization
description: >
  原子構造とトラジェクトリの可視化を扱うスキルです。
  show_gui, pfcc_extras, nglview, show_ase, show_asetraj, view_ngl,
  ase.visualize.view, POV-Ray, 可視化, visualization,
  構造確認, 目視チェック, trajectory 表示
  に関するコード生成時に使用してください。
---

# 可視化

## 概要

Matlantis の JupyterLab 環境における原子構造とトラジェクトリの可視化ガイドです。シミュレーション結果の検証や構造確認は、計算ワークフロー全体の品質を左右する重要なステップです。

Matlantis 環境では主に以下の可視化手段を利用できます。

- **pfcc_extras.show_gui**: Matlantis 標準の可視化ツール。単一構造とトラジェクトリの両方に対応
- **nglview**: 高機能な分子可視化ライブラリ。インタラクティブな操作と GUI パネルを提供
- **ASE visualization**: ASE 組み込みの可視化機能。nglview バックエンドを指定可能
- **POV-Ray**: 高品質な静止画レンダリング用

可視化は計算の補助的な役割ですが、構造の妥当性確認を怠ると「原子が重なっている」「PBC が正しく設定されていない」といった初歩的なミスに気づかず、長時間の計算を無駄にすることがあります。計算投入前の可視化チェックは必須の手順です。

本ガイドでは、可視化の目的に応じた手法の選択基準と、各手法の具体的な使い方を説明します。

## ワークフロー

可視化の標準的なワークフローは以下の通りです。

1. **目的の明確化**: 何を確認したいかを決める（構造の妥当性、最適化の前後比較、トラジェクトリの追跡など）
2. **手法の選択**: 目的に応じた可視化手法を選ぶ
3. **構造・トラジェクトリの表示**: 選択した手法で表示する
4. **視覚的検証**: 結合、周期境界、原子の重なり等をチェックする
5. **必要に応じた調整**: フレーム間引き、表示オプションの変更など

### 可視化手法の判断チェックリスト

以下のチェックリストに基づいて、適切な可視化手法を選択してください。

| 確認事項 | 推奨手法 |
| :--- | :--- |
| 構造をざっと確認したい | `pfcc_extras.show_gui` |
| 最適化前後の比較が必要 | `nglview` で 2 構造を並べて表示 |
| トラジェクトリを再生したい | `nglview.show_asetraj()` またはフレーム間引き表示 |
| 距離や角度の測定が必要 | `nglview` の GUI パネル（gui=True） |
| 高品質な静止画が必要 | POV-Ray レンダリング |
| 大規模系（1000+ 原子）の表示 | フレーム間引き + 軽量表示設定 |
| static image として保存したい | POV-Ray または matplotlib |
| interactive viewer が必要 | `nglview` |

## 実装パターン

### パターン A: pfcc_extras.show_gui による可視化

Matlantis 環境で最も手軽に使える可視化手法です。単一の Atoms オブジェクトとトラジェクトリ（List[Atoms]）の両方に対応しています。Matlantis が提供する `pfcc_extras` パッケージに含まれています。

```python
from pfcc_extras import show_gui

# 単一構造の可視化
show_gui(atoms)

# トラジェクトリの可視化
show_gui(trajectory)  # List[Atoms] を渡す
```

#### 用途

- 計算投入前の構造確認（結合の妥当性、PBC のチェック、原子の重なり検出）
- 構造最適化後の結果確認
- MD トラジェクトリの簡易再生

#### 入力形式

`show_gui` は以下の形式を受け付けます。

- `ase.Atoms`: 単一構造
- `List[ase.Atoms]`: トラジェクトリ（最適化過程、MD 過程など）
- ASE Trajectory オブジェクト

#### 注意事項

- JupyterLab 環境でのみ動作します
- 大規模系では表示が重くなる場合があります
- Notebook セル内で呼び出す必要があります（スクリプト実行では表示されません）

### パターン B: nglview による単一構造の可視化

nglview は高機能な分子可視化ライブラリで、インタラクティブな操作が可能です。回転、ズーム、パン操作に加え、原子の選択や距離・角度の測定ができます。

```python
import nglview as nv

# 単一構造の表示（GUI パネル付き）
view = nv.show_ase(atoms, gui=True)
view
```

`gui=True` を指定すると、距離・角度測定やレプリゼンテーション変更などの GUI パネルが表示されます。

#### 表示オプションの調整

```python
view = nv.show_ase(atoms, gui=True)

# 表示レプリゼンテーションの変更
view.clear_representations()
view.add_ball_and_stick()

# セルの表示（結晶系で有用）
view.add_unitcell()

# ズームの調整
view.center()
```

#### 利用可能なレプリゼンテーション

nglview では以下のレプリゼンテーションを利用できます。

| メソッド | 説明 | 推奨用途 |
| :--- | :--- | :--- |
| `add_ball_and_stick()` | ボール・アンド・スティック表示 | 分子系、小規模結晶 |
| `add_spacefill()` | 空間充填モデル | 大規模系（軽量） |
| `add_licorice()` | ボンド強調表示 | 有機分子の結合確認 |
| `add_surface()` | 表面表示 | タンパク質等の大規模分子 |

### パターン C: nglview によるトラジェクトリの可視化

MD や構造最適化のトラジェクトリをアニメーションとして再生できます。

```python
import nglview as nv

# トラジェクトリの表示
view = nv.show_asetraj(trajectory, gui=True)
view
```

#### フレーム間引き

大規模なトラジェクトリではフレームを間引いて表示することで、メモリ使用量と表示速度を改善できます。

```python
# 10フレームごとに間引いて表示
view = nv.show_asetraj(trajectory[::10], gui=True)
view
```

#### 間引き間隔の選び方

| トラジェクトリの規模 | 推奨間引き間隔 | 表示フレーム数の目安 |
| :--- | :--- | :--- |
| 構造最適化（~100 steps） | 間引き不要 | ~100 |
| MD（~1,000 steps） | `[::5]` | ~200 |
| MD（~10,000 steps） | `[::50]` | ~200 |
| MD（~100,000 steps） | `[::500]` | ~200 |

目安として、表示フレーム数が 200 以下になるように間引き間隔を設定してください。

### パターン D: ASE visualization（ngl バックエンド）

> **注意**: `ase.visualize.view` よりも `pfcc_extras.show_gui` または `nglview` の直接使用を推奨します。`show_gui` はより簡潔で、Matlantis 環境に最適化されています。

ASE 組み込みの `view` 関数に nglview バックエンドを指定する方法です。

```python
from ase.visualize import view

# nglview バックエンドで表示
view(atoms, viewer="ngl")
```

この方法は ASE のインターフェースに慣れている場合に便利ですが、nglview を直接使う場合と比べて細かな制御は制限されます。レプリゼンテーションの変更や GUI パネルの表示などが必要な場合は、`nv.show_ase()` を直接使用してください。

### パターン E: POV-Ray による高品質レンダリング

論文や発表用の高品質な静止画を生成する場合に使用します。ラスタライズベースのレンダリングにより、光源、影、反射を含む写実的な画像を出力できます。

```python
from ase.io import write

# POV-Ray 入力ファイルの生成
write(
    'structure.pov',
    atoms,
    rotation='10x,20y,30z',   # 視点の回転角度
    show_unit_cell=2,          # ユニットセルの表示方法
    colors=None,               # デフォルトのJmol色を使用
)

# POV-Ray でレンダリング（要 povray インストール）
# !povray structure.pov
```

POV-Ray は Matlantis 環境にデフォルトでインストールされていない場合があります。必要に応じてインストールしてください。

#### POV-Ray のカスタマイズオプション

```python
write(
    'structure.pov',
    atoms,
    rotation='0x,0y,0z',
    show_unit_cell=2,
    radii=0.85,                # 原子の半径スケール
    bond_type='default',       # 結合の表示タイプ
    celllinewidth=0.05,        # セル境界線の太さ
)
```

### パターン F: 最適化前後の比較

構造最適化の前後で構造を比較する場合の手法です。

```python
import nglview as nv

# 最適化前
view_before = nv.show_ase(atoms_before, gui=True)
print("Before optimization:")
display(view_before)

# 最適化後
view_after = nv.show_ase(atoms_after, gui=True)
print("After optimization:")
display(view_after)
```

または、トラジェクトリとして最適化過程全体を表示します。

```python
# 最適化トラジェクトリの表示
from ase.io import read
opt_traj = read("optimization.traj", index=":")
view = nv.show_asetraj(opt_traj, gui=True)
view
```

### パターン G: 大規模系の可視化

1000 原子以上の大規模系では、以下の点に注意してください。

```python
import nglview as nv

# 大規模系ではレプリゼンテーションを軽量なものに変更
view = nv.show_ase(large_atoms)
view.clear_representations()
view.add_spacefill(radius_type='vdw', scale=0.3)  # 軽量な表示
view
```

大規模系の可視化で考慮すべき点:

- **メモリ**: 原子数が多いと Jupyter カーネルのメモリを大量に消費します
- **描画速度**: Ball-and-stick 表示よりも Spacefill 表示の方が軽量です
- **トラジェクトリ**: フレーム間引きを積極的に行ってください
- **部分表示**: 関心のある領域だけを切り出して表示することも検討してください

#### 部分表示の例

```python
# 特定の元素のみ表示
view = nv.show_ase(atoms)
view.clear_representations()
view.add_ball_and_stick(selection="Fe")  # Feのみ表示

# インデックス範囲で選択
view.add_ball_and_stick(selection="0-99")  # 最初の100原子のみ
```

### パターン H: 距離・角度の測定

nglview の GUI パネルを使って、原子間距離や角度をインタラクティブに測定できます。

```python
import nglview as nv

view = nv.show_ase(atoms, gui=True)
# GUI パネルの "Measurements" タブで:
# - 2原子をクリック: 距離を測定
# - 3原子をクリック: 角度を測定
# - 4原子をクリック: 二面角を測定
view
```

プログラム的に測定する場合は ASE の関数を使用します。

```python
from ase.geometry import get_distances

# 原子0と原子1の距離（PBC考慮）
d = atoms.get_distance(0, 1, mic=True)  # mic=True で最小イメージ規約を使用
print(f"Distance: {d:.3f} A")

# 3原子間の角度
angle = atoms.get_angle(0, 1, 2)
print(f"Angle: {angle:.1f} degrees")

# 4原子間の二面角
dihedral = atoms.get_dihedral(0, 1, 2, 3)
print(f"Dihedral: {dihedral:.1f} degrees")
```

### パターン I: トラジェクトリファイルの読み込みと表示

`.traj` ファイルからトラジェクトリを読み込んで表示する一連の手順です。

```python
from ase.io import read
import nglview as nv

# .traj ファイルからの全フレーム読み込み
trajectory = read("md_trajectory.traj", index=":")
print(f"Total frames: {len(trajectory)}")

# 適切に間引いて表示
step = max(1, len(trajectory) // 200)
view = nv.show_asetraj(trajectory[::step], gui=True)
print(f"Displaying {len(trajectory[::step])} frames (skip: {step})")
view
```

## ベストプラクティス

### pfcc-extras の可視化を優先する

ASE と pfcc-extras の両方で可視化できる場合は、pfcc-extras を優先してください。`pfcc_extras.show_gui` は単一構造とトラジェクトリの両方に 1 行で対応でき、Matlantis 環境に最適化されています。`ase.visualize.view` は細かな制御が制限されるため、特別な理由がない限り使用を避けてください。

### 計算前の視覚チェック

計算投入前に必ず構造を可視化して確認してください。以下の問題は可視化で即座に発見できます。

- **結合の妥当性**: 不自然な結合がないか
- **周期境界のチェック**: PBC が正しく設定されているか、意図した通りに周期構造が繰り返されているか
- **原子の重なり**: 原子同士が異常に接近していないか（`Too many neighbors` エラーの原因）
- **真空層の確認**: 表面モデルの真空層が十分か（最低 10 A、推奨 15-20 A）

### 可視化は主タスクの補助

可視化のコードは簡潔に書き、メインの計算ワークフローの流れを妨げないようにしてください。1-2 行の可視化コードで十分な場合がほとんどです。

### 単一構造とトラジェクトリの使い分け

- 単一構造の確認: `show_gui(atoms)` または `nv.show_ase(atoms)`
- トラジェクトリの再生: `nv.show_asetraj(traj)` またはフレーム間引き
- 最適化過程の確認: `.traj` ファイルを読み込んでトラジェクトリ表示

### Notebook 固有の注意点

- nglview の表示は Notebook のセル出力として表示されます。セルを再実行すると前の表示は消えます
- 大量の可視化セルを含む Notebook はファイルサイズが大きくなります。`.ipynb` ファイルに画像データが埋め込まれるためです
- Notebook を保存する前に、不要な可視化出力をクリアすることを検討してください
- `view` オブジェクトを変数に保持しておくと、後から表示オプションを変更できます

### static vs interactive の選択

- **論文・レポート用**: POV-Ray による静止画レンダリングを使用してください
- **探索・デバッグ用**: nglview のインタラクティブ表示を使用してください
- **チーム共有用**: Notebook に nglview の出力を含めるか、POV-Ray 画像を Group Drive に保存してください

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| nglview が表示されない | Jupyter 拡張が未有効化 | `jupyter-nbextension enable nglview --py --sys-prefix` を実行する |
| 構造が空白で表示される | Atoms オブジェクトが空 | `len(atoms)` で原子数を確認する |
| トラジェクトリ表示が重い | フレーム数が多すぎる | `traj[::N]` でフレームを間引く |
| 原子が分断されて見える | PBC による wrap が原因 | `atoms.wrap()` を試す、または PBC を考慮した表示設定にする |
| Jupyter カーネルがクラッシュ | 大規模系でメモリ不足 | 表示する原子数を減らす、レプリゼンテーションを軽量にする |
| POV-Ray でレンダリングできない | povray がインストールされていない | `pip install povray` または環境に応じたインストールを行う |
| GUI パネルが表示されない | `gui=True` が指定されていない | `nv.show_ase(atoms, gui=True)` を使用する |
| 距離測定が PBC を考慮していない | `mic=False` がデフォルト | `atoms.get_distance(i, j, mic=True)` で PBC を考慮する |
| nglview の出力が Notebook に保存されない | Widget 状態が保存されていない | Notebook を保存する前に全セルを再実行し、Widget 状態を有効にする |
| 表示が更新されない | キャッシュの問題 | `view._remote_call('setSize', target='Widget', args=['100%', '400px'])` で再描画を試みる |

## 関連ガイド

- **反応経路探索** (reaction/SKILL.md): 反応経路の構造を可視化で確認
- **物性計算** (property-analysis/SKILL.md): 計算結果の構造確認
- **外部データ連携** (external-data/SKILL.md): 取得した構造の視覚検証
- **Group Drive ストレージ** (storage/SKILL.md): 可視化画像の保存・共有
