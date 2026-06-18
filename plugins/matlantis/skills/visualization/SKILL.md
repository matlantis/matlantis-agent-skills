---
name: mt-visualization
description: >
  原子構造とトラジェクトリの可視化を扱うスキルです。
  view_ngl, show_gui, AddEditor, pfcc_extras, nglview, show_ase, show_asetraj,
  ase.visualize.view, POV-Ray, 可視化, visualization,
  構造確認, 目視チェック, trajectory 表示
  に関するコード生成時に使用してください。
---

# 可視化

## 概要

Matlantis の JupyterLab 環境における原子構造とトラジェクトリの可視化ガイドです。シミュレーション結果の検証や構造確認は、計算ワークフロー全体の品質を左右する重要なステップです。

pfcc_extras パッケージが提供する以下の3ツールが主要な選択肢です。

- **pfcc_extras.view_ngl**: 構造・トラジェクトリの確認に最も手軽なビューア
- **pfcc_extras.show_gui**: view_ngl に加え、Calculator 連携によるエネルギー・力の表示が可能
- **pfcc_extras.visualization.AddEditor**: 原子の追加・削除・置換・最適化まで行える構造エディタ

上記3ツールでカバーできない場合（nglview の直接制御が必要な場合、高品質な静止画が必要な場合など）は、nglview や POV-Ray を補助的に使用します。  
可視化は計算の補助的な役割ですが、構造の妥当性確認を怠ると「原子が重なっている」「PBC が正しく設定されていない」といった初歩的なミスに気づかず、長時間の計算を無駄にすることがあります。計算投入前の可視化チェックは必須の手順です。

### 可視化手法の判断チェックリスト

| | view_ngl | show_gui | AddEditor |
| :--- | :--- | :--- | :--- |
| 主な目的 |構造・軌跡の確認 |構造確認 + エネルギー表示| 構造の編集・構築 |
| トラジェクトリ表示 | ✓ | ✓ | × | 
| エネルギー・力の表示 | × | ✓ | ✓ |
| 距離や角度の測定 | ✓ | ✓ | ✓ |
| 構造を図(pngなど)として保存 | × | ✓ | ✓ |
| 原子の追加・削除・置換 | × | × | ✓ |

## 実装パターン

### パターン A: pfcc_extras.view_ngl による可視化

Matlantis 環境で最も手軽に使える可視化手法です。単一の Atoms オブジェクトとトラジェクトリ（List[Atoms]）の両方に対応しています。  
Matlantis が提供する `pfcc_extras` パッケージに含まれています。

```python
from pfcc_extras import view_ngl

# 標準的な可視化
view_ngl(atoms, representations=["ball+stick"], replace_structure=True)  # Atoms オブジェクトもしくは ASE Trajectory オブジェクトを渡す

# 結合の描画をしない場合
view_ngl(atoms, replace_structure=True)
```
#### 用途

- 計算投入前の構造確認（結合の妥当性、PBC のチェック、原子の重なり検出）
- 構造最適化後の結果確認
- MD トラジェクトリの簡易再生

#### 入力形式
- `ase.Atoms`: 単一構造
- `List[ase.Atoms]`: トラジェクトリ（最適化過程、MD 過程など）
- ASE Trajectory オブジェクト

#### 表示オプションの調整
- `representations` オプションで表示スタイルを指定可能（例: `["ball+stick"]`）。省略した場合は結合は表示されません。
- `replace_structure=True` を指定すると、frame ごとに構造を丸ごと置き換えます。デフォルト（`False`）では frame 0 の構造情報を保持したまま更新するため、MD や NEB のように結合が変化するトラジェクトリでは意図しない結合が表示され続ける場合があります。

#### 注意事項

- JupyterLab 環境でのみ動作します
- 大規模系では表示が重くなる場合があります
- Notebook セル内で呼び出す必要があります（スクリプト実行では表示されません）

### パターン B: pfcc_extras.show_gui による可視化

`view_ngl` と同じく `pfcc_extras` が提供する可視化ツールです。
入力形式・注意事項は`view_ngl` と同様です。
Calculator を設定した状態で使うと、エネルギー・最大力の表示や簡易最適化まで GUI 上で実行できます。

```python
from pfcc_extras import show_gui

show_gui(atoms, representations=["ball+stick"], replace_structure=True)
```

#### 固有オプション

- **show_axes**: 座標軸の表示
- **show_atom_index**: 原子インデックスの表示
- **ball_size**: 原子球の表示サイズ

#### 画像の保存

Viewer 右側の control box にある "download image" / "save image" ボタンから、
現在の表示を PNG 画像として保存できます。

### パターン C: AddEditor による構造編集

nglview ベースのインタラクティブ構造エディタです。
可視化だけでなく、GUI 上で原子の追加・削除・置換・移動・回転と
構造最適化まで実行できます。モデリング作業に使用してください。

```python
from pfcc_extras.visualize.addeditor import AddEditor

# Atoms オブジェクトから起動
editor = AddEditor(atoms)
editor.display()

# SMILES 文字列から起動（引数なしの場合は CH4 が可視化されます）
editor = AddEditor("c1ccccc1")
editor.display()
```

atoms に Calculator が設定されていない場合は、自動設定されます。

#### 主な編集操作（タブ別）

| タブ | 主な機能 |
| :--- | :--- |
| Viewer | 視点変更、エネルギー・最大力の表示 |
| Range | XYZ 範囲・元素によるインデックス一括選択 |
| Editing | 選択原子の移動・回転（軸指定含む） |
| Addition | 原子・分子の追加、元素置換、削除 |
| Opt | 構造最適化（FIRELBFGS / LBFGS / FIRE） |
| Other | インデックス交換、拘束条件設定、Calculator 切替、画像保存 |

#### 注意事項

- トラジェクトリには対応していません。単一 Atoms のみ受け付けます。

### パターン D: nglview による単一構造の可視化

`view_ngl` や `show_gui` の内部で使用されているライブラリです。  
表示スタイルをより細かく制御したい場合や、`view_ngl` でカバーできない操作が必要な場合に直接利用します。

```python
import nglview as nv

# 単一構造の表示
view = nv.show_ase(atoms)
view
```

```
# トラジェクトリの表示
view = nv.show_asetraj(trajectory)
view
```

#### 表示オプションの調整

```python
view = nv.show_ase(atoms)

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

#### フレーム間引き

大規模なトラジェクトリではフレームを間引いて表示することで、メモリ使用量と表示速度を改善できます。

```python
# 10フレームごとに間引いて表示
view = nv.show_asetraj(trajectory[::10])
view
```

### パターン F: ASE visualization（ngl バックエンド）

> **注意**: `ase.visualize.view` よりも `pfcc_extras.view_ngl` または `nglview` の直接使用を推奨します。`view_ngl` はより簡潔で、Matlantis 環境に最適化されています。

ASE 組み込みの `view` 関数に nglview バックエンドを指定する方法です。

```python
from ase.visualize import view

# nglview バックエンドで表示
view(atoms, viewer="ngl")
```

この方法は ASE のインターフェースに慣れている場合に便利ですが、nglview を直接使う場合と比べて細かな制御は制限されます。レプリゼンテーションの変更や GUI パネルの表示などが必要な場合は、`nv.show_ase()` を直接使用してください。

### パターン G: POV-Ray による高品質レンダリング

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

### パターン H: 最適化前後の比較

構造最適化の前後で構造を比較する場合の手法です。

```python
from pfcc_extras import view_ngl
from ipywidgets import HBox

# 最適化前
view_before = view_ngl(atoms_before, representations=["ball+stick"], replace_structure=True)

# 最適化後
view_after = view_ngl(atoms_after, representations=["ball+stick"], replace_structure=True)

# 2つのビューを横並びで表示
HBox([view_before.v, view_after.v])
```

または、トラジェクトリとして最適化過程全体を表示します。

```python
# 最適化トラジェクトリの表示
from ase.io import read
opt_traj = read("optimization.traj", index=":")
view = view_ngl(opt_traj, representations=["ball+stick"], replace_structure=True)
view
```

### パターン I: 大規模系の可視化

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
# インデックス範囲で選択
view = nv.show_ase(atoms)
view.clear_representations()
sel = "@" + ",".join(str(i) for i in range(100))  # 最初の100原子のみ表示
view.add_ball_and_stick(selection=sel)
view
```

### パターン J: 距離・角度の測定

`view_ngl`, `show_gui`, `AddEditor`, `nglview` では、原子を右クリックで選択し、選択した最後の原子を再度クリックすることで距離、角度、二面角を測定できます。

```python
import nglview as nv

view = nv.show_ase(atoms, gui=True)
# 原子を右クリックして選択し、最後に選択した原子を再クリックすると測定結果が表示されます。
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

### パターン K: トラジェクトリファイルの読み込みと表示

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

ASE と pfcc-extras の両方で可視化できる場合は、pfcc-extras を優先してください。`pfcc_extras.view_ngl` は単一構造とトラジェクトリの両方に 1 行で対応でき、Matlantis 環境に最適化されています。

### 計算前の視覚チェック

計算投入前に必ず構造を可視化して確認してください。以下の問題は可視化で即座に発見できます。

- **結合の妥当性**: 不自然な結合がないか
- **周期境界のチェック**: PBC が正しく設定されているか、意図した通りに周期構造が繰り返されているか
- **原子の重なり**: 原子同士が異常に接近していないか（`Too many neighbors` エラーの原因）
- **真空層の確認**: 表面モデルの真空層が十分か（最低 10 A、推奨 15-20 A）

### 可視化は主タスクの補助

可視化のコードは簡潔に書き、メインの計算ワークフローの流れを妨げないようにしてください。1-2 行の可視化コードで十分な場合がほとんどです。

### 単一構造とトラジェクトリの使い分け

- 構造の確認: `view_ngl(atoms)`、`show_gui(atoms)` または `nv.show_ase(atoms)`
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
