---
name: ase-basics
description: >
  ASE (Atomic Simulation Environment) の基礎操作を扱うスキルです。
  Atoms オブジェクトの生成・編集、ase.build.molecule, ase.build.bulk, repeat, supercell,
  PBC (periodic boundary conditions), cell, positions,
  get_chemical_formula, get_distance, translate, center,
  分子と結晶の違い、構造ファイルの読み書き(ase.io.read/write)
  に関するコード生成時に使用してください。
---
# ASE 基礎操作

## 概要

ASE (Atomic Simulation Environment) は原子スケールシミュレーションのための Python ライブラリで、Matlantis のすべての操作の基盤となります。このガイドでは、Atoms オブジェクトの生成・編集、分子と結晶の違い、構造ファイルの入出力、および基本的な操作パターンを解説します。

ASE を理解すれば Matlantis の操作方法がわかります。Matlantis は ASE の Calculator インターフェースを通じて PFP を利用する設計になっているためです。

### 中心的な概念

| 概念 | 説明 |
|------|------|
| `Atoms` | 原子配置を保持するオブジェクト。座標、元素、セル、周期境界条件を持つ |
| `Calculator` | エネルギーや力を返す計算手法。Atoms にアタッチして使用する |
| 分子 (Molecule) | 有限系。通常 `pbc=False` |
| 結晶 (Crystal) | 周期系。通常 `pbc=True`。セル (unit cell) を持つ |

分子と結晶の違いは、`pbc` (周期境界条件) と cell 最適化の有無で区別します。分子最適化では原子位置のみを扱い、結晶最適化では原子位置と cell の両方を扱うことが多いです。

## ワークフロー

```
1. Create  - Atoms オブジェクトを生成（builder または手動座標指定）
2. Edit    - 構造を編集（平行移動、真空層追加、距離確認など）
3. Save    - 構造ファイルに保存（.cif, .xyz, .traj 等）
4. Load    - 保存した構造を再読み込み
```

## 実装パターン

### パターン A: 分子の生成 (ase.build.molecule)

ASE に組み込まれた分子データベースから構造を取得します。

```python
from ase.build import molecule

# 水分子を作成
atoms = molecule("H2O")

# 基本情報の確認
print(len(atoms))                    # 原子数: 3
print(atoms.get_chemical_formula())  # 化学式: H2O
print(atoms.pbc)                     # 周期境界: [False, False, False]
```

molecule() で利用可能な分子名の例: `"H2O"`, `"CH4"`, `"CO2"`, `"NH3"`, `"C2H6"` など。

### パターン B: 結晶の生成 (ase.build.bulk)

結晶構造を格子定数と結晶型を指定して生成します。

```python
from ase.build import bulk

# FCC 銅を生成
atoms = bulk("Cu", "fcc", a=3.6)

# セル情報の確認
print(atoms.cell.lengths())  # セルの各辺の長さ
print(atoms.pbc)             # [True, True, True]
```

bulk() の主な結晶型: `"fcc"`, `"bcc"`, `"hcp"`, `"diamond"`, `"rocksalt"`, `"zincblende"` など。

### パターン C: 手動座標指定による Atoms 生成

任意の原子種と座標を直接指定して構造を作成します。

```python
from ase import Atoms

# CO2 分子を手動で作成
atoms = Atoms(
    "CO2",
    positions=[
        [0.0, 0.0, 0.0],    # C
        [0.0, 0.0, 1.16],   # O
        [0.0, 0.0, -1.16],  # O
    ],
)

print(atoms.get_chemical_formula())  # CO2
```

元素記号の文字列と positions リストの順序は一致している必要があります。

### パターン D: スーパーセルの作成 (repeat)

単位胞を繰り返してスーパーセルを作成します。大規模系のシミュレーションや、欠陥・ドーピングの導入に必要です。

```python
from ase.build import bulk

# Si ダイヤモンド構造の単位胞
atoms = bulk("Si", "diamond", a=5.43)
print(f"Unit cell: {len(atoms)} atoms")  # 2 atoms

# 2x2x2 スーパーセル
supercell = atoms.repeat((2, 2, 2))
print(f"Supercell: {len(supercell)} atoms")  # 16 atoms
```

`repeat((nx, ny, nz))` は各軸方向に nx, ny, nz 回繰り返します。セルサイズも自動的に拡大されます。

### パターン E: 周期境界条件 (PBC) の理解と設定

PBC (Periodic Boundary Conditions) は系の端で原子が反対側に繋がるかどうかを制御します。

```python
from ase.build import molecule, bulk

# 分子: PBC なし（有限系）
mol = molecule("H2O")
print(mol.pbc)  # [False, False, False]

# 結晶: PBC あり（周期系）
crystal = bulk("Cu", "fcc", a=3.6)
print(crystal.pbc)  # [True, True, True]

# PBC の手動設定
mol.pbc = False      # 全方向で非周期
crystal.pbc = True   # 全方向で周期
```

結晶やスラブ系では必ず `pbc=True` を設定してください。PBC を設定し忘れると物理的に無意味な結果になります。

### パターン F: セル情報の確認

```python
from ase.build import bulk

atoms = bulk("Cu", "fcc", a=3.6)

# セル各辺の長さ
lengths = atoms.cell.lengths()
print(f"Cell lengths: {lengths}")

# セル各辺の長さと角度
cell_params = atoms.get_cell_lengths_and_angles()
print(f"a, b, c, alpha, beta, gamma = {cell_params}")

# セルベクトル (3x3 行列)
print(atoms.cell[:])
```

### パターン G: 構造の編集

原子位置の平行移動、真空層の追加、原子間距離の測定など基本的な操作です。

```python
# 平行移動
atoms.translate((0.0, 0.0, 2.0))

# 真空層を追加してセル中央に配置（分子系で必須）
atoms.center(vacuum=8.0)

# 2 原子間の距離を取得
d = atoms.get_distance(0, 1)
print(f"Distance between atom 0 and 1: {d:.3f} A")
```

`center(vacuum=N)` は分子系でセルを設定する際に重要です。十分な真空層 (推奨 10A 以上) を確保することで、周期イメージとの不要な相互作用を防ぎます。

### パターン H: 構造ファイルの読み込み

ASE の `read` 関数は拡張子からフォーマットを自動判別します。

```python
from ase.io import read

# 単一構造の読み込み（最後のフレーム）
atoms = read("structure.cif")

# 最初のフレームを読み込み
atoms = read("structure.cif", index=0)

# 最後のフレームを明示的に読み込み
atoms = read("structure.cif", index=-1)

# 全フレームをリストとして読み込み（トラジェクトリ）
frames = read("opt.traj", index=":")
print(f"Total frames: {len(frames)}")
last_frame = frames[-1]
```

`index` パラメータ:
- `-1`: 最後のフレーム（デフォルト）
- `0`: 最初のフレーム
- `":"`: 全フレームをリストで取得

### パターン I: 構造ファイルの書き出し

```python
from ase.io import write

# CIF 形式で保存
write("output.cif", atoms)

# XYZ 形式で保存
write("output.xyz", atoms)

# VASP POSCAR 形式で保存
write("POSCAR", atoms, format="vasp")

# ASE Trajectory 形式で保存（エネルギー・力・応力も保持）
write("output.traj", atoms)
```

拡張子からフォーマットが自動判定されます。POSCAR のように拡張子がない場合は `format` 引数で明示します。

### パターン J: トラジェクトリの読み込みと操作

最適化や MD で生成された `.traj` ファイルから全フレームを読み込み、解析に使用します。

```python
from ase.io import read

# 全フレームを読み込み
traj = read("opt.traj", ":")

# フレーム数の確認
print(f"Total frames: {len(traj)}")

# 最後のフレーム（最適化後の構造）
last = traj[-1]

# 各フレームのエネルギーを取得
for i, frame in enumerate(traj):
    if frame.calc is not None:
        e = frame.get_potential_energy()
        print(f"Frame {i}: {e:.4f} eV")
```

### パターン K: 基本情報の確認

```python
# 原子数
print(f"Number of atoms: {len(atoms)}")

# 化学式
print(f"Formula: {atoms.get_chemical_formula()}")

# 元素リスト
print(f"Symbols: {atoms.get_chemical_symbols()}")

# 座標
print(f"Positions shape: {atoms.positions.shape}")

# 周期境界条件
print(f"PBC: {atoms.pbc}")
```

## ベストプラクティス

1. **pfcc-extras が提供する機能は ASE よりも pfcc-extras を優先する**: `pfcc_extras` は Matlantis 環境に最適化された高レベル API を提供しています。ASE と pfcc-extras の両方で実現できる操作について、`pip show pfcc-extras` を実行して pfcc-extras がインストールしてあることを確認した場合は pfcc-extras を優先して使用してください。主な優先対象は以下の通りです。
   - **可視化**: `ase.visualize.view` ではなく `pfcc_extras.show_gui` / `view_ngl` を使用
   - **PBC ラッピング**: `atoms.wrap()` ではなく `pfcc_extras.structure.wrap_molecule` を使用（分子断裂を防止）
   - **衝突検出**: 手動の距離比較ではなく `pfcc_extras.structure.collision.CollisionDetector` を使用
   - **近傍・結合解析**: 手動の距離計算ではなく `pfcc_extras.structure.get_neighbors` / `get_connectivity_matrix` を使用
   - **部分占有**: 手動処理ではなく `pfcc_extras.structure.PartialOccupancy` を使用
   - **軌跡変換**: MDTraj / MDAnalysis への変換には `pfcc_extras.trajectory.asetraj_to_mdtraj` / `asetraj_to_mdanalysis` を使用

2. **分子と結晶を混ぜない**: PBC、cell、stress の有無を明確に分けて扱ってください。分子系に stress を計算したり、結晶系で PBC を設定し忘れたりすると、物理的に無意味な結果になります。

2. **真空層の確保**: 分子系やスラブモデルでは `center(vacuum=10.0)` 以上の真空層を設けてください。15-20A が推奨です。

3. **Trajectory 形式の活用**: 計算結果の保存には `.traj` 形式が便利です。構造だけでなくエネルギー、力、応力情報もバイナリで保持されます。

4. **構造確認の習慣化**: 計算投入前に必ず構造を確認してください。`pfcc_extras.show_gui` を使えば Jupyter 環境上で 3D 可視化が可能です。

5. **スーパーセルサイズの検討**: シミュレーション目的に応じて適切なスーパーセルサイズを選択してください。原子数が多すぎると計算コストが増大し、少なすぎると有限サイズ効果が出ます。

6. **フォーマット自動判定の活用**: `read` / `write` では拡張子からフォーマットが自動判定されるため、`format` 引数は通常不要です。

### 補足: PBC が必要な系の判定基準

| 系の種類 | PBC 設定 | 理由 |
|---------|---------|------|
| 孤立分子 | `pbc=False` | 有限系。周期イメージとの相互作用は不要 |
| バルク結晶 | `pbc=True` (全方向) | 無限周期系の近似 |
| スラブ (表面) | `pbc=[True, True, False]` または `pbc=True` | 面内は周期、法線方向は真空層で分離 |
| ナノチューブ | `pbc=[False, False, True]` | 軸方向のみ周期 |

Matlantis では多くの場合 `pbc=True`（全方向周期）で設定し、非周期方向には十分な真空層を設けることで擬似的に有限系を表現します。

### 補足: Atoms オブジェクトの複製

構造を変更する前にコピーを取っておくと、元の構造を保持できます。

```python
# 浅いコピー（Calculator は共有される）
atoms_copy = atoms.copy()

# Calculator も含めたコピー
atoms_copy = atoms.copy()
atoms_copy.calc = atoms.calc
```

`atoms.copy()` は Calculator をコピーしません。Calculator も複製する場合は明示的に再設定してください。

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| `read()` でファイルが読めない | ファイルパスの誤り、またはフォーマット未対応 | パスを確認し、必要に応じて `format` 引数を明示してください |
| `pbc` が意図と異なる | `molecule()` は `pbc=False`、`bulk()` は `pbc=True` がデフォルト | 生成後に `atoms.pbc` を確認・設定してください |
| スーパーセル後に原子数が想定と異なる | `repeat()` の引数の理解不足 | `repeat((2,2,2))` は各軸 2 倍で原子数は 8 倍になります |
| `center()` 後にセルが大きすぎる | `vacuum` パラメータが大きすぎる | 分子系は `vacuum=10.0` 程度から始めてください |
| Trajectory 読み込みで 1 フレームしか得られない | `index` パラメータの指定漏れ | `read("file.traj", index=":")` で全フレームを取得してください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): PFP の設定方法
- **構造モデリング** (modeling/SKILL.md): 表面・分子生成の高度な手法
- **構造最適化** (optimization/SKILL.md): 構造最適化
- **分子動力学** (dynamics/SKILL.md): 分子動力学
