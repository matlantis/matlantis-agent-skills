---
name: mt-modeling
description: >
  構造モデリングを扱うスキルです。表面(slab)切り出し、分子生成、原子置換、複雑系構造の構築を包括します。
  generate_all_slabs, pymatgen, get_orthogonal_c_slab, miller_index, fcc111, slab,
  SMILES, RDKit, Chem.MolFromSmiles, EmbedMolecule, 分子生成,
  原子置換, ドーピング, substitute,
  LiquidGenerator, EMC, ポリマー, 液体, アモルファス, HEA, 高エントロピー合金,
  PartialOccupancy, wrap_molecule, get_mol_list, CollisionDetector, 衝突検出,
  FixAtoms, 拘束条件, 表面モデル, 真空層
  に関するコード生成時に使用してください。
---
# 構造モデリング

## 概要

シミュレーションの目的に合わせた原子構造モデルを構築するためのガイドです。バルク結晶からの表面 (slab) 切り出し、SMILES 文字列からの分子生成、原子置換 (ドーピング)、液体・アモルファス・高エントロピー合金 (HEA) などの複雑系構造の作成、衝突検出、部分占有サイトの処理、外部データベースからの構造取得までを包括的に扱います。

モデリング品質がシミュレーション結果の信頼性を決定します。生成直後の構造は平衡状態にないため、物性値計算の前に必ず構造最適化を行ってください。

## ワークフロー

```
1. Select   - 目的に応じた構造構築ツールを選択
2. Generate - 初期構造を生成（表面切り出し、分子生成、液体充填など）
3. Modify   - 構造を修正（原子置換、拘束設定、真空層追加など）
4. Validate - 衝突検出、可視化による確認
5. Save     - 構造ファイルとして保存
```

### ツール選択ガイド

| 対象系 | 推奨ツール | 備考 |
|--------|-----------|------|
| バルク結晶・単純分子 | `ase.build.bulk`, `ase.build.molecule` | 基本の builder |
| 表面 (slab) | `ase.build.fcc111` 等の slab builder、pymatgen `generate_all_slabs` | 簡易 or 高度な切り出し |
| 低分子・短鎖の液体 | `LiquidGenerator` (`pfcc_extras`) | 凝集構造の自動生成 |
| 長鎖ポリマー | `EMC` | 長鎖用の専用ツール |
| 高エントロピー合金 (HEA) | 手組み HEA / special builder | ランダム配置 |
| アモルファス | Melt-Quench MD | MD による急冷 |
| クラスター | 手動配置 or Wulff construction | 表面エネルギーベース |

## 実装パターン

### パターン A: 表面 (Slab) の簡易生成 (ASE builder)

ASE 組み込みの slab builder で素早く表面モデルを作成します。

```python
from ase.build import fcc111

# Pt(111) 表面: 2x2 ユニットセル、4 層、真空 10A
atoms = fcc111("Pt", size=(2, 2, 4), vacuum=10.0)

print(f"Atoms: {len(atoms)}")
print(f"Cell: {atoms.cell.lengths()}")
```

`size=(nx, ny, nz)` は面内繰り返し (nx, ny) と層数 (nz) を指定します。`vacuum` は表面上下の真空層の厚さ (A) です。

### パターン B: 表面 (Slab) の高度な生成 (pymatgen)

任意のミラー指数に対応し、直交化とスーパーセル化を自動で行う堅牢な実装です。

```python
from pymatgen.core import Structure
from pymatgen.core.surface import generate_all_slabs
from ase import Atoms
import numpy as np

def create_slabs(
    bulk_atoms: Atoms,
    miller_index: tuple,
    min_slab_thick: float = 10.0,
    min_vacuum: float = 15.0,
    min_ab_size: float = 10.0
) -> list[Atoms]:
    """
    指定されたミラー指数のスラブ構造を作成する。

    Args:
        bulk_atoms: バルク結晶の Atoms オブジェクト
        miller_index: 切り出す面 (h, k, l) 例: (1, 1, 1)
        min_slab_thick: スラブ部分の最小厚さ (A)
        min_vacuum: 真空層の最小厚さ (A)
        min_ab_size: 表面 (ab 面) の最小サイズ (A)
    """
    # ASE -> Pymatgen 変換
    structure = Structure.from_ase_atoms(bulk_atoms)

    # 全終端面候補を生成
    slabs_gen = generate_all_slabs(
        structure,
        max_index=max(miller_index),
        min_slab_size=min_slab_thick,
        min_vac_size=min_vacuum
    )

    # 指定ミラー指数のみ抽出
    target_slabs = [s for s in slabs_gen if s.miller_index == miller_index]

    refined_slabs = []
    for slab in target_slabs:
        # C 軸の直交化
        slab = slab.get_orthogonal_c_slab()

        # Pymatgen -> ASE 変換
        atoms = slab.to_ase_atoms()

        # 面内サイズの拡張
        cell = atoms.get_cell()
        a_len = np.linalg.norm(cell[0])
        b_len = np.linalg.norm(cell[1])
        repeat_a = int(np.ceil(min_ab_size / a_len))
        repeat_b = int(np.ceil(min_ab_size / b_len))
        if repeat_a > 1 or repeat_b > 1:
            atoms *= (repeat_a, repeat_b, 1)

        # Z 軸方向のセンタリング
        atoms.center(axis=2)

        refined_slabs.append(atoms)

    return refined_slabs
```

### パターン C: SMILES からの分子生成

RDKit を使用して SMILES 文字列から 3 次元分子構造を生成します。

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms

def smiles_to_atoms_rdkit(smiles: str, vacuum: float = 10.0) -> Atoms:
    """SMILES 文字列から 3 次元分子構造 (Atoms) を生成する。"""
    # 分子オブジェクト生成
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # 水素付加
    mol = AllChem.AddHs(mol)

    # 3 次元配座の生成 (ETKDG 法)
    res = AllChem.EmbedMolecule(mol)
    if res == -1:
        # 生成失敗時はランダム座標で再試行
        res = AllChem.EmbedMolecule(mol, useRandomCoords=True)

    # 座標と元素情報の抽出
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # ASE Atoms 生成（真空層付き）
    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.center(vacuum=vacuum)

    return atoms
```

使用例:

```python
# エタノール
ethanol = smiles_to_atoms_rdkit("CCO")

# ベンゼン
benzene = smiles_to_atoms_rdkit("c1ccccc1")
```

### パターン D: 原子置換 (ドーピング)

結晶内の特定原子を別の元素に置換します。

```python
from ase import Atoms

def substitute_atom(atoms: Atoms, index: int, new_symbol: str) -> Atoms:
    """指定インデックスの原子を別の元素に置換する。"""
    new_atoms = atoms.copy()
    symbols = new_atoms.get_chemical_symbols()
    symbols[index] = new_symbol
    new_atoms.set_chemical_symbols(symbols)
    return new_atoms

# 使用例: 0 番目の原子を Pt に置換
# doped_atoms = substitute_atom(bulk_atoms, 0, "Pt")

# numbers 属性を直接変更する方法（原子番号で置換）
def substitute_atom_by_number(atoms: Atoms, index: int, new_atomic_number: int) -> Atoms:
    """指定インデックスの原子を原子番号で置換する。"""
    new_atoms = atoms.copy()
    new_atoms.numbers[index] = new_atomic_number
    return new_atoms

# 使用例: 0 番目の原子を Pt (78) に置換
# doped_atoms = substitute_atom_by_number(bulk_atoms, 0, 78)
```

置換位置が複数考えられる場合（サイト依存性）、すべての候補について構造を作成しエネルギー比較を行う必要があります。

### パターン E: 固定層の設定 (FixAtoms)

スラブモデルの下層を固定して上層のみ緩和する場合に使用します。

```python
import numpy as np
from ase.constraints import FixAtoms

# Z 座標が一定値以下の原子を固定
z_threshold = atoms.positions[:, 2].min() + 3.0  # 下から 3A 以内を固定
indices_fix = np.where(atoms.positions[:, 2] < z_threshold)[0]

constraint = FixAtoms(indices=indices_fix)
atoms.set_constraint(constraint)

# NumPy の boolean 配列でマスク指定する方法
mask_fix = atoms.positions[:, 2] < z_threshold  # dtype=bool, shape=(n_atoms,)
print(np.where(mask_fix)[0])  # マスクをインデックス配列に変換して表示
constraint = FixAtoms(mask=mask_fix)
atoms.set_constraint(constraint)
```

### パターン F: 液体構造の生成 (LiquidGenerator)

低分子や短鎖の凝集構造を自動生成します。

```python
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator

# 分子の Atoms オブジェクトを用意（例: 水分子）
from ase.build import molecule
atoms_mol = molecule("H2O")

# packmol バックエンドで充填（composition は dict のリスト）
generator = LiquidGenerator(
    "packmol",
    cell={"lx": 30, "ly": 30, "lz": 30},
    composition=[
        {"name": "water", "structure": atoms_mol, "number": 15}
    ]
)
atoms_liquid = generator.run()

# torch バックエンドで充填（composition は Atoms のリスト）
generator = LiquidGenerator(
    "torch",
    cell=[[30, 0, 0], [0, 30, 0], [0, 0, 30]],
    composition=[atoms_mol] * 15
)
atoms_liquid = generator.run()
```

### パターン G: 液体・アモルファスの Seed 構造

結晶構造を繰り返して高温 MD の出発点とする方法です。

```python
from ase.build import bulk

# Al の FCC 構造を 3x3x3 に拡大
atoms = bulk("Al", "fcc", a=4.05).repeat((3, 3, 3))
# この後、高温 MD -> 急冷でアモルファス化
```

乱数で初期配置を作って Seed とする方法もあります。

```python
import numpy as np
from ase import Atoms

# 30A 立方セルに Si を 64 原子ランダム配置
seed = 42
rng = np.random.default_rng(seed)
cell = np.diag([30.0, 30.0, 30.0])
positions = rng.random((64, 3)) @ cell
atoms = Atoms(symbols=["Si"] * 64, positions=positions, cell=cell, pbc=True)

# この後、衝突チェック -> 高温 MD -> 急冷でアモルファス化
```

### パターン H: 衝突検出

モデリング操作後に原子の異常な接近がないか検証します。

```python
from pfcc_extras.structure.collision import (
    exists_colliding_atom,
    generate_colliding_info_dataframe,
    CollisionDetector,
)

# 衝突の有無を真偽値で判定
has_collision = exists_colliding_atom(atoms)

# 衝突ペアの詳細情報をデータフレームで取得
if has_collision:
    df = generate_colliding_info_dataframe(atoms)
    print(df)  # 衝突ペアのインデックスと距離

# PBC を跨ぐ衝突の検出（複数分子系向け）
detector = CollisionDetector(atoms)
colliding_indices = detector.get_colliding_indices()
colliding_atoms = detector.extract_colliding_atoms()
```

### パターン I: 部分占有サイトの処理

部分占有 (partial occupancy) を含む結晶構造を実体化します。

```python
from pfcc_extras.structure.partial_occupancy import PartialOccupancy

# 部分占有サイトを明示的な原子配置に変換
po = PartialOccupancy(atoms)
realized_atoms = po.assign_atom_positions()
```

電荷中性を維持するために `charge_values` を指定して検証し、組成比がずれる場合は `cell_repeat` を増やして整数比近似を改善します。

### パターン J: 分子単位の PBC 整形

複数分子をセル内に配置した系で、分子が周期境界で分断されないようにします。

```python
from pfcc_extras.structure.molecule_wrap import wrap_molecule, get_mol_list

# 分子リストの取得
mol_list = get_mol_list(atoms)

# 分子単位で PBC ラッピング
wrapped_atoms = wrap_molecule(atoms)
```

原子単位ではなく分子単位で PBC 整形することで、分子断裂を防ぎます。

## ベストプラクティス

1. **pfcc-extras を優先する**: 構造モデリングにおいて ASE と pfcc-extras の両方で実現できる操作は、pfcc-extras を優先してください。特に衝突検出（`CollisionDetector`）、PBC ラッピング（`wrap_molecule`）、部分占有処理（`PartialOccupancy`）、近傍解析（`get_neighbors`）、SMILES からの分子生成（`pfcc_extras.smiles_to_atoms`）は pfcc-extras の専用 API がより安全・簡潔です。

2. **構造緩和は必須**: 生成直後の構造（特に Slab 切り出しや SMILES 生成後）は平衡状態にありません。物性値計算の前に必ず構造最適化を行ってください。

3. **真空層の厚さ**: 表面モデルや孤立分子では、周期イメージとの相互作用を防ぐため、最低 10A 以上、推奨 15-20A の真空層を確保してください。

4. **スラブの対称性**: 表面エネルギーを計算する場合、スラブの上下が同じ終端面である（対称性がある）ことが望ましいです。非対称な場合は双極子補正が必要になることがあります。

5. **衝突チェック**: モデリング操作後は必ず衝突検出を行ってください。見た目で妥当に見える構造でも接触原子が残ることがあり、計算時に `Too many neighbors` エラーの原因になります。

6. **可視化による確認**: 構造を生成・編集したら必ず可視化して確認してください。周期境界の不整合、意図しない原子配置、分子の断裂などを早期に発見できます。

7. **ドーピング位置の比較**: 置換位置が複数ある場合は、すべての候補構造を作成してエネルギー比較を行ってください。

8. **PBC の設定**: 結晶・スラブ系では `pbc=True` を忘れないでください。分子系では `pbc=False` が適切です。

9. **内包表記より NumPy ベクトル化を優先**: 数値配列（座標、距離、マスク、インデックス抽出など）に対する処理は、リスト内包表記より `NumPy` のベクトル演算（例: `np.where`, ブロードキャスト, boolean mask）を優先してください。可読性と性能の両面で有利です。

### 外部データベース

モデリングの出発点となる構造データは以下のデータベースから取得できます。

| データベース | 対象 | 形式 |
|-------------|------|------|
| Materials Project | 無機結晶全般 | `.cif`, `.vasp` |
| AFLOW | 合金・無機結晶系 | 各種 |
| PubChem | 有機分子 | `.sdf`, `.mol` |
| COD (Crystallography Open Database) | 結晶構造 | `.cif` |

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| `Too many neighbors` | 原子間距離が近すぎる（重なり） | 衝突検出で問題箇所を特定し、座標を修正してください |
| SMILES から 3D 座標生成失敗 | 複雑な分子構造 | `useRandomCoords=True` で再試行してください |
| スラブの終端面が意図と異なる | `generate_all_slabs` の結果から適切な候補を選択していない | `miller_index` でフィルタリングし、各候補を可視化して確認してください |
| 部分占有で組成がずれる | 部分占有の丸め誤差 | `cell_repeat` を増やして整数比近似を改善してください |
| 分子が PBC で分断される | 原子単位でラッピングしている | `wrap_molecule` で分子単位のラッピングを使用してください |
| `Cell is too small` | モデルのセルサイズが PFP のカットオフ半径に対して不足 | スーパーセル化してください |
| 表面モデルで上下の相互作用 | 真空層が不足 | `vacuum` を 15A 以上に設定してください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): Calculator 初期化
- **ASE 基礎操作** (ase-basics/SKILL.md): Atoms 生成・ファイル I/O
- **構造最適化** (optimization/SKILL.md): モデリング後の必須工程
- **分子動力学** (dynamics/SKILL.md): アモルファス作成の Melt-Quench 等
