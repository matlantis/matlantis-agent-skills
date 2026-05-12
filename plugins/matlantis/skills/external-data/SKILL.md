---
name: mt-external-data
description: >
  外部データベースからの構造取得とファイル形式変換を扱うスキルです。
  PubChem, Materials Project, mp-api, AFLOW, COD, Gaussian, VASP,
  外部データ, external data, データベース, database,
  CID, SMILES, SDF, CIF, POSCAR,
  API key, pubchem_atoms_search,
  構造取得, 構造インポート, DFT入力生成
  に関するコード生成時に使用してください。
---

# 外部データ連携

## 概要

Matlantis で使用する原子構造を外部データベースやファイル形式から取得し、PFP による計算に接続するためのガイドです。シミュレーションの出発点となる構造データは、目的に応じて適切なデータソースから取得し、ASE の `Atoms` オブジェクトに変換してから利用します。

本ガイドでは、以下の主要なデータソースとの連携方法を扱います。

- **PubChem**: 有機分子の 3D 構造（CID、名前、SMILES による検索）
- **Materials Project**: 無機結晶の構造データ（API キーが必要）
- **AFLOW**: 合金・無機結晶系の探索用データベース
- **COD (Crystallography Open Database)**: 結晶構造
- **Gaussian 16 / VASP**: DFT 計算コードへの入力生成

取得した構造は必ず PFP で再最適化することを推奨します。外部データベースの構造は異なるレベルの計算で最適化されているため、PFP のポテンシャル面上での平衡構造とは一致しない場合があります。

## ワークフロー

外部データ連携の標準的なワークフローは以下の通りです。

1. **データソースの特定**: 対象系（分子 or 結晶）に応じたデータベースを選択する
2. **構造の検索・取得**: API やファイルダウンロードで構造データを取得する
3. **ASE Atoms への変換**: 取得した構造を `ase.Atoms` オブジェクトに変換する
4. **構造の検証**: 可視化して構造の妥当性を確認する
5. **PFP による再最適化**: Matlantis で構造最適化を実行する

### データソース判断ツリー

| 対象系 | 推奨データソース | ファイル形式 |
| :--- | :--- | :--- |
| 有機分子 | PubChem | `.sdf`, `.mol` |
| 無機結晶（一般） | Materials Project | `.cif`, `.vasp` |
| 合金・無機結晶（探索） | AFLOW | `.cif` |
| 結晶構造（実験データ） | COD | `.cif` |
| DFT 検証用入力 | Gaussian 16 / VASP | `.gjf`, `POSCAR` |

## 実装パターン

### パターン A: PubChem からの分子構造取得

PubChem は有機分子の 3D 構造を提供するデータベースです。名前、化学式、SMILES、CID（Compound ID）で検索できます。

#### ASE ヘルパーによる取得

最もシンプルな方法は、ASE に組み込まれた `pubchem_atoms_search` 関数を使用することです。

```python
from ase.data.pubchem import pubchem_atoms_search

# 名前で検索
atoms = pubchem_atoms_search(name="aspirin")

# CID で検索
atoms = pubchem_atoms_search(cid=2244)

# SMILES で検索
atoms = pubchem_atoms_search(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
```

#### SDF ファイルからの読み込み

PubChem のウェブサイトから SDF ファイルをダウンロードした場合の読み込み方法です。

```python
import ase.io

# SDF ファイルの読み込み
atoms = ase.io.read("Conformer3D_CID_1140.sdf")
```

#### 孤立分子としての設定

PubChem から取得した分子は孤立系として計算するため、十分な真空層を設定します。

```python
# 真空層を設定（周期境界条件による隣接イメージとの相互作用を防ぐ）
atoms.center(vacuum=15.0)

# PBC を非周期に設定（孤立分子の場合）
atoms.pbc = False

# または PBC を有効にして十分な真空層を確保（PFP は PBC=True を推奨）
atoms.pbc = True
atoms.center(vacuum=15.0)
```

### パターン B: Materials Project からの結晶構造取得

Materials Project は無機結晶の構造データを体系的に提供するデータベースです。利用には API キーの取得が必要です。

#### API キーの取得

1. [Materials Project](https://materialsproject.org/) にアカウントを作成する
2. Dashboard から API キーを取得する
3. 環境変数またはコード内で API キーを設定する

#### mp-api による取得

```python
from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor

# Materials Project から構造を取得
with MPRester("YOUR_API_KEY") as mpr:
    # material ID で取得
    structure = mpr.get_structure_by_material_id("mp-149")  # Si

# Pymatgen Structure -> ASE Atoms に変換
atoms = AseAtomsAdaptor.get_atoms(structure)

# PBC を設定
atoms.pbc = True
```

#### 化学組成による検索

```python
from mp_api.client import MPRester

with MPRester("YOUR_API_KEY") as mpr:
    # 化学組成で検索
    docs = mpr.summary.search(
        chemsys="Li-Fe-O",
        fields=["material_id", "formula_pretty", "energy_above_hull"]
    )

    # 最も安定な構造を選択
    stable = min(docs, key=lambda x: x.energy_above_hull)
    structure = mpr.get_structure_by_material_id(stable.material_id)

atoms = AseAtomsAdaptor.get_atoms(structure)
```

### パターン C: CIF / VASP ファイルからの読み込み

外部データベースからダウンロードしたファイルを ASE で読み込む汎用的な方法です。

```python
from ase.io import read

# CIF ファイル
atoms_cif = read("structure.cif")

# VASP POSCAR
atoms_vasp = read("POSCAR")

# XYZ ファイル
atoms_xyz = read("structure.xyz")

# 複数構造を含むファイル（全フレーム読み込み）
all_frames = read("trajectory.xyz", index=":")
```

ASE の `read` 関数はファイル拡張子からフォーマットを自動判別します。

### パターン D: AFLOW データベースからの取得

AFLOW は合金や無機結晶系の構造データを豊富に持つデータベースです。

```python
from ase.io import read

# AFLOW からダウンロードした CIF ファイルを読み込み
atoms = read("aflow_structure.cif")
atoms.pbc = True
```

AFLOW のウェブインターフェースから構造を検索し、CIF 形式でダウンロードした後、`ase.io.read` で読み込むことを推奨します。

### パターン E: COD（Crystallography Open Database）からの取得

COD は実験的に決定された結晶構造のデータベースです。

```python
from ase.io import read

# COD からダウンロードした CIF ファイルを読み込み
atoms = read("cod_structure.cif")
atoms.pbc = True
```

### パターン F: Gaussian 16 入力ファイルの生成

PFP で探索した構造を Gaussian 16 で精密検証するための入力ファイル生成です。位置づけとしては「PFP で探索 → 外部 DFT で精密検証」というワークフローの一部です。

```python
from ase.io import write

# Gaussian 入力ファイルの生成
write(
    "input.gjf",
    optimized_atoms,
    method="B3LYP",
    basis="6-31G(d)",
)
```

### パターン G: VASP 入力ファイルの生成

PFP で得た構造を VASP 計算に渡す場合の入力ファイル生成です。

```python
from ase.io import write

# POSCAR の生成
write("POSCAR", optimized_atoms, format="vasp")

# sort=True で元素記号順にソート
write("POSCAR", optimized_atoms, format="vasp", sort=True)
```

### パターン H: PFP による再最適化

外部データベースから取得した構造は、PFP のポテンシャル面上で再最適化することを推奨します。

```python
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator
from ase.optimize import BFGS

def reoptimize_with_pfp(atoms, calc_mode="PBE", fmax=0.01):
    """
    外部データベースから取得した構造をPFPで再最適化する。

    再最適化の理由:
    - 外部DBの構造は異なるレベルの計算で最適化されている
    - PFPのポテンシャル面上での平衡構造と一致しない場合がある
    - 後続の計算（物性、反応経路等）の一貫性を確保するため
    """
    calc = ASECalculator(
        Estimator(model_version="v8.0.0", calc_mode=calc_mode)
    )
    atoms.calc = calc

    opt = BFGS(atoms)
    opt.run(fmax=fmax)

    return atoms
```

#### calc_mode の選択指針

| 対象系 | 推奨 calc_mode |
| :--- | :--- |
| 金属・合金・半導体 | `PBE` |
| 遷移金属酸化物 | `PBE_U` |
| 有機分子・分子性結晶 | `WB97XD` または `PBE_PLUS_D3` |
| 表面吸着系 | `PBE_PLUS_D3` |
| 高精度が必要な結晶 | `R2SCAN` または `R2SCAN_PLUS_D3` |

### パターン I: データ形式変換ユーティリティ

各種データソースから ASE Atoms への変換を統一的に行うユーティリティです。

```python
from ase.io import read
from ase import Atoms

def load_from_any_source(source, source_type="file", **kwargs):
    """
    各種データソースからASE Atomsを取得する統一インターフェース。

    Args:
        source: ファイルパス、CID、material_id など
        source_type: "file", "pubchem_name", "pubchem_cid",
                     "pubchem_smiles", "materials_project"
    """
    if source_type == "file":
        return read(source, **kwargs)

    elif source_type == "pubchem_name":
        from ase.data.pubchem import pubchem_atoms_search
        return pubchem_atoms_search(name=source)

    elif source_type == "pubchem_cid":
        from ase.data.pubchem import pubchem_atoms_search
        return pubchem_atoms_search(cid=int(source))

    elif source_type == "pubchem_smiles":
        from ase.data.pubchem import pubchem_atoms_search
        return pubchem_atoms_search(smiles=source)

    elif source_type == "materials_project":
        from mp_api.client import MPRester
        from pymatgen.io.ase import AseAtomsAdaptor
        api_key = kwargs.get("api_key", "YOUR_API_KEY")
        with MPRester(api_key) as mpr:
            structure = mpr.get_structure_by_material_id(source)
        return AseAtomsAdaptor.get_atoms(structure)

    else:
        raise ValueError(f"Unknown source_type: {source_type}")
```

## ベストプラクティス

### 取得方法と後処理を分離する

外部データの取得コードと、取得後の加工・計算コードを明確に分離してください。これにより、データソースの変更や取得方法の更新が容易になります。

### API キーの管理

Materials Project などの API キーをコードにハードコードしないでください。環境変数や設定ファイルで管理することを推奨します。

```python
import os
api_key = os.environ.get("MP_API_KEY", "YOUR_API_KEY")
```

### PFP での再最適化を必ず行う

外部データベースから取得した構造をそのまま物性計算や反応経路探索に使用しないでください。PFP のポテンシャル面上で再最適化することで、後続の計算の一貫性と精度を確保できます。

### 外部サービス固有の説明を長くしすぎない

Matlantis のドキュメントやチュートリアルでは、外部サービス（PubChem、Materials Project 等）の詳細な使い方ではなく、Matlantis との連携に焦点を当ててください。

### 構造取得後の可視化確認

取得した構造は必ず可視化して妥当性を確認してください。特に以下の点をチェックします。

- 原子の欠損がないか
- 結合が正しいか
- PBC の設定が適切か
- 原子の重なりがないか

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| `pubchem_atoms_search` で構造が見つからない | 名前のスペルミス、または PubChem に 3D 構造が登録されていない | CID や SMILES で再検索する。PubChem ウェブサイトで確認する |
| Materials Project API でエラー | API キーが無効、またはレート制限 | API キーを確認する。時間を置いて再試行する |
| CIF ファイルの読み込みエラー | CIF ファイルのフォーマットが不正 | `pymatgen` の `CifParser` で読み込みを試みる |
| 取得した構造の原子数が予想と異なる | 部分占有サイトがある、または非対称単位のみ | `pymatgen` で対称操作を適用する |
| PFP での再最適化が収束しない | 取得した構造が PFP のポテンシャル面から大きくずれている | fmax を緩めて段階的に最適化する。calc_mode を変更する |
| `WB97XD` で元素がサポートされていない | WB97XD は 9 元素（H,C,N,O,F,P,S,Cl,Br）のみ対応 | `PBE` または `PBE_PLUS_D3` に変更する |
| SDF ファイルの読み込みで座標がおかしい | 2D 構造データを読み込んでいる | PubChem から 3D Conformer を再取得する |
| `Too many neighbors` エラー | 取得した構造に原子の重なりがある | 構造を可視化して確認し、問題のある原子を修正する |

## 関連ガイド

- **可視化** (visualization/SKILL.md): 取得した構造の視覚的検証
- **反応経路探索** (reaction/SKILL.md): 反応物・生成物の構造取得
- **物性計算** (property-analysis/SKILL.md): 取得した構造からの物性計算
- **Group Drive ストレージ** (storage/SKILL.md): 取得した構造データの保存・共有
