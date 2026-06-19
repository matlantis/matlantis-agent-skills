---
name: pfcc-extras
description: >
  pfcc_extras (Matlantis/ASE workflow utilities) を使うコード生成時に使用します。
  可視化 (show_gui, view_ngl, SurfaceEditor)、構造加工 (smiles_to_atoms,
  LiquidGenerator, make_alloy, makesurface, 各種 interface)、MD スケジューラー
  (DepositionScheduler, ApplyUniformEfield など)、物性解析 (MonteCarlo, GCMC,
  calculate_bulk_modulus, 粘度・熱伝導率, trajectory 変換)、吸着・構造探索
  (adsorption_structure_search)、ジョブ実行制御 (run_jobs,
  ResourceAwareJobScheduler, parameter sweep) を扱います。
---

# pfcc_extras ユーティリティ

## 概要

`pfcc_extras` は Matlantis / ASE ワークフローを補完する実務ユーティリティ群です。可視化・構造加工・MDスケジューラー・物性解析・ジョブ実行制御といった、ASE 単体では不足しやすい機能をカバーします。以下のサンプルコード中で事前に生成していないオブジェクトに関しては事前に生成しておくようにすること。

| カテゴリ | 主な機能 |
|---|---|
| 可視化 | `show_gui`, `view_ngl`, `SurfaceEditor`, `AddEditor`, `traj_to_gif`, `traj_to_apng`, `pov_to_png`, `png_to_gif`, `save_traj_pov` |
| 構造加工 | `smiles_to_atoms`, `generate_conformers`, `LiquidGenerator`, `PartialOccupancy`, `wrap_molecule`, `CollisionDetector`, `make_alloy`, `select_subset`, `fill_interstitial_sites`, `makesurface`, `make_surfaces_pmg`, `make_rectangular_slab`, `make_mol_surface`, `get_hkl_bulk_structure` |
| MD スケジューラー | `DepositionScheduler`, `DeleteMoleculeScheduler`, `ElasticVirtualWall`, `ApplyUniformEfield`, `LinearTemperatureScheduler`, `CellDeformationScheduler`, `EarlyStopScheduler`, `TemperatureScaleScheduler` |
| 物性解析 | `IsosurfaceCalculator`, `calculate_bulk_modulus`, `ChemicalWindow`, `MonteCarlo` |
| 吸着・構造探索 | `adsorption_structure_search`, `IsotropicFilter` |
| ジョブ実行制御 | `run_jobs`, `ResourceAwareJobScheduler`, `QueueScheduler`, `PapermillScheduler`, `ParameterizedJob` |
| Light-PFP 構造ビルダー | `cut_sphere`, `cut_cube`, `solid_molecule_interface`, `liquid_liquid_interface`, `soak`, `solid_solid_interface_random` |
| Calculator ユーティリティ | `get_pfp_calculator`, `WrappedCalculator`, `TimeProfileHook` |
| アルゴリズム | `farthest_point_sampling`, `opt`, `opt_cell_size`, `opt_with_symmetry` |
| 軌跡変換 | `asetraj_to_mdtraj`, `asetraj_to_mdanalysis`, `modify_trajectory_by_mic` |

---

## 可視化

### show_gui

```python
from pfcc_extras.visualize import show_gui
atoms = read("input.traj")
show_gui(atoms, show_axes=True, show_atom_index=True, ball_size=0.3)
```

### view_ngl / view_ngl_traj

```python
from pfcc_extras.visualize.view import view_ngl

view_ngl(atoms, representations=["ball+stick"])  # representationsで結合表示も指定可能
```

### SurfaceEditor（表面系インタラクティブ編集）

Calculator を設定した状態で使うと、エネルギー・最大力の表示や簡易最適化まで GUI 上で実行できます。

```python
from pfcc_extras.visualize.surface_editor import SurfaceEditor

atoms.calc = calculator   # Calculator 設定が必要
editor = SurfaceEditor(atoms, w=600, h=500)
```

**主な機能**: カラースキーム・描画軸変更 / 選択原子の削除・置換・移動・回転 / エネルギー・最大力表示 / 簡易最適化（Run mini opt）/ 画像保存 / ポインターによる ID・座標表示

### AddEditor（分子系インタラクティブ編集）

```python
from pfcc_extras.visualize.addeditor import AddEditor
editor = AddEditor(atoms)
editor.display()
```


### POVRay アニメーション

`traj_to_gif` / `traj_to_apng` が最もシンプルなワンライナーです。内部的には `save_traj_pov` → `pov_to_png` → `png_to_gif` / `png_to_apng` のパイプラインになっており、個別関数を呼ぶと細かく制御できます。

```python
from pfcc_extras.visualize.povray import traj_to_gif, traj_to_apng

# GIF（ワンライナー）
traj_to_gif(
    atoms_list,
    gif_filepath="anim.gif",
    rotation="30x,30y",       # 固定回転（str）または lambda index, atoms: "..." で可変
    width=400,                 # 解像度（ピクセル）
    delay=100,                 # フレーム間隔（ms）
    n_jobs=8,                  # POV レンダリングの並列数
    clean=True,                # 中間ファイル（pov/png ディレクトリ）を自動削除
)

# APNG（GIF より高画質・大ファイル）
traj_to_apng(atoms_list, apng_filepath="anim.png", width=400, delay=100)
```

個別ステップで制御する場合:

```python
from pfcc_extras.visualize.povray import save_traj_pov, pov_to_png, png_to_gif

save_traj_pov(traj, outdir="pov", width=400, rotation="30x,30y")
pov_to_png(povdir="pov", pngdir="png", n_jobs=16)
png_to_gif(pngdir="png", gif_filepath="anim.gif", delay=100)
```

- `povray` が事前にシステムにインストールされている必要があります
- `rotation` に callable を渡すとフレームごとに回転角を変えられます（例: `lambda i, atoms: f"30x,{i}y"`）

---

## 構造加工

### SMILES / RDKit 変換（ase_rdkit_converter）

```python
from pfcc_extras.structure.ase_rdkit_converter import (
    smiles_to_atoms, atoms_to_smiles, smiles_to_rdmol, rdmol_to_atoms,
)

atoms = smiles_to_atoms("c1ccccc1", randomSeed=42)   # SMILES → Atoms
smiles = atoms_to_smiles(atoms)                       # Atoms → SMILES
```

- `randomSeed`（または `random_seed`）を固定して再現性を確保してください（デフォルト 1）
- 水素原子は自動付加（`AddHs`）されます
- 3D 座標の埋め込みに失敗する場合は `useRandomCoords=True` を指定するか、`maxAttempts` を増やしてください

```python
atoms = smiles_to_atoms("c1ccccc1", randomSeed=42, useRandomCoords=True, maxAttempts=1000)
```

### 配座探索（generate_conformers）

```python
from pfcc_extras.structure.ase_rdkit_converter import generate_conformers

mol, conf_ids = generate_conformers(
    smiles_or_atoms="CC(C)C",   # SMILES 文字列または ASE Atoms
    num_conformers=50,
    pruneRmsThresh=0.1,         # RMS 閾値（低いほど多様な配座を保持）
    energy_tolerance=10.0,      # MMFF エネルギーカットオフ（kcal/mol）
    seed=1234,
)
```

### 溶液構造生成（LiquidGenerator）

Packmol または Torch バックエンドで分子をパッキングした溶液系を作成します。

```python
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator

composition = [
    {"smiles": "O",   "number": 100},
    {"smiles": "CCO", "number": 20},
]
gen = LiquidGenerator(engine="packmol", composition=composition, density=1.0, tolerance=2.0)
liquid_atoms = gen.run()
```

`engine="packmol"` を使用してください。`packmol_bin` が存在しない場合は自動インストールされます。エラーが出る場合は、`density` を小さくしてみてください。

### 部分占有構造の離散化（PartialOccupancy）

```python
from pfcc_extras.structure.partial_occupancy import PartialOccupancy

po = PartialOccupancy(
    input_structure=atoms,
    csv_path="occupancy.csv",
    occupancy_header="occupancy",
    random_seed=42,
)
po.set_parameters(cell_repeat=(2, 2, 2))
realized_atoms = po.assign_atom_positions()
```

`cell_repeat` を大きくすると占有率の整数比近似が改善されます。

### 分子ラップと分子リスト（wrap_molecule / get_mol_list）

```python
from pfcc_extras.structure.molecule import wrap_molecule, get_mol_list

wrapped = wrap_molecule(atoms)       # 分子単位で PBC ラップ
```

原子単位で `wrap` すると分子がセル境界で分断されます。必ず分子単位で適用してください。

### 合金・元素置換（composition）

```python
from pfcc_extras.structure.composition import make_alloy, make_chemical_formula

formula = make_chemical_formula(atoms, elem_ratio={"Au": 1, "Pt": 3})
alloy_atoms = make_alloy(atoms, chemformula=formula, seed=42)
```

### 合金最安定構造探索（select_subset）

ランダム・貪欲法・アニーリング法で目的関数（エネルギー等）を最小化する原子配置を探索します。

```python
from pfcc_extras.structure.select_subset import (
    select_subset_random,
    select_subset_greedy,
    select_subset_annealing,
)
import numpy as np

indices = np.arange(len(atoms))

def objective(subset_indices):
    return get_energy(atoms, subset_indices)

# ランダム法（多数のランダム試行から最良を選択）
best, history, costs = select_subset_random(
    indices, n=20, objective_function=objective, max_iters=500, seed=42, n_jobs=-1,
)

# 貪欲法（1 サイトずつ最良の元素を確定）
best, history, costs = select_subset_greedy(
    indices, n=20, objective_function=objective, seed=42,
)

# アニーリング法（局所最適を脱出）
best, history, costs = select_subset_annealing(
    indices, n=20, objective_function=objective, initial_temp=500, final_temp=100, alpha=0.9,
)
```

戻り値: `(best_subset, subset_history, cost_history)`。`n_jobs=-1` で joblib 並列実行可能です。

### 衝突判定（CollisionDetector）

```python
from pfcc_extras.structure.connectivity import CollisionDetector

detector = CollisionDetector(atoms, collision_mult=0.8, connection_mult=1.05)

if detector.is_colliding().any():
    print(detector.get_info_as_dataframe())

clean_atoms = detector.extract_not_colliding_molecules()
```

- `collision_mult=0.8`: 共有結合距離 × 0.8 より近ければ「衝突」と判定
- `extract_colliding_atoms()` / `extract_not_colliding_molecules()` でクリーニングできます
- 手動配置後や `LiquidGenerator` 後の構造確認に使用してください

### 欠陥・格子間サイト生成（fill_interstitial_sites）

```python
from pfcc_extras.structure.defects import fill_interstitial_sites

atoms_with_interstitials, site_list, species, indices = fill_interstitial_sites(
    host_atoms=host,
    insert_species=["Li", "H"],
    mode="pymatgen",   # "pymatgen"（推奨）/ "voronoi" / "random"
    host_min_distance=1.5,
    seed=42,
)
```

### 表面スラブの簡易切り出し（makesurface）

ASE の `ase.build.surface` をベースにした最もシンプルなスラブ生成関数です。Miller 指数・層数・繰り返し・真空層を指定してバルクから 1 枚のスラブを切り出します。

```python
from pfcc_extras.structure.surface import makesurface

slab = makesurface(
    bulk_atoms,
    miller_indices=(1, 1, 1),  # Miller 指数
    layers=6,                  # スラブの層数
    rep=(4, 4, 1),             # 面内の繰り返し
    vacuum=30.0,               # 真空層の厚さ（Å、両側合計）
)
```

- スラブは Z 軸方向に切り出され、`vacuum/2` ずつ上下に真空層が付加されます
- 複数の Miller 面を一括生成したい場合や直交セルが必要な場合は、後述の `make_surfaces_pmg` / `make_rectangular_slab` を使用してください
- バルクが切断面でずれて切り出される場合は、事前に原子位置を微小シフトしてください

### 表面スラブ構築（make_surfaces_pmg / make_rectangular_slab）

```python
from pfcc_extras.structure.surface import make_surfaces_pmg
from pfcc_extras.structure.solid import make_rectangular_slab

# pymatgen ベースのスラブ（複数表面を一括生成）
slabs = make_surfaces_pmg(
    atoms=bulk_atoms, miller_index=(1, 1, 0),
    min_slab_size=10.0, min_vacuum_size=20.0, ab_rep=(2, 2), center_slab=True,
)

# MD 用直交セルスラブ
slab = make_rectangular_slab(
    atoms=bulk_atoms, miller_index=(1, 1, 1),
    min_slab_size=6.0, min_vacuum_size=10.0, max_atoms=3000, min_length=5.0,
)
```

### 有機分子結晶の表面構築（make_mol_surface）

バルクの有機分子結晶から表面スラブを構築します。構成分子を `mols` に渡すと、分子として完結しないフラグメントを自動的に除去します。

```python
from pfcc_extras.structure.surface import make_mol_surface

slab = make_mol_surface(
    atoms=bulk_atoms,          # 有機分子結晶のバルク構造
    mols=[mol_a, mol_b],       # バルクを構成する全分子（Atoms のリスト）
    hkl=(0, 0, 1),             # 表面法線方向（Miller 指数）
    slab_thickness=10,         # スラブ厚さ（Å）
    vacuum=10,                 # 真空層の厚さ（Å）
)
```

- `mols` にはバルク構造を構成する全分子を渡してください。漏れがあると表面に欠損が生じます
- 表面に出現するフラグメント（分子として不完全な原子群）は自動除去されます

### NPT 用セルの上三角変換（convert_atoms_to_upper）

ASE の NPT モジュールはセルが上三角形式（`cell[1,0] == cell[2,0] == cell[2,1] == 0`）であることを要求します。この条件を満たさない構造を NPT に渡す前に使用してください。

```python
from pfcc_extras.structure.rotate import convert_atoms_to_upper

atoms_upper = convert_atoms_to_upper(atoms)
# atoms_upper.cell は上三角形式に回転済み
```

### バルク構造の hkl 変換（get_hkl_bulk_structure）

指定した Miller 指数面を Z 軸とするバルク構造に変換します。界面 MD の初期構造生成に使います。

```python
from pfcc_extras.structure.solid import get_hkl_bulk_structure

bulk_list = get_hkl_bulk_structure(atoms=bulk_atoms, miller_index=(1, 1, 1), max_atoms=3000)
hkl_bulk = bulk_list[0]
```

### 近傍・結合解析

```python
from pfcc_extras.structure.connectivity import get_neighbors, get_connectivity_matrix

neighbors = get_neighbors(atoms, r=3.0) # r が None であれば、covalent radii を基に計算されます
connectivity = get_connectivity_matrix(atoms, mult=1.05)
```

`cutoff` の閾値はワークフロー全体で統一し、結合判定の揺らぎを防いでください。

---

## MD スケジューラー

### DepositionScheduler（分子堆積）

```python
from pfcc_extras.molecular_dynamics.scheduler import (
    DepositionScheduler, convert_kinetic_energy_to_velocity,
)

scheduler = DepositionScheduler(
    atoms=slab_atoms,
    dyn=dyn,
    molecules=[molecule_atoms],
    fractions=[1.0],
    num_total_steps=10000,
    incident_energy=1.0,     # eV（または incident_velocities で直接指定）
    initial_height=10.0,
    axis=2,
    seed=42,
)
dyn.attach(scheduler, interval=1)
```

- `incident_energy` → `incident_velocities` 変換は `convert_kinetic_energy_to_velocity(atoms, energy)` を使うとケース間の比較が容易になります
- `fractions` の総和は自動正規化されます

### DeleteMoleculeScheduler（セル外分子の自動削除）

```python
from pfcc_extras.molecular_dynamics.scheduler import DeleteMoleculeScheduler

deleter = DeleteMoleculeScheduler(atoms=atoms, dyn=dyn, axis=2)
dyn.attach(deleter, interval=100)
```

### ElasticVirtualWall（弾性仮想壁）

非周期系の MD でセル境界に弾性壁を設け、原子・分子を反射させます。

```python
from pfcc_extras.molecular_dynamics.scheduler import ElasticVirtualWall

wall = ElasticVirtualWall(
    atoms=atoms, dyn=dyn,
    wall_position=30.0,    # 壁の位置（Å）
    axis=2,
    wall_direction="top",  # "top" または "bottom"
    distance_tol=2.0,      # 分子グループ化の距離閾値（None で原子単位）
    seed=42,
)
dyn.attach(wall, interval=1)
```

### ApplyUniformEfield（一様電場）

```python
from pfcc_extras.molecular_dynamics.uniform_electric_field import ApplyUniformEfield

# charges を渡さない場合：MD ステップごとに atoms.get_charges() を呼び出して動的に更新
efield_constraint = ApplyUniformEfield(atoms=atoms, efield=(0.0, 0.0, 0.01))

# charges を渡す場合：固定値を使用（len(charges) == len(atoms) が必要）
efield_constraint = ApplyUniformEfield(atoms=atoms, efield=(0.0, 0.0, 0.01), charges=[...])

atoms.set_constraint(efield_constraint)
```

- `charges` を省略すると各 MD ステップで `atoms.get_charges()` を呼び出して電荷を動的に取得します。ポテンシャルモデルが電荷を計算できる場合（例: ReaxFF 系）に有効です
- `charges` を固定する場合は `len(charges) == len(atoms)` を確認してください

### LinearTemperatureScheduler（線形温度スケジュール）

昇温・降温・急冷に使います。`estimate_rate` / `estimate_steps` で必要なレートやステップ数を事前に計算できます。

```python
from pfcc_extras.molecular_dynamics.scheduler import (
    LinearTemperatureScheduler, estimate_rate, estimate_steps,
)

# 必要なレート・ステップ数を推定
rate  = estimate_rate(start=300, end=1500, nsteps=10000, timestep_fs=2.0)  # K/fs
steps = estimate_steps(start=300, end=1500, rate=0.12, timestep_fs=2.0)

scheduler = LinearTemperatureScheduler(
    dyn=dyn, temp_start=300, temp_end=1500,
    nsteps=10000,   # または temp_rate=rate で K/fs 指定も可
)
scheduler.attach(interval=1)
```

### CellDeformationScheduler（セル変形スケジュール）

引張・圧縮 MD に使います。

```python
from pfcc_extras.molecular_dynamics.scheduler import CellDeformationScheduler

target_cell = atoms.get_cell().copy()
target_cell[2, 2] *= 1.5   # Z 軸を 1.5 倍に伸長

scheduler = CellDeformationScheduler(
    dyn=dyn, atoms=atoms, target_cell=target_cell,
    nsteps=5000,                    # または strain_rate=1e8（1/s）で指定
    mask=[False, False, True],      # Z 軸のみ変形
)
dyn.attach(scheduler, interval=1)
```

### EarlyStopScheduler（条件付き早期終了）

```python
from pfcc_extras.molecular_dynamics.scheduler import EarlyStopScheduler

def stop_condition(dyn, atoms):
    return atoms.get_temperature() > 1500

stopper = EarlyStopScheduler(
    dyn=dyn, atoms=atoms, func=stop_condition, post_stop_steps=100,
)
dyn.attach(stopper, interval=10)
```

### TemperatureScaleScheduler（温度スケーリング）

```python
from pfcc_extras.molecular_dynamics.scheduler import TemperatureScaleScheduler

scheduler = TemperatureScaleScheduler(atoms, dyn=dyn, temperature=300.0, indices=None, mask=None)
dyn.attach(scheduler, interval=100)
```

---

## 物性解析・化学ポテンシャル

### イオン伝導パスの等値面計算（IsosurfaceCalculator）

MD 軌跡から特定元素の確率分布グリッドを計算し、VESTA / OVITO / NGLViewer で読み込める Gaussian cube ファイルを出力します。

```python
from pfcc_extras.isosurface.isosurface import IsosurfaceCalculator
from ase.io import Trajectory

iso = IsosurfaceCalculator(Trajectory("md.traj"), grid_size=0.5)
iso.calculate(symbol="Li")
iso.export("Li_isosurface.cube")   # Gaussian cube 形式で出力
```

`grid_size` はいずれの格子定数より小さい値を設定してください。

### バルク弾性率・状態方程式（calculate_bulk_modulus）

```python
from pfcc_extras.analysis.eos import calculate_bulk_modulus

def get_calc():
    return get_pfp_calculator(calc_mode="crystal_u0", model_version="latest")

V0, E0, B = calculate_bulk_modulus(
    atoms=atoms,
    volume_ratio_range=(0.9, 1.1),
    get_calculator=get_calc,
    num_data_points=10,
    fmax=0.005,
    n_jobs=4,           # -1 で全 CPU
)
# V0: 平衡体積（Å³）, E0: 最小エネルギー（eV）, B: バルク弾性率（GPa）
```

`get_calculator` は毎回新しい calculator を返すファクトリ関数を渡してください。

### 化学ポテンシャル窓（ChemicalWindow）

Materials Project データと PFP エネルギーを組み合わせ、安定相図上の化学ポテンシャル許容域を計算・可視化します。

```python
from pfcc_extras.chemical_potentials import ChemicalWindow

cw = ChemicalWindow(elements=["Li", "Ni", "O"])

# データの取得（Materials Project からのダウンロード、または自前の PFP 計算エネルギーの読み込み）
cw.download_data_from_materials_project(api_key="YOUR_MP_API_KEY", thermo_types=["GGA_GGA+U"])
# cw.read_csv("energy_per_atoms.csv")   # 代わりに CSV から読み込む場合

# 化学ポテンシャル窓の計算
cw.calculate_chemical_potentials(x_element="Li", y_element="Ni")
cw.get_chemical_potentials_under_constraints(constraints={...})   # 必要に応じて制約を指定
cw.calculate_formation_energy()

# 可視化（plot_type は 'heatmap' / 'contour'）
cw.plot_2d(width=500, height=500, plot_type="contour", title="Chemical window")
```

### モンテカルロシミュレーション（MonteCarlo）

格子 MC（LMC）・オフ格子 MC・グランドカノニカル MC（GCMC）を統一インターフェースで提供します。

```python
from pfcc_extras.monte_carlo.mc import MonteCarlo

mc = MonteCarlo(
    base_str=base_str,
    calculator=calculator,
    total_mc_iterations=10000,
    pressure=1.0,           # atm
    temperature=300.0,      # K
    mc_type_ratios={"lattice": 1.0},   # "lattice" / "off_lattice" / "gcmc" の混合比（合計 1.0）
    random_seed=42,
    output_data_path="output_data",
)

# Lattice MC
mc.setup_lattice_mc(atom_ids_to_swap={"Fe": 0.5, "Ni": 0.5})

# GCMC
mc.setup_gcmc(
    adsorbate_species=adsorbate_species,
    gcmc_move_ratios={"insert": 0.25, "remove": 0.25, "translate": 0.25, "rotate": 0.25},
    mu_TP0=mu_TP0,
)

mc.run()
```

`mc_type_ratios` の合計が 1.0 になるよう設定してください。`random_seed` を固定して再現性を確保してください。

### シフトエネルギー補正

PFP のエネルギーを VASP / Gaussian の参照基準に揃えます。`pfcc_extras.shift_energies` は **v0.8.1 以降 deprecated** で、代替として `Estimator.get_shift_energy_table()` を使用してください。

```python
import numpy as np
shift_table = estimator.get_shift_energy_table()   # {原子番号: shift_energy}
shift_energy = np.sum([shift_table[n] for n in atoms.get_atomic_numbers()])
vasp_aligned_energy = atoms.get_potential_energy() + shift_energy
```

---

## 吸着・構造探索

### Optuna ベース吸着構造探索（adsorption_structure_search）

Optuna のベイズ最適化で吸着分子の位置・姿勢を網羅的に探索します。スラブ・クラスター・多孔質結晶の各バリアントがあります。

```python
from pfcc_extras.adsorption.adsorption_structure_search import (
    adstructure_search_for_slab,
    adstructure_search_for_cluster,
    adstructure_search_for_porus,
    cluster_position_search_on_slab,
)

# スラブへの吸着
study = adstructure_search_for_slab(
    calc_mode="crystal_u0_plus_d3", model_version="latest",
    molec_path="molec_opt.traj", slab_path="Slab.traj",
    TH_max_f=5.0, TH_min_f=0.0005, tol=0.5,
    output_path="output/", study_name="slab_study", n_trials=100,
)

# クラスター上への吸着
study = adstructure_search_for_cluster(
    molec_path="opt_molecule.traj", clus_path="opt_cluster.traj",
    fix_ctr_cluster=4.0, output_path="output/", n_trials=100,
)

# 多孔質結晶への吸着
study = adstructure_search_for_porus(
    molec_path="opt_molecule.traj", crystal_path="porus_crystal.traj",
    output_path="output/", n_trials=100,
)

# スラブ上クラスターの位置探索
study = cluster_position_search_on_slab(
    clus_path="opt_cluster.traj", slab_path="Slab.traj",
    output_path="output/", n_trials=100,
)
```

- `TH_max_f` / `TH_min_f` で無効な初期配置を早期スキップし、計算コストを削減できます
- ストレージは `JournalFileBackend` を使い、中断再開が可能です（`load_if_exists=True`）
- Optuna の `Terminator` を使うと、改善が見られない場合に自動停止できます

```python
import optuna
from optuna.terminator import Terminator, TerminatorCallback, BestValueStagnationEvaluator

improvement_evaluator = BestValueStagnationEvaluator(max_stagnation_trials=150)
terminator = Terminator(improvement_evaluator=improvement_evaluator)

study = adstructure_search_for_slab(
    calc_mode="crystal_u0_plus_d3",
    ...,
    n_trials=1000,                              # 大きめに設定しておく
    sampler=optuna.samplers.TPESampler(),
    callbacks=[TerminatorCallback(terminator)], # 150 回改善なしで自動停止
)
```

### 等方圧フィルター（IsotropicFilter）

NPT でセルを等方的に変形させます。剪断成分を除去したい場合に使います。

```python
from pfcc_extras.filters.isotropic_filter import IsotropicFilter
from ase.optimize import LBFGS

opt = LBFGS(IsotropicFilter(atoms))
opt.run(fmax=0.005)
```

応力テンソルの剪断成分（`stress[3:]`）を 0 に固定し、対角成分を平均化して等方圧力として扱います。

---

## ジョブ実行制御

### run_jobs（統一インターフェース）

`run_jobs` は Notebook ジョブを実行する統一関数です。`jobs` には **Notebook のパス文字列** を渡します（`Job` オブジェクトではありません）。フラットリスト・グループ化リスト・パラメータスイープの 3 形態に対応します。

```python
from pfcc_extras.job_scheduler.runner import run_jobs

# 1) フラットリスト：複数 Notebook を並列投入
results = run_jobs(
    jobs=["01.ipynb", "02.ipynb", "03.ipynb"],
    max_workers=2,        # 並列投入数
    sequential=False,     # True で各ジョブの完了を待つ
    my_limit=0.4,         # 個人のトークンレート上限（0.0-1.0、None で無効）
    tenant_limit=0.7,     # テナントのトークンレート上限
    polling_interval_sec=5,
    post_submit_wait_sec=5,
    status_check_interval_sec=5,
)

# 2) グループ化リスト：グループ単位で逐次（優先度付き）実行
results = run_jobs(
    jobs=[["01.ipynb", "02.ipynb"], ["03.ipynb"]],  # 前のグループ完了後に次へ
    max_workers=2,
    sequential=True,
)

# 3) パラメータスイープ：同一 Notebook を複数パラメータで実行
results = run_jobs(
    jobs="template.ipynb",          # 単一の Notebook パス
    parameters=[
        {"temperature": 300},
        {"temperature": 500},
    ],
    max_workers=2,
    sequential=False,
)
```

`run_jobs` は `JobResult` のリストを返します。`JobResult.status` で `FAILED` と `SKIPPED`（出力ファイルが既に存在）を分離処理してください。

### スケジューラクラスを直接使う

`run_jobs` を使わず、スケジューラクラスを直接生成して `scheduler.run(...)` で実行することもできます。共通の引数は `my_limit` / `tenant_limit` / `mem_limit_mb` / `cpu_limit_percent` / `max_workers` / `polling_interval_sec` / `post_submit_wait_sec` です（`resource_limits` や `monitor_interval`、`max_concurrent` という引数は存在しません）。

```python
from pfcc_extras.job_scheduler.scheduler import (
    ResourceAwareJobScheduler, QueueScheduler, PapermillScheduler,
)

# 並列投入（リソース監視つき）
scheduler = ResourceAwareJobScheduler(
    my_limit=0.4, tenant_limit=0.7,
    max_workers=2, polling_interval_sec=5, post_submit_wait_sec=5,
)
results = scheduler.run(["01.ipynb", "02.ipynb"])

# キュー実行（完了待ち。グループ化リストも可）
scheduler = QueueScheduler(
    my_limit=0.4, tenant_limit=0.7,
    max_workers=2, polling_interval_sec=5, post_submit_wait_sec=5,
    status_check_interval_sec=5,
)
results = scheduler.run(["01.ipynb", "02.ipynb"])

# パラメータスイープ（papermill）
scheduler = PapermillScheduler(
    my_limit=0.4, tenant_limit=0.7,
    max_workers=2, sequential=False,
    polling_interval_sec=5, post_submit_wait_sec=5, status_check_interval_sec=5,
)
results = scheduler.run(
    infile="template.ipynb",
    parameters_list=[{"learning_rate": 0.01}, {"learning_rate": 0.001}],
)
```

---

## LightPFP 向け構造ビルダー

LightPFP 学習データ生成・大規模 MD の初期構造構築向けのユーティリティです。LightPFP 以外の用途でも使えます。

### 球形・立方体クラスターの切り出し（cut_sphere / cut_cube）

```python
from pfcc_extras.light_pfp.cluster import cut_sphere, cut_cube
from ase.build import bulk

cluster_sphere = cut_sphere(bulk("Fe") * (10, 10, 10), radius=10.0, vaccum=5.0)
cluster_cube   = cut_cube(bulk("Fe") * (10, 10, 10),   length=20.0, vaccum=5.0)
```

### 固体-液体界面（solid_molecule_interface）

```python
from pfcc_extras.light_pfp.solid_liquid_interface import solid_molecule_interface
from ase.build import fcc111

atoms = solid_molecule_interface(
    solid=fcc111("Cu", size=(5, 5, 3), vacuum=10.0),
    molecules=["c1ccccc1", "CCO"],    # SMILES または Atoms
    n_molecules=[1, 5],               # 分子比率
    molecule_layer_thickness=15.0,    # 液体層の厚さ（Å）
    molecule_density=0.8,             # 密度（g/cm³）
    spacing=1.5,
)
```

### 液体-液体界面（liquid_liquid_interface）

```python
from pfcc_extras.light_pfp.liquid_liquid_interface import liquid_liquid_interface

atoms = liquid_liquid_interface(
    molecules_left=["O", "CCO"],               # 下層
    molecules_right=["CCCCCCC", "c1ccccc1"],   # 上層
    n_molecules_left=[2, 1],
    n_molecules_right=[1, 1],
    xy_length=15.0,
    molecule_layer_thickness_left=10.0,
    molecule_layer_thickness_right=10.0,
)
```

### 固体を液体に浸す（soak）

```python
from pfcc_extras.light_pfp.soak import soak

soaked = soak(solid=cluster_atoms, fluid=liquid_atoms, spacing=2.5)
```

### 固体-固体界面（solid_solid_interface_random）

```python
from pfcc_extras.light_pfp.solid_solid_interface import solid_solid_interface_random

atoms_gb = solid_solid_interface_random(
    atoms_left=bulk("Fe") * (5, 5, 5),
    atoms_right=bulk("Cu") * (5, 5, 5),
    length=15.0, spacing=1.5,
)
```

### Materials Project からの構造取得・元素置換

```python
from pfcc_extras.light_pfp.solid import search_materials, replace_elements

atoms_list = search_materials(["mp-149", "mp-6930"], api_key="YOUR_KEY")

new_atoms = replace_elements(
    atoms_list[0] * (3, 3, 3),
    replace_scheme={"O": [("N", 0.20), ("F", 0.10)], "Fe": [("Ni", 0.30)]},
)
```

### rNEMD（粘度・熱伝導率計算）

MD 実行中のモメンタム交換（`RNEMDExtension`）と、終了後のログ解析（`get_viscosity` / `get_thermal_conductivity`）がセットになっています。

```python
from pfcc_extras.light_pfp.rnemd_utils import (
    RNEMDExtension,
    VelocityLogExtension,
    TemperatureLogExtension,
    get_viscosity,
    get_thermal_conductivity,
    plot_profile,
)

# MD 中のモメンタム交換 + プロファイルログ記録
rnemd = RNEMDExtension(
    md=md_dynamics,
    rnemd_type="viscosity",          # "viscosity" or "thermal_conductivity"
    n_slab=10,                       # スラブ数（偶数必須）
    n_swap=1, swap_interval=10,
    logfile="rnemd.log",
)
vel_log  = VelocityLogExtension(md_dynamics,   n_slab=10, logfile="velocity.log")
temp_log = TemperatureLogExtension(md_dynamics, n_slab=10, logfile="temperature.log")

md_dynamics.attach(rnemd,    interval=1)
md_dynamics.attach(vel_log,  interval=100)
md_dynamics.attach(temp_log, interval=100)

md_dynamics.run(100000)

# MD 終了後の解析
viscosity = get_viscosity(
    rnemd_log="rnemd.log", profile_log="velocity.log",
    timestep=2.0, init_steps=10000,   # 平衡化ステップを除外
)
print(f"粘度: {viscosity:.4f} mPa·s")

thermal_cond = get_thermal_conductivity(
    rnemd_log="rnemd.log", profile_log="temperature.log",
    timestep=2.0, init_steps=10000,
)
print(f"熱伝導率: {thermal_cond:.4f} W/(m·K)")

# 速度/温度プロファイルの可視化
plot_profile("velocity.log", figname="velocity_profile.png", init_steps=10000)
```

- `n_slab` は偶数必須です
- `init_steps` で平衡化ステップを除外して定常状態のデータのみ使用してください

---

## 軌跡変換

ASE の軌跡を MD 解析ライブラリ（MDTraj / MDAnalysis）に変換します。

```python
from pfcc_extras.structure.ase_traj_converter import asetraj_to_mdtraj, asetraj_to_mdanalysis

md_traj  = asetraj_to_mdtraj(traj, set_PBC=True, bond_cutoff=1.2)
universe = asetraj_to_mdanalysis(traj, set_PBC=True, bond_cutoff=1.2)
```

- 全フレームで原子数・元素順を一致させてください
- `bond_cutoff` は系ごとに調整して過結合・未結合の両方を防いでください

### 軌跡への MIC 補正（modify_trajectory_by_mic）

PBC 系の MD 軌跡では原子がセル境界をまたぐと座標が不連続に飛びます。`modify_trajectory_by_mic` は最小像規則（MIC）を各フレーム間に適用し、原子が連続的に移動するように補正します。拡散係数の計算や軌跡の可視化の前処理として使用してください。

```python
from pfcc_extras.structure.boundary import modify_trajectory_by_mic

traj_corrected = modify_trajectory_by_mic(traj)
```

---

## Calculator ユーティリティ

### get_pfp_calculator（PFP / EMT フォールバック）

```python
from pfcc_extras.calculators.pfp_calculator import get_pfp_calculator

calc = get_pfp_calculator(calc_mode="crystal_u0", model_version="latest", priority=50)

# config dict で一括指定
calc = get_pfp_calculator(config={"calc_mode": "crystal_u0", "model_version": "v8.0.0"})

# EMT で強制フォールバック（テスト用）
calc = get_pfp_calculator(model_version="emt")
```

`pfp_api_client` が利用可能なら PFP を、そうでなければ EMT を返します。ファクトリ関数として他の関数に渡す用途に適しています。

### WrappedCalculator と TimeProfileHook

```python
from pfcc_extras.calculators.wrapped_calculator import WrappedCalculator
from pfcc_extras.calculators.hooks.time_profile_hook import TimeProfileHook

wrapped = WrappedCalculator(base_calculator)
hook = TimeProfileHook()
wrapped.register_hook(hook)

atoms.calc = wrapped
atoms.get_potential_energy()

# 計算ステップの所要時間を集計（report() メソッドは存在しない）
print(f"呼び出し回数: {hook.n_call_count()}")
print(f"合計時間: {hook.total_time():.4f} s")
print(f"平均時間: {hook.avg_time():.4f} s")
```

`reset()` でフックの計測値をリセットできます。`WrappedCalculator` 側では `reset_hook()` で全フックをクリアしてから再登録してください。

---

## アルゴリズムユーティリティ

### farthest_point_sampling（最遠点サンプリング）

Light-PFP 学習データの多様性確保に使います。距離行列から互いに最も遠い点を貪欲法で選択します。

```python
from pfcc_extras.algorithm.sampling import farthest_point_sampling

selected_indices = farthest_point_sampling(
    distance_matrix=dist_mat,   # (N, N) の距離行列
    min_distance=0.2,
    seed=42,
)
selected_atoms = [all_atoms[i] for i in selected_indices]
```

### opt / opt_cell_size（段階収束型最適化）

収束を監視しながら段階的に `maxstep` を縮小するラッパー最適化です。

```python
from pfcc_extras.algorithm.optimization import opt, opt_cell_size

atoms = opt(atoms, calculator, sn=10, constraintatoms=[0, 1, 2])
atoms = opt_cell_size(atoms, calculator, sn=10, fix_symmetry=True)
```

- `sn`: 最大段階数（各段階で `maxstep` を `0.9^s` 倍に縮小）
- `opt_cell_size` は `FixSymmetry` で対称性を保ちながらセル最適化します

### opt_with_symmetry（UnitCellFilter + 対称性拘束）

`UnitCellFilter` ベースのシンプルなセル最適化です。`opt_cell_size` と異なり、単一ステップで収束まで実行します。

```python
from pfcc_extras.algorithm.optimization import opt_with_symmetry

atoms = opt_with_symmetry(
    atoms_in=atoms,
    calculator=calculator,
    fix_symmetry=True,          # FixSymmetry で対称性を保持
    hydrostatic_strain=False,   # True にすると等方圧縮のみ許可
    fmax=0.001,
)
```

---

## よくあるエラーと対処

| エラー / 問題 | 原因 | 対処 |
|---|---|---|
| `run_jobs` で `FAILED` になる | Notebook 内のエラー | `JobResult` の詳細を確認し、Notebook 単体でデバッグしてください |
| `SKIPPED` と `FAILED` の混同 | 出力ファイルが既に存在 | `SKIPPED` は正常動作（冪等性）です。`FAILED` とは区別して処理してください |
| `DepositionScheduler` で分子が重なる | `tol` が小さすぎる | `calculate_default_distance_tol(molecules)` で適切な距離閾値を算出してください |
| `ApplyUniformEfield` でエラー | `charges` の長さが原子数と不一致 | `len(charges) == len(atoms)` を確認してください |
| 部分占有の離散化誤差が大きい | `cell_repeat` が小さい | `cell_repeat` を増やして誤差を低減してください |
| 軌跡変換でフレーム数が合わない | 原子数・元素順の不一致 | 全フレームで原子数と元素順が一致しているか確認してください |
| 分子が PBC 境界で分断される | 原子単位の wrap | `wrap_molecule` で分子単位に wrap してください |
| `DeprecationWarning: shift_energy has been removed` | pfcc_extras 0.8.1 以降 | `Estimator().get_shift_energy_table()` を使用してください |
| MC の `mc_type_ratios` が合計 1.0 にならない | 比率の設定ミス | 各値の合計が 1.0 になるよう正規化してください |
| `IsosurfaceCalculator` で `grid_size` エラー | グリッドが格子定数より大きい | `grid_size` を最小格子定数より小さく設定してください |
| `cut_sphere` で中心が見つからない | `radius` がセルに対して大きすぎる | `radius` を小さくするかスーパーセルを大きくしてください |
| `WrappedCalculator` でフックが重複登録される | `register_hook` の重複呼び出し | `reset_hook()` でクリアしてから再登録してください |

## 関連ガイド

- **統合ワークフロー** (`integrated-workflows/SKILL.md`): pfcc_extras を使ったエンドツーエンドのパイプライン構築
- **バックグラウンドジョブ** (`background-job/SKILL.md`): 長時間ジョブの実行管理
- **Light-PFP** (`light-pfp/SKILL.md`): 大規模系での高速推論との組み合わせ
