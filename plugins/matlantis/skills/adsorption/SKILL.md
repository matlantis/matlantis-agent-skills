---
name: mt-adsorption
description: >
  表面への分子吸着を扱うスキルです。
  吸着エネルギー, adsorption energy, 吸着構造, adsorption structure,
  吸着サイト探索, adsorption site search, add_adsorbate,
  ontop, bridge, fcc, hcp, 吸着サイト, slab+molecule,
  FixAtoms, 下層固定, D3分散補正, pbe_plus_d3,
  Optuna, ブラックボックス最適化, 被覆率, coverage
  に関するコード生成時に使用してください。
---
# 吸着構造と吸着エネルギー

## 概要

表面 (slab) への分子吸着を扱うガイドです。吸着構造の構築、吸着エネルギーの計算、安定吸着サイトの探索までを包括的に扱います。

吸着エネルギーは以下で定義されます。

```
E_ads = E_slab+mol - (E_slab + E_mol)
```

- `E_slab+mol`: 吸着後の系（表面 + 吸着分子）のエネルギー
- `E_slab`: 清浄表面のエネルギー
- `E_mol`: 孤立分子のエネルギー

`E_ads` が負であるほど吸着が安定です。3 つのエネルギーは必ず同一の `calc_mode` / `model_version` / 拘束条件で計算してください。

本スキルのスコープは**単一系**の吸着構造構築・吸着エネルギー計算・サイト探索です。多数の候補構造を一括処理する規模のバッチスクリーニング（早期棄却フィルタ、CSV ランキング出力）は、統合ワークフロー (integrated-workflows/SKILL.md) のパターン B「吸着キャンペーンパイプライン」を参照してください。

## ワークフロー

吸着計算の標準的なワークフローは以下の通りです。

1. **バルクの最適化**: 格子定数を最適化する（`StrainFilter` によるセル緩和）
2. **表面 (slab) の構築**: `fcc111` 等で slab を作成し、十分な真空層を確保する
3. **分子の最適化**: 吸着させる分子を単独で構造最適化する
4. **吸着構造の構築**: `add_adsorbate` で分子を表面に配置する
5. **吸着系の最適化**: slab + 分子の系を構造最適化する
6. **吸着エネルギーの算出**: `E_ads = E_slab+mol - (E_slab + E_mol)` を計算する
7. **サイト探索**: 複数の配置候補を比較し、最安定構造を特定する

### サイト探索のパターン選択ガイド

| 状況 | 推奨パターン |
| :--- | :--- |
| 対称性の高い表面で既知の数サイト（ontop/bridge/fcc/hcp）を比較するだけ | パターン C（単純比較） |
| 配置位置・分子の向きまで含めて安定構造を探索する（一般のサイト探索） | **パターン D（Optuna、推奨）** |
| 多数の候補を一括スクリーニングし CSV でランキングする | integrated-workflows パターン B |

## 実装パターン

### パターン A: 吸着構造の構築（add_adsorbate）

`ase.build.add_adsorbate` で分子を表面に配置します。ASE の slab builder（`fcc111` 等）で作成した表面では、名前付きサイトを指定できます。

- `ontop`: 1 つの原子の直上
- `bridge`: 2 つの原子の中間
- `fcc`: 3 つの原子の間で、2 層目に原子が存在しない位置
- `hcp`: 3 つの原子の間で、2 層目の原子の直上の位置

```python
from ase.build import fcc111, molecule, add_adsorbate
from ase.constraints import FixAtoms

# Pt(111) 表面の作成（4x4 ユニットセル、4 層、真空 40A）
slab = fcc111("Pt", size=(4, 4, 4), vacuum=40.0, periodic=True)

# 下層の固定: バルク的な層を FixAtoms で固定する
z_threshold = slab.positions[:, 2].min() + 2.0
slab.set_constraint(FixAtoms(mask=slab.positions[:, 2] < z_threshold))

# CO 分子を fcc サイトに配置（height は表面からの初期高さ [A]）
mol = molecule("CO")
add_adsorbate(slab, mol, height=2.0, position="fcc")

# 名前付きサイトの代わりに (x, y) 座標 [A] での明示指定も可能
# add_adsorbate(slab, mol, height=2.0, position=(3.2, 1.8))
```

小分子の吸着構造での初期高さ `height` は 1.5〜3.0 A 程度が目安です。低すぎると表面にめり込み、高すぎると最適化中に脱離したまま収束します。

### パターン B: 吸着エネルギー計算

3 つのエネルギー（slab+mol, slab, mol）を同一条件(同一のversion、 calc_mode、収束条件 (`fmax`))で計算し、吸着エネルギーを算出します。吸着系では分子-表面間の分散力が重要であるため、pfpバージョンによらずD3分散補正付きのcalc_mode(`calc_mode="r2scan_plus_d3"`、`calc_mode="pbe_u_plus_d3"`もしくは`calc_mode="pbe_plus_d3"`)を使用してください。

```python
from ase import Atoms
from ase.build import bulk, fcc111, molecule, add_adsorbate
from ase.filters import StrainFilter, FrechetCellFilter
from ase.optimize import LBFGS
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

def get_opt_energy(
    atoms: Atoms,
    calculator: ASECalculator,
    fmax: float = 0.001,
    opt_mode: str = "normal"
) -> float:
    """構造を局所最適化してエネルギーを返す。

    Args:
        atoms: 対象構造
        calculator: PFP calculator
        fmax: 力の収束判定値 [eV/A]
        opt_mode: "normal"(原子位置のみ) / "scale"(格子定数) / "all"(セル+原子)
    Returns:
        最適化後のポテンシャルエネルギー [eV]
    """
    atoms.calc = calculator
    if opt_mode == "scale":
        opt = LBFGS(StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0]), logfile=None)
    elif opt_mode == "all":
        # ExpCellFilter は deprecated（FrechetCellFilter が推奨）
        atoms_filter = FrechetCellFilter(atoms)
        opt = LBFGS(atoms_filter, logfile=None)
    else:
        opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    return atoms.get_potential_energy()

def calc_adsorption_energy(
    slab: Atoms,
    mol: Atoms,
    site: str = "fcc",
    height: float = 2.0,
    model_version: str = "v9.0.0",
    calc_mode: str = "r2scan_plus_d3"
) -> float:
    """吸着エネルギー E_ads = E_slab+mol - (E_slab + E_mol) を計算する。"""
    calculator = ASECalculator(
        Estimator(model_version=model_version, calc_mode=calc_mode)
    )

    # 清浄表面と孤立分子のエネルギー（同一 calculator で統一）
    slab_opt = slab.copy()
    E_slab = get_opt_energy(slab_opt, calculator)

    mol_opt = mol.copy()
    E_mol = get_opt_energy(mol_opt, calculator)

    # 吸着構造の構築と最適化
    system = slab_opt.copy()
    add_adsorbate(system, mol_opt, height=height, position=site)
    E_slab_mol = get_opt_energy(system, calculator)

    E_ads = E_slab_mol - (E_slab + E_mol)
    print(f"E_slab+mol = {E_slab_mol:.3f} eV")
    print(f"E_slab     = {E_slab:.3f} eV")
    print(f"E_mol      = {E_mol:.3f} eV")
    print(f"E_ads      = {E_ads:.3f} eV")
    return E_ads
```

`E_slab` の計算と `E_slab+mol` の計算で `FixAtoms` 等の拘束条件を一致させてください。片方だけ下層を固定すると、拘束の差分がそのまま吸着エネルギーの誤差になります。

### パターン C: 吸着サイトの単純比較

既知の数サイトを総当たりで比較する簡易法です。対称性の高い表面で候補が自明な場合に使用します。配置・配向まで含めた本格的な探索にはパターン D（Optuna）を推奨します。

```python
from typing import List, Tuple
from ase import Atoms
from ase.build import add_adsorbate
from ase.optimize import LBFGS
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

def search_adsorption_sites(
    slab: Atoms,
    adsorbate: Atoms,
    sites: List[str] = ["ontop", "bridge", "fcc", "hcp"],
    height: float = 2.0,
    fmax: float = 0.05,
    model_version: str = "v9.0.0",
    calc_mode: str = "r2scan_plus_d3"
) -> Tuple[List[dict], dict]:
    """候補サイトごとに吸着構造を最適化し、エネルギーを比較する。

    Returns:
        results: 各候補サイトのエネルギーと構造のリスト
        best: 最小エネルギー構造の情報
    """
    results = []
    for site in sites:
        # 各候補サイトに adsorbate を配置
        system = slab.copy()
        add_adsorbate(system, adsorbate.copy(), height=height, position=site)

        # PFP で構造最適化（候補ごとに calculator を生成）
        estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
        system.calc = ASECalculator(estimator)
        opt = LBFGS(system, logfile=None)
        opt.run(fmax=fmax)

        results.append({
            'site': site,
            'atoms': system,
            'energy': system.get_potential_energy()
        })

    # 最小エネルギー構造を特定
    best = min(results, key=lambda x: x['energy'])
    print(f"Best site: {best['site']} (E = {best['energy']:.3f} eV)")
    return results, best
```

名前付きサイトの代わりに `(x, y)` 座標のリストを総当たりすることもできますが、候補数が増える場合はパターン D の方が少ない試行数で安定構造に到達できます。

### パターン D: Optuna によるブラックボックス探索（推奨）

吸着サイト探索の**推奨パターン**です。[Optuna](https://optuna.org/) によるブラックボックス最適化で、吸着位置 (x, y, height) と分子の配向（オイラー角 phi, theta, psi）を同時に探索します。グリッド総当たりと異なり、TPE サンプラーが有望な領域を集中的に探索するため、少ない試行数で安定構造に到達できます。

slab と分子は `Study.user_attrs` に JSON 文字列として保存し、トライアルごとの再構築・再最適化を回避します。

```python
import io
import numpy as np
import optuna
from ase import Atoms
from ase.build import add_adsorbate
from ase.io import read, write

def atoms_to_json(atoms: Atoms) -> str:
    """Atoms を JSON 文字列に変換する（Study/Trial への保存用）。"""
    f = io.StringIO()
    write(f, atoms, format="json")
    return f.getvalue()

def json_to_atoms(atoms_str: str) -> Atoms:
    """JSON 文字列から Atoms を復元する。"""
    return read(io.StringIO(atoms_str), format="json")

def objective(trial: optuna.Trial) -> float:
    """吸着エネルギーを目的関数としたブラックボックス最適化。"""
    # Study に保存した最適化済み構造を復元
    slab = json_to_atoms(trial.study.user_attrs["slab"])
    E_slab = trial.study.user_attrs["E_slab"]
    mol = json_to_atoms(trial.study.user_attrs["mol"])
    E_mol = trial.study.user_attrs["E_mol"]

    # 分子の配向（オイラー角）: theta は縮退回避のため 0-180 度に変換
    phi = 180.0 * trial.suggest_float("phi", -1, 1)
    theta = np.arccos(trial.suggest_float("theta", -1, 1)) * 180.0 / np.pi
    psi = 180.0 * trial.suggest_float("psi", -1, 1)
    mol.euler_rotate(phi=phi, theta=theta, psi=psi)

    # 吸着位置: セルの分数座標で指定し、対称性の範囲 (0-0.5) に制限
    x_pos = trial.suggest_float("x_pos", 0, 0.5)
    y_pos = trial.suggest_float("y_pos", 0, 0.5)
    z_hig = trial.suggest_float("z_hig", 1, 5)
    xy_position = np.matmul([x_pos, y_pos, 0], slab.cell)[:2]

    add_adsorbate(slab, mol, z_hig, xy_position)
    E_slab_mol = get_opt_energy(slab, calculator, fmax=1e-3)

    # 最適化後の構造を後から参照できるように保存
    trial.set_user_attr("structure", atoms_to_json(slab))

    return E_slab_mol - E_slab - E_mol

# slab / mol は事前に最適化して Study に保存（トライアルごとの再計算を回避）
study = optuna.create_study()
study.set_user_attr("slab", atoms_to_json(slab))
study.set_user_attr("E_slab", E_slab)
study.set_user_attr("mol", atoms_to_json(mol))
study.set_user_attr("E_mol", E_mol)

study.optimize(objective, n_trials=30)

print(f"Best trial: #{study.best_trial.number}")
print(f"  E_ads = {study.best_value:.3f} eV")
print(f"  params = {study.best_params}")

# 最安定構造の復元
best_atoms = json_to_atoms(study.best_trial.user_attrs["structure"])
```

#### 簡易版: 名前付きサイトのカテゴリカル探索

探索対象を既知の 4 サイトに限定する場合は、カテゴリカル変数として探索できます。

```python
def objective(trial: optuna.Trial) -> float:
    slab = json_to_atoms(trial.study.user_attrs["slab"])
    mol = json_to_atoms(trial.study.user_attrs["mol"])

    site = trial.suggest_categorical("site", ["ontop", "bridge", "fcc", "hcp"])
    add_adsorbate(slab, mol, 3.0, site)
    E_slab_mol = get_opt_energy(slab, calculator, fmax=1e-3)

    return E_slab_mol - trial.study.user_attrs["E_slab"] - trial.study.user_attrs["E_mol"]
```

#### 補足

- `optuna.visualization.plot_optimization_history(study)` で探索過程を可視化できます
- Study のストレージはデフォルトでメモリ上のため、永続化が必要な場合は `optuna.create_study(storage=...)` で RDB やファイルを指定してください
- フルノートブックは `matlantis/references/community-examples/adsorption_structure_search_basic/` を参照してください

## ベストプラクティス

### サイト探索は Optuna を第一候補にする

配置・配向の自由度を含む一般のサイト探索では、グリッド総当たりより Optuna（パターン D）の方が少ない試行数で安定構造に到達できます。単純比較（パターン C）は候補が数個で自明な場合に限定してください。

### D3 分散補正を使用する

吸着では分子-表面間の分散力（ファンデルワールス力）が重要です。pfpバージョンによらずD3分散補正付きのcalc_mode(`calc_mode="r2scan_plus_d3"`、`calc_mode="pbe_u_plus_d3"`もしくは`calc_mode="pbe_plus_d3"`)を使用してください。R2SCANやPBEなどの分散力補正なしでは吸着エネルギーを過小評価する傾向があります。

### 3 つのエネルギーの計算条件を統一する

`E_slab+mol`, `E_slab`, `E_mol` は必ず同一の `calc_mode` / `model_version` / 収束条件 (`fmax`) で計算してください。条件が混在すると吸着エネルギーが誤差を持ちます。

### 拘束条件を一致させる

下層を `FixAtoms` で固定する場合、`E_slab` と `E_slab+mol` の両方で同じ拘束を適用してください。

### 十分な真空層を確保する

周期境界条件下では、真空層が薄いと slab の上下イメージ間で吸着分子が相互作用します。真空層は最低でも 10 A、吸着分子が大きい場合はさらに厚く設定してください。

### NEB の端点に使う場合は厳密に収束させる

吸着構造を反応経路探索（NEB / ReactionStringFeature）の端点として使用する場合は、`fmax < 0.01` eV/A まで最適化を完了させてください（reaction/SKILL.md 参照）。

### 多数候補のバッチ処理は統合ワークフローへ

候補を大量に生成して早期棄却 → 最適化 → CSV ランキングする場合は、integrated-workflows/SKILL.md のパターン B（吸着キャンペーンパイプライン）を使用してください。

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| E_ads が正の値になる | サイトが不安定 / 最適化未収束 / 参照エネルギーの条件不一致 | 他のサイトを試す。fmax 収束を確認する。3 項の calc_mode / 拘束条件を統一する |
| 吸着分子が表面から脱離する | 初期 height が高すぎる | height を 1.5〜3.0 A 程度に下げる |
| 吸着分子が表面にめり込む | 初期 height が低すぎる | height を上げる。構造最適化初期のForceが大きい場合は配置を見直す |
| サイト間の比較結果が不自然 | calc_mode / model_version の混在 | 全候補を同一条件で計算する |
| 吸着エネルギーが拘束に依存して変わる | E_slab と E_slab+mol で FixAtoms 設定が異なる | 両者で同一の拘束を適用する |
| `KeyError` (named site) | ASE builder 以外で作成した slab に名前付きサイトを指定した | `(x, y)` 座標で明示指定するか、ASE builder で slab を作成する |
| `MultiCalculatorUseDetected` | 複数候補の並列計算で calculator を共有した | 候補ごとに Estimator + Calculator を生成する |
| Optuna が安定構造に到達しない | 試行数不足 / 探索範囲が広すぎる | `n_trials` を増やす。対称性を利用して探索範囲を狭める（分数座標 0-0.5 等） |

## 関連ガイド

- **構造モデリング** (modeling/SKILL.md): slab 構築の詳細（pymatgen による任意ミラー指数の切り出し等）
- **構造最適化** (optimization/SKILL.md): LBFGS / FIRE の使い分け、セル緩和
- **反応経路探索** (reaction/SKILL.md): 吸着構造を端点とした NEB / 遷移状態探索
- **統合ワークフロー** (integrated-workflows/SKILL.md): 吸着キャンペーンパイプライン（バッチスクリーニング）
- **PFP** (pfp/SKILL.md): calc_mode の選択指針（D3 分散補正を含む）
