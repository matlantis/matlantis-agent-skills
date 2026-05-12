---
name: mt-reaction
description: >
  反応経路探索・遷移状態探索を扱うスキルです。
  NEB, Nudged Elastic Band, NEBFeature, NEBTools, climbing image, CI-NEB,
  ReactionStringFeature, fmax_rp, fmax_ts, fmax_rd, fmax_eq,
  RestScanFeature, ScanDistance, ScanDihedral,
  modify_trajectory_by_mic, 原子マッピング, atom mapping,
  反応経路, reaction path, 遷移状態, transition state, 活性化エネルギー, activation energy,
  エネルギー障壁, energy barrier, IDPP, allow_shared_calculator,
  吸着サイト探索, adsorption site search
  に関するコード生成時に使用してください。
---

# 反応経路探索

## 概要

化学反応の反応経路（Reaction Path）と遷移状態（Transition State, TS）、および活性化エネルギー（Activation Energy, Ea）を探索するためのガイドです。Matlantis では、ASE 標準の NEB（Nudged Elastic Band）法と、より高機能で自動化された `ReactionStringFeature` の両方を利用できます。

NEB 法は反応経路上の中間構造（image）をバンドで接続し、最小エネルギー経路（MEP: Minimum Energy Path）を求める手法です。Climbing Image NEB（CI-NEB）を有効にすることで、バンド上の最高エネルギー image が正確な遷移状態に近づきます。

`ReactionStringFeature` は matlantis-features が提供する高レベル API であり、NEB 法をベースにしつつ自動パラメータ調整やタイムアウト管理を備えています。網羅的な全反応探索が必要な場合は、GRRM 連携ワークフローを別途検討してください。

本ガイドでは、端点の最適化から NEB / ReactionString の実行、エネルギー障壁の解析、遷移状態の検証までを包括的に扱います。

## ワークフロー

反応経路探索の標準的なワークフローは以下の通りです。

1. **端点の最適化**: 始状態（Reactant）と終状態（Product）を構造最適化する（`fmax < 0.05` eV/A 必須）
2. **原子対応の確認**: 始状態と終状態で原子インデックスの対応関係を検証する
3. **構造の位置合わせ**: 周期境界条件（PBC）を考慮した連続的な軌跡補正を行う
4. **中間 image の生成**: 線形補間または IDPP 補間で中間構造を作成する
5. **経路探索の実行**: NEB 法または ReactionStringFeature を実行する
6. **エネルギー障壁の解析**: 活性化エネルギー（Ea）と反応エネルギー（dE）を算出する
7. **遷移状態の検証**: 必要に応じて虚振動解析を行う

各ステップを省略すると計算が破綻する可能性があるため、必ず順序通りに実行してください。

## 実装パターン

### パターン A: ReactionStringFeature による経路探索

matlantis-features の `ReactionStringFeature` を用いた高レベル API による実装です。タイムアウト管理やバッチ処理に対応しています。

#### PBC 補正関数

周期境界条件をまたぐ原子移動がある場合、補間が破綻するため、必ず `modify_trajectory_by_mic` による補正を行います。これは FAQ 記載の推奨実装です。

```python
from ase.geometry import find_mic
from ase import Atoms
from typing import List
import numpy as np

def modify_trajectory_by_mic(traj: List[Atoms]) -> List[Atoms]:
    """
    周期境界条件(PBC)を考慮して、軌跡が連続的になるように原子位置を補正する。
    端点最適化後に wrap() で座標を折り返すと補間が破綻しやすいため、
    この関数で連続化してから中間像を生成する。
    """
    output: List[Atoms] = []
    previous_pos: np.ndarray = None

    for i in range(len(traj)):
        if i == 0:
            output.append(traj[i].copy())
            previous_pos = traj[i].get_positions()
        else:
            pos_diff = find_mic(
                traj[i].get_positions() - previous_pos,
                traj[i].get_cell(),
                pbc=True
            )
            current_atoms = traj[i].copy()
            shift = pos_diff if isinstance(pos_diff, tuple) else pos_diff
            current_atoms.set_positions(previous_pos + shift)
            output.append(current_atoms)
            previous_pos = current_atoms.get_positions()

    return output
```

#### ReactionStringFeature の実行

```python
from matlantis_features.features.reaction import ReactionStringFeature
from matlantis_features.utils.calculators import pfp_estimator_fn
from datetime import timedelta

def run_reaction_string(
    start_atoms: Atoms,
    end_atoms: Atoms,
    calc_mode: str = "PBE",
    timeout_sec: int = 1800
) -> List[Atoms]:
    """
    ReactionStringFeatureを用いて反応経路を探索する。
    """
    # 入力トラジェクトリの作成とPBC補正
    traj_input = [start_atoms, end_atoms]
    traj_use = modify_trajectory_by_mic(traj_input)

    # ReactionStringFeature の設定
    # fmax_rp: 反応経路の力の収束判定値
    # fmax_ts: 遷移状態における力の収束判定値
    # fmax_rd: 律速段階の力の収束判定値
    # fmax_eq: 平衡構造の力の収束判定値
    rp = ReactionStringFeature(
        fmax_rp=0.05,
        fmax_ts=0.01,
        fmax_rd=0.005,
        fmax_eq=0.05,
        timeout=timedelta(seconds=timeout_sec),
        estimator_fn=pfp_estimator_fn(
            model_version="v8.0.0",
            calc_mode=calc_mode
        )
    )

    # 入力はリストのリスト（バッチ処理対応）
    rp_result = rp([traj_use])

    # 結果の展開
    images = rp_result.reaction_string_images
    pathway_atoms = []
    if images and len(images) > 0:
        pathway_atoms = [img.ase_atoms for img in images]

    return pathway_atoms
```

#### バッチ処理

複数の反応経路を一括で探索する場合は、リストのリストとして入力します。

```python
# 複数反応のバッチ実行
result = rs([
    [initial1, final1],
    [initial2, final2],
])
```

### パターン B: ASE NEB 法による経路探索

ASE 標準の NEB 法を直接利用する実装です。パラメータの細かな制御が可能です。

#### NEB の構築と実行

```python
from ase.mep import NEB
from ase.optimize import FIRE
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

def run_neb(
    react: Atoms,
    prod: Atoms,
    n_images: int = 7,
    k: float = 0.1,
    fmax: float = 0.05,
    steps: int = 500,
    calc_mode: str = "PBE"
) -> list:
    """
    ASE NEB法による反応経路探索。
    各imageに個別のCalculatorを設定する（共有禁止）。
    """
    # image リストの構築
    images = [react.copy()]
    images += [react.copy() for _ in range(n_images)]
    images += [prod.copy()]

    # 各imageに個別のCalculatorを設定
    # allow_shared_calculator=False の理由:
    #   NEB像間でcalculator を共有すると、
    #   像ごとの独立計算ができず競合や不整合が発生する
    for image in images:
        image.calc = ASECalculator(
            Estimator(
                model_version="v8.0.0",
                calc_mode=calc_mode
            )
        )

    # NEB オブジェクトの構築
    neb = NEB(
        images,
        k=k,                              # バネ定数
        climb=True,                        # CI-NEB で正確な TS を取得
        allow_shared_calculator=False,     # Calculator 共有を禁止
        parallel=True                      # 並列計算（個別Calculator必須）
    )
    neb.interpolate()  # 線形補間で初期経路を生成

    # 最適化の実行
    opt = FIRE(neb)
    opt.run(fmax=fmax, steps=steps)

    return images
```

#### Calculator 分離の原則

`parallel=True` を使用する場合は、各 image に個別の Calculator を設定する必要があります。1 つの Estimator を複数の ASECalculator で共有すると `MultiCalculatorUseDetected` エラーが発生します。

```python
# 正しいパターン: 各imageに個別のEstimator + Calculator
for image in images:
    image.calc = ASECalculator(
        Estimator(model_version="v8.0.0", calc_mode=calc_mode)
    )

# 誤ったパターン: Estimatorの共有（エラーになる）
# estimator = Estimator(...)
# for image in images:
#     image.calc = ASECalculator(estimator)  # NG
```

### パターン C: matlantis-features NEBFeature による経路探索

matlantis-features の `NEBFeature` は NEB 法の高レベルラッパーです。

```python
from matlantis_features.features.common.neb import NEBFeature
from matlantis_features.features.common.opt import FireASEOptFeature
from matlantis_features.utils.calculators import pfp_estimator_fn

neb = NEBFeature(
    optimizer=FireASEOptFeature(),
    n_images=5,
    climb=True,        # CI-NEB で正確な TS を取得
    idpp=False,        # 結晶系ではIDPP を無効にする場合がある
    estimator_fn=pfp_estimator_fn(
        model_version="v8.0.0",
        calc_mode="PBE"
    ),
)
neb_result = neb(initial_atoms, final_atoms)
```

### パターン D: 拘束スキャン（RestScanFeature）

反応座標に沿ったポテンシャルエネルギー表面のスキャンを行います。

```python
from matlantis_features.features.reaction import RestScanFeature, ScanDistance, ScanDihedral

# 距離スキャン
scan = RestScanFeature(
    scan_params=[ScanDistance(atom_indices=[0, 1], values=[1.0, 1.5, 2.0, 2.5, 3.0])],
    estimator_fn=estimator_fn,
)

# 二面角スキャン
scan = RestScanFeature(
    scan_params=[ScanDihedral(atom_indices=[0, 1, 2, 3], values=[0, 30, 60, 90, 120, 150, 180])],
    estimator_fn=estimator_fn,
)
```

### パターン E: エネルギー障壁解析

探索された経路のエネルギープロファイルを計算し、活性化エネルギーを算出します。

```python
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator
import matplotlib.pyplot as plt
import numpy as np

def analyze_reaction_profile(pathway: List[Atoms], calc_mode: str = "PBE"):
    """
    反応経路のエネルギープロファイルを計算・表示する。
    Forward/backward barrierの両方を出力する。
    """
    estimator = Estimator(calc_mode=calc_mode, model_version="v8.0.0")
    calc = ASECalculator(estimator)

    energies = []
    distances = [0.0]

    for i, atoms in enumerate(pathway):
        atoms.calc = calc
        e = atoms.get_potential_energy()
        energies.append(e)
        if i > 0:
            d = np.linalg.norm(
                atoms.get_positions() - pathway[i-1].get_positions()
            )
            distances.append(distances[-1] + d)

    # 始状態基準の相対エネルギー
    base_energy = energies[0]
    rel_energies = [e - base_energy for e in energies]

    # Forward barrier: 始状態からTSまでの障壁
    forward_barrier = max(rel_energies)
    # Backward barrier: 終状態からTSまでの障壁
    backward_barrier = max(rel_energies) - rel_energies[-1]
    # 反応エネルギー
    reaction_energy = rel_energies[-1]

    print(f"Forward Activation Energy (Ea_fwd):  {forward_barrier:.3f} eV")
    print(f"Backward Activation Energy (Ea_bwd): {backward_barrier:.3f} eV")
    print(f"Reaction Energy (dE):                {reaction_energy:.3f} eV")

    plt.plot(distances, rel_energies, 'o-')
    plt.xlabel("Reaction Coordinate (A)")
    plt.ylabel("Relative Energy (eV)")
    plt.title(f"Reaction Profile (Ea={forward_barrier:.2f} eV)")
    plt.grid(True)
    plt.show()

    return rel_energies
```

### パターン F: NEBTools による収束評価

`NEBTools` を使用してバリアと最大力を同時に評価し、未収束の誤判定を防ぎます。

```python
from ase.mep import NEBTools

def evaluate_neb_convergence(images):
    """
    NEBToolsを用いてバリアとfmaxの二軸で収束を評価する。
    """
    tools = NEBTools(images)
    barrier = tools.get_barrier()
    fmax = tools.get_fmax()

    print(f"Energy Barrier: {barrier[0]:.3f} eV (forward), {barrier[1]:.3f} eV (backward)")
    print(f"Maximum Force:  {fmax:.4f} eV/A")

    if fmax > 0.05:
        print("WARNING: NEB has not converged (fmax > 0.05)")
    else:
        print("NEB converged successfully")

    return barrier, fmax
```

### パターン G: 吸着サイト探索

表面上の吸着サイトを網羅的に探索する方法です。

```python
def search_adsorption_sites(
    slab: Atoms,
    adsorbate: Atoms,
    positions: list,
    calc_mode: str = "PBE"
):
    """
    パラメータ化した初期配置を複数作り、
    各候補をPFPで最適化してエネルギーを比較する。
    """
    results = []
    for pos in positions:
        # 各候補位置にadsorbateを配置
        system = slab.copy()
        ads = adsorbate.copy()
        ads.translate(pos)
        system += ads

        # PFPで構造最適化
        system.calc = ASECalculator(
            Estimator(model_version="v8.0.0", calc_mode=calc_mode)
        )
        # 最適化後のエネルギーを記録
        results.append({
            'position': pos,
            'atoms': system,
            'energy': system.get_potential_energy()
        })

    # 最小エネルギー構造を特定
    best = min(results, key=lambda x: x['energy'])
    return results, best
```

## ベストプラクティス

### 端点の最適化は必須

始状態と終状態は必ず事前に構造最適化（`fmax < 0.05` eV/A）を完了させてください。最適化されていない構造をつなぐと経路が不安定になり、非物理的な結果を得る原因になります。

### Atom Mapping の検証

始状態と終状態で原子のインデックス順序が完全に一致している必要があります。原子番号 1 番が始状態で「炭素 A」なら、終状態でも「炭素 A」でなければなりません。ここがずれていると、原子が物理的にありえない移動を行い、エネルギーが発散します。

### PBC 補正の徹底

周期境界条件をまたぐ原子移動がある場合、`modify_trajectory_by_mic` による補正を必ず行ってください。「原子が突然消えた」「結合が異常に伸びた」ように見える場合は、PBC の補正ミスが原因です。

### 補間法の使い分け

| 補間法 | 特徴 | 推奨用途 |
| :--- | :--- | :--- |
| **Linear** | 単純な直線補間。高速だが初期経路が粗くなりやすい | 端点が近く単純移動の場合 |
| **IDPP** | Image Dependent Pair Potential。原子間距離を考慮した補間 | 初期経路の粗さが強い場合 |

### 最適化器の選択

NEB 本計算では `FIRE` を使うと、バンド全体の力が大きい初期段階で LBFGS より安定するケースがあります。

### 終了判定の二軸化

`NEBTools.get_barrier()` だけでなく `get_fmax()` も同時評価し、「見かけ上バリアが出たが未収束」の誤判定を防いでください。

### 虚振動解析による TS 確認

厳密な遷移状態確認には、得られた TS 構造に対して振動解析を行い、虚の振動数が 1 つだけ存在することを確認するのが一般的です。ただし計算コストが非常に高いため、まずはエネルギープロファイルの形状で判断し、必要に応じて実施してください。

```python
from matlantis_features.features.common.vibration import VibrationFeature

vib = VibrationFeature(delta=0.01, estimator_fn=estimator_fn)
vib_result = vib(ts_atoms)
# 虚振動数が1つだけであることを確認
```

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| エネルギーが発散する | Atom Mapping のずれ（始状態と終状態で原子インデックスの対応が不一致） | 始状態と終状態の原子順序を確認し、必ず同じ原子が同じインデックスに対応するように修正する |
| 原子が突然消える、結合が異常に伸びる | PBC の補正ミス | `modify_trajectory_by_mic` を必ず適用する。`wrap()` を補間前に使わない |
| `MultiCalculatorUseDetected` | 1 つの Estimator を複数の Calculator で共有した | 各 image に個別の Estimator + Calculator を生成する |
| NEB が収束しない | 初期経路が不適切、または steps 不足 | IDPP 補間を試す、image 数を増やす、steps を増やす |
| Timeout で計算が中断される | ReactionStringFeature のデフォルトタイムアウト超過 | `timeout` 引数で時間を延長する、中間構造を増やして初期パスを改善する |
| 「見かけ上バリアが出たが未収束」 | `get_barrier()` のみで判定し `get_fmax()` を確認していない | `NEBTools.get_fmax()` も同時に評価し、fmax < 0.05 を確認する |
| 非物理的な経路が得られる | 端点が未最適化（fmax > 0.05） | 端点を `fmax < 0.05` まで最適化してから NEB を実行する |
| `ConcurrentUseDetected` | Calculator を複数スレッドで共有した | 各スレッドで新しく Calculator を生成する |

## 関連ガイド

- 構造最適化: 端点の最適化手順の詳細
- **物性計算** (property-analysis/SKILL.md): 振動解析による TS 検証
- **可視化** (visualization/SKILL.md): 反応経路の構造確認
- **外部データ連携** (external-data/SKILL.md): 反応物・生成物の構造取得
