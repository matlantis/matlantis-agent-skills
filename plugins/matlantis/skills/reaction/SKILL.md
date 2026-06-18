---
name: mt-reaction
description: >
  反応経路探索・遷移状態探索を扱うスキルです。
  NEB, Nudged Elastic Band, NEBFeature, NEBTools, climbing image, CI-NEB,
  ReactionStringFeature, fmax_rp, fmax_ts, fmax_rd, fmax_eq,
  RestScanFeature, ScanDistance, ScanDihedral,
  NEBFeatureResult, ReactionStringFeatureResult, RestScanFeatureResult,
  plot_energy, energy_log, get_images, num_segments, transition_state_indices,
  modify_trajectory_by_mic, 原子マッピング, atom mapping,
  反応経路, reaction path, 遷移状態, transition state, 活性化エネルギー, activation energy,
  エネルギー障壁, energy barrier, IDPP, allow_shared_calculator,
  Sella, TS最適化, TS optimization, 鞍点, saddle point, 虚振動, imaginary frequency,
  反応プロファイル, reaction profile
  に関するコード生成時に使用してください。
---

# 反応経路探索

## 概要

化学反応の反応経路（Reaction Path）と遷移状態（Transition State, TS）、および活性化エネルギー（Activation Energy, Ea）を探索するためのガイドです。

反応経路探索は **「準備 → 経路探索 → 解析」の 3 段**で構成されます。準備段階には 2 つのルート（a: 端点を手動準備して補間、b: RestScan で終構造と初期経路を生成）があり、経路探索には 3 つの手法（A: ASE NEB、B: NEBFeature、C: ReactionStringFeature）があります。**解析コードは使用した手法によって異なる**ため、本ガイドの解析セクションは手法別に分かれています。

NEB（Nudged Elastic Band）法は反応経路上の中間構造（image）をバンドで接続し、最小エネルギー経路（MEP: Minimum Energy Path）を求める手法です。ASE標準のNEBではclimb=TrueとすることでClimbing Image NEB（CI-NEB）を有効にし、バンド上の最高エネルギー image が正確な遷移状態に近づきます。`ReactionStringFeature` は matlantis-features が提供する高レベル API であり、String 法をベースにしつつ自動パラメータ調整やタイムアウト管理を備え、TS の最適化まで自動で行います。

NEB法とReactionStringFeatureはいずれも端点（始構造(Reactant) と終構造(Product)）の情報を必要とする、Double-Ended Methods に分類されます。端点の構造は事前に最適化されている必要があり、原子のインデックス順序も完全に一致している必要があります。

網羅的な全反応探索、終構造が不明な場合の反応経路自動探索が必要な場合は、GRRM 連携ワークフローを別途検討してください。

## ワークフロー

```
1. 準備段階（a / b のいずれかを選択）
   準備 a: 始構造・終構造を手動準備（または別手順で作成した初期構造を使用）
           初期推定経路は補間（linear / IDPP）で作成
   準備 b: 始構造のみ準備し、終構造と初期推定経路を RestScan で生成
   共通の前処理: 端点の最適化・原子対応の確認・PBC 補正
   ↓
2. 経路探索（A / B / C のいずれかを選択）
   手法 A: ASE の NEB
   手法 B: Matlantis Features の NEBFeature
   手法 C: Matlantis Features の ReactionStringFeature
   ↓
3. 解析（使用した手法に対応する解析コードを使用）
   A → NEBTools・反応プロファイル plot
   B → NEBFeatureResult
   C → ReactionStringFeatureResult
   ↓
4. 遷移状態の精密化と検証（必要に応じて）
   Sella による鞍点最適化・虚振動解析
```

### 手法選択ガイド

| 手法 | 特徴 | 向くケース |
| :--- | :--- | :--- |
| **A: ASE NEB** | パラメータを細かく制御できる ASE 標準実装 | バネ定数・optimizer・収束条件等を自分で調整したい場合 |
| **B: NEBFeature** | NEB の高レベル API。結果オブジェクトに plot / 出力ヘルパーを持つ | 定型的な NEB を簡潔に実行したい場合 |
| **C: ReactionStringFeature** | 自動パラメータ調整・TS 収束・多段反応のセグメント分割・タイムアウト管理 | 多段反応、TS まで自動で収束させたい場合、RestScan 経路を初期値にする場合 |

## 準備段階

### 共通の前処理

ルート a / b のいずれでも、経路探索に入る前に以下を必ず実施します。

**1. 端点の最適化**: 始状態・終状態を構造最適化します。最低でも `fmax < 0.01` eV/A、理想的には `fmax < 0.001` eV/A まで収束させてください。

```python
from ase.optimize import LBFGS
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

for atoms in (react, prod):
    atoms.calc = ASECalculator(
        Estimator(model_version="v9.0.0", calc_mode="R2SCAN_PLUS_D3")
    )
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=0.001)  # 理想は 0.001、最低でも 0.01
```

**2. 原子対応（Atom Mapping）の確認**: 始状態と終状態で原子のインデックス順序が完全に一致している必要があります。ずれていると原子が物理的にありえない移動を行い、エネルギーが発散します。

**3. PBC 補正**: 周期境界条件をまたぐ原子移動がある場合、線形補間や IDPP 補間が破綻するため、必ず `modify_trajectory_by_mic` による補正を行います。ASE の `find_mic`（minimum image convention）を用いて、隣接構造間の変位が最小になるように原子位置を連続化します。原典は pfcc-extras の `examples/reaction_path/modify_trajectory_by_mic.ipynb` です。

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
            # find_mic は (最小イメージ変位ベクトル, その長さ) のタプルを返す
            shift, _ = find_mic(
                traj[i].get_positions() - previous_pos,
                traj[i].get_cell(),
                pbc=True
            )
            current_atoms = traj[i].copy()
            current_atoms.set_positions(previous_pos + shift)
            output.append(current_atoms)
            previous_pos = current_atoms.get_positions()

    return output
```

```python
# 端点に適用してから経路探索に渡す
react, prod = modify_trajectory_by_mic([react, prod])
```

### 準備 a: 端点の手動準備と補間による初期経路

始構造・終構造を手動で構築するか、別手順で作成した初期構造（構造モデリング、外部データベースからの取得、吸着構造探索の結果など）を使用するルートです。端点を共通の前処理（最適化・atom mapping・PBC 補正）にかけたうえで、初期推定経路は補間で作成します。

補間の実施方法は経路探索の手法ごとに異なります。

| 手法 | 補間の指定方法 |
| :--- | :--- |
| A: ASE NEB | `neb.interpolate()`（線形）/ `neb.interpolate(method="idpp")` |
| B: NEBFeature | コンストラクタの `idpp` フラグ（False で線形） |
| C: ReactionStringFeature | 端点 2 点を渡すと内部で線形補間。中間構造を含む軌跡を渡すことも可能 |

### 準備 b: RestScan による終構造・初期推定経路の生成

終構造が不明だが反応座標（結合距離・二面角など）の見当がつく場合のルートです。`RestScanFeature` で反応座標に拘束力をかけて構造を段階的に変化させ、終構造の候補と初期推定経路を生成します。

注意: RestScan の経路は多くの場合、最適な反応経路（MEP）上を通らず、TS（一次の鞍点）も通りません。活性化障壁が必要な場合は、得られた大まかな経路を初期値として経路探索（特に手法 C）で最適化してください。

```python
from matlantis_features.features.reaction.rest_scan import RestScanFeature, ScanDistance, ScanDihedral
from matlantis_features.features.common.opt import LBFGSASEOptFeature
from matlantis_features.utils.calculators import pfp_estimator_fn

estimator_fn = pfp_estimator_fn(model_version="v9.0.0", calc_mode="R2SCAN_PLUS_D3")

# 始構造を最適化
opt_feature = LBFGSASEOptFeature(fmax=0.01, n_run=1000, estimator_fn=estimator_fn)
opt_result = opt_feature(atoms)

# スキャン定義: Scan のリストのリストを渡す（順に実行される）
# ScanDistance / ScanDihedral の引数:
#   indices: 対象原子のインデックス（距離は2つ、二面角は4つ）
#   direction: 増加方向は 1、減少方向は -1
#   destination: この値（距離 [A] / 角度 [deg]）に達したらスキャン停止
#   exceed: destination を超えてよい量
#   fmax: 拘束スキャン力
scans = [
    [ScanDistance(indices=[0, 1], direction=1, destination=2.0, exceed=0.1, fmax=0.05)],
    # 二面角の例:
    # [ScanDihedral(indices=[2, 0, 1, 5], direction=-1, destination=90.0, exceed=1.0, fmax=0.05)],
]

scan_feature = RestScanFeature(
    scans=scans,
    trajectory="scan.traj",
    estimator_fn=estimator_fn,
)
scan_result = scan_feature(opt_result.atoms.ase_atoms, fmax=0.01, n_run=1000)
```

#### スキャン結果の利用

`RestScanFeatureResult` は `converged`（正常終了か）、`warnflag`（途中終了の理由: 1=最大ステップ到達, 2=ExitEnergy 到達）、`trajectory_logger`、`traj_path` を持ちます。

```python
# エネルギープロファイルの確認（kink_dx 以下の変位の image は同一とみなして省略）
fig = scan_result.plot(plt_name=None, kink_dx=0.1)

# 始状態と終状態の間の image 列を取得
images = scan_result.get_images(kink_dx=0.1)

# トラジェクトリとして保存
scan_result.to_traj("scan_path.traj")

# 初期推定経路として ReactionString（手法 C）へ渡す
scan_images = [atoms.copy() for atoms in scan_result.get_images()]
```

スキャン終点を終構造として使う場合は、必ず拘束を外して構造最適化（共通の前処理）してから端点に使用してください。

## 経路探索

準備した端点・初期推定経路に対し、手法 A / B / C のいずれかで経路探索を実行します（選択基準はワークフローの手法選択ガイドを参照）。

### 手法 A: ASE NEB

ASE 標準の NEB 法を直接利用する実装です。パラメータの細かな制御が可能です。

```python
from ase import Atoms
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
    model_version: str = "v8.0.0",
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
        estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
        calc = ASECalculator(estimator)
        image.calc = calc

    # NEB オブジェクトの構築
    neb = NEB(
        images,
        k=k,                              # バネ定数
        climb=True,                        # CI-NEB で正確な TS を取得
        method="improvedtangent",          # ASE 3.28 でデフォルトが変更されたため明示指定
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
    estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    calc = ASECalculator(estimator)
    image.calc = calc

# 誤ったパターン: Estimatorの共有（エラーになる）
# estimator = Estimator(...)
# for image in images:
#     image.calc = ASECalculator(estimator)  # NG
```

### 手法 B: NEBFeature

matlantis-features の `NEBFeature` は NEB 法の高レベルラッパーです。結果は `NEBFeatureResult` オブジェクトで返されます（解析は「手法 B の解析」参照）。

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
        model_version="v9.0.0",
        calc_mode="R2SCAN_PLUS_D3"
    ),
)
neb_result = neb(initial_atoms, final_atoms)
```

### 手法 C: ReactionStringFeature

matlantis-features の `ReactionStringFeature` です。反応経路の最適化と TS の最適化が自動で進行し、タイムアウト管理やバッチ処理に対応しています。TS の収束条件は `fmax_ts` / `fmax_rd` で指定します。手法とパラメータの詳細は Matlantis Guidebook（`matlantis/references/guidebook/matlantis-features/reaction_string.md`）、および API リファレンス（`matlantis/references/api/matlantis-features/generated/` の `ReactionStringFeature`）を参照してください。

入力は `List[List[MatlantisAtoms | Atoms]]` です。端点 2 点だけを渡すと内部で線形補間され、RestScan などで得た中間構造を含む軌跡を渡すとそれが初期推定経路になります。

```python
from matlantis_features.features.reaction import ReactionStringFeature
from matlantis_features.utils.calculators import pfp_estimator_fn
from datetime import timedelta
from ase import Atoms

def run_reaction_string(
    start_atoms: Atoms,
    end_atoms: Atoms,
    model_version: str = "v9.0.0",
    calc_mode: str = "R2SCAN_PLUS_D3",
    timeout_sec: int = 1800
):
    """
    ReactionStringFeatureを用いて反応経路を探索する。
    Returns:
        ReactionStringFeatureResult（解析は「手法 C の解析」参照）
    """
    # 入力トラジェクトリの作成とPBC補正
    traj_use = modify_trajectory_by_mic([start_atoms, end_atoms])

    # ReactionStringFeature の設定
    # fmax_rp: 反応経路の力の収束判定値
    # fmax_ts: 遷移状態における力の収束判定値
    # fmax_rd: 律速段階の力の収束判定値
    # fmax_eq: 平衡構造（端点等の極小点）の力の収束判定値
    rp = ReactionStringFeature(
        fmax_rp=0.05,
        fmax_ts=0.01,
        fmax_rd=0.005,
        fmax_eq=0.001,
        timeout=timedelta(seconds=timeout_sec),
        estimator_fn=pfp_estimator_fn(
            model_version=model_version,
            calc_mode=calc_mode
        )
    )

    # 入力はリストのリスト（バッチ処理対応）
    rp_result = rp([traj_use])
    return rp_result
```

#### RestScan 経路を初期値にする場合

準備 b で得たスキャン軌跡をそのまま初期推定経路として渡せます。スキャン経路の端点は最適化されていないため、`optimize_is` / `optimize_fs` で端点最適化を有効にします。

```python
scan_images = [atoms.copy() for atoms in scan_result.get_images()]

rp = ReactionStringFeature(
    optimize_is=True,   # 始状態を最適化してから経路探索
    optimize_fs=True,   # 終状態を最適化してから経路探索
    fmax_rp=0.05,
    fmax_ts=0.01,
    fmax_rd=0.005,
    fmax_eq=0.001,
    estimator_fn=estimator_fn,
)
rp_result = rp([scan_images])
```

#### バッチ処理

複数の反応経路を一括で探索する場合は、リストのリストとして入力します。

```python
# 複数反応のバッチ実行
rp_result = rp([
    [initial1, final1],
    [initial2, final2],
])
```

## 解析（使用した手法に対応）

解析コードは経路探索にどの手法を使ったかで異なります。A は生の Atoms リスト、B / C はそれぞれ専用の結果オブジェクトを返すためです。

### 手法 A の解析: NEBTools と反応プロファイル

NEB 計算直後の image は calculator と計算済みエネルギーを保持しているため、`get_potential_energy()` はキャッシュを返し再計算は発生しません。

#### NEBTools による収束評価

`NEBTools` を使用してバリアと最大力を同時に評価し、未収束の誤判定を防ぎます。

```python
from ase.mep import NEBTools
from typing import Tuple

def evaluate_neb_convergence(images) -> Tuple[Tuple[float, float], float]:
    """
    NEBToolsを用いてバリアとfmaxの二軸で収束を評価する。
    Arguments:
        images: NEB計算後のimageリスト
    Returns:
        barrier: (forward_barrier, backward_barrier)
        fmax: 最大力
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

`NEBTools.plot_band()` を使うと、スプライン補間付きのバンドプロファイルを 1 行で描画できます。

```python
fig = NEBTools(images).plot_band()
```

#### 反応プロファイルの plot

エネルギープロファイルから活性化エネルギーと反応エネルギーを算出します。手法 B / C の結果も `.ase_atoms` で Atoms リスト化すれば適用できます（calculator を持たない場合は 1 つの calculator を順次設定して評価してください）。

```python
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from typing import List

def analyze_reaction_profile(images: List[Atoms]) -> List[float]:
    """
    反応経路のエネルギープロファイルを計算・表示する。
    Forward/backward barrierの両方を出力する。
    """
    energies = [atoms.get_potential_energy() for atoms in images]

    # 反応座標: PBC を考慮した累積変位（MIC 補正）
    distances = [0.0]
    for i in range(1, len(images)):
        diff, _ = find_mic(
            images[i].get_positions() - images[i-1].get_positions(),
            images[i].get_cell(),
            pbc=True
        )
        distances.append(distances[-1] + np.linalg.norm(diff))

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

#### image ごとのエネルギー・最大力の plot

収束診断には、image（replica）番号に対するエネルギーと最大力の plot も有効です。特定の image だけ力が大きい場合、その近傍で経路が破綻しています。

```python
energies = [img.get_potential_energy() for img in images]
max_forces = [np.abs(img.get_forces()).max() for img in images]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(range(len(energies)), energies, 'o-')
axes[0].set_xlabel("replica")
axes[0].set_ylabel("energy (eV)")
axes[0].grid(True)
axes[1].plot(range(len(max_forces)), max_forces, 'o-')
axes[1].set_xlabel("replica")
axes[1].set_ylabel("max force (eV/A)")
axes[1].grid(True)
plt.show()
```

### 手法 B の解析: NEBFeatureResult

`NEBFeature` の結果は `NEBFeatureResult` オブジェクトです。主な属性は `converged`、`energy_log`（最適化ステップごとのエネルギー列のリスト）、`neb_images`（`.ase_atoms` を持つ MatlantisAtoms のリスト）です。

```python
import numpy as np

# 収束確認
print(f"converged: {neb_result.converged}")

# エネルギー: energy_log の最終要素が収束後の image ごとのエネルギー列
energies = neb_result.energy_log[-1]
barrier = np.max(energies) - np.min(energies)
print(f"barrier: {barrier:.3f} eV")

# 組み込みのバンドプロファイル plot（plotly）
fig = neb_result.plot()

# image の取り出し（MatlantisAtoms → ASE Atoms）
neb_images = [img.ase_atoms for img in neb_result.neb_images]

# 出力・保存
neb_result.to_traj("neb_path.traj")        # ASE トラジェクトリ
neb_result.to_gif("neb_path.gif")          # アニメーション GIF
neb_result.to_dataframe("neb_result.csv")  # pandas DataFrame / CSV
neb_result.save("neb_result.json")         # 結果オブジェクトの保存

# 保存した結果の読み込み
# from matlantis_features.features.common.neb import NEBFeatureResult
# neb_result = NEBFeatureResult.load("neb_result.json")
```

### 手法 C の解析: ReactionStringFeatureResult

`ReactionStringFeature` の結果は `ReactionStringFeatureResult` オブジェクトです。経路は途中のすべての安定点（EQ）で分割され、セグメントのリストとして格納されます。

主な属性:

| 属性 | 内容 |
| :--- | :--- |
| `num_segments` | 経路中の山（障壁）の数 |
| `reaction_string_images` | EQ で区切られたセグメント構造 `List[List[MatlantisAtoms]]`。各セグメントは `[EQ, ..., TS, ..., EQ]` で、セグメント末尾と次セグメント先頭は同一構造 |
| `transition_state_indices` | 各セグメント内の TS の index。`None` は TS 未計算（`fmax_ts=inf` 指定、タイムアウト、障壁が許容値未満など） |
| `imaginary_eigenvectors` | 計算された TS の虚振動モードの固有ベクトル |

```python
# 収束・セグメント情報
print(f"converged: {rp_result.converged}")
print(f"num_segments: {rp_result.num_segments}")
print(f"TS indices: {rp_result.transition_state_indices}")

# 組み込みのエネルギープロファイル plot（plotly; EQ と TS がマーカー表示される）
rp_result.plot_energy()

# セグメント構造を 1 本の経路に展開する
# （セグメント境界で構造が重複するため、先頭セグメント以外は 2 番目以降を連結する）
rs = rp_result.reaction_string_images
atomslist = [rs[0][0].ase_atoms]
atomslist += [matatoms.ase_atoms for reac in rs for matatoms in reac[1:]]

# TS 構造の取得（例: 最初のセグメント）
ts_idx = rp_result.transition_state_indices[0]
if ts_idx is not None:
    ts_atoms = rs[0][ts_idx].ase_atoms

# 出力・保存
rp_result.to_traj("reaction_string.traj")
rp_result.to_gif("reaction_string.gif")
rp_result.save("rp_result.json")

# 保存した結果の読み込み
# from matlantis_features.features.reaction.reaction_string import ReactionStringFeatureResult
# rp_result = ReactionStringFeatureResult.load("rp_result.json")
```

## 遷移状態の精密化と検証

### Sella による TS 精密化

手法 A / B（NEB / CI-NEB）で得られる最高エネルギー image は遷移状態の近似です。厳密な活性化エネルギーが必要な場合は、[Sella](https://github.com/zadorlab/sella) を用いて TS 候補を一次の鞍点（虚振動が 1 つだけの点）に収束させます。

手法 C（ReactionStringFeature）の TS は `fmax_ts` を収束条件として最適化済みです（内部実装の記述は手法 C のセクション参照）。さらに高精度化したい場合や別アルゴリズムで検証したい場合は、`transition_state_indices` で TS を取り出して Sella を適用できます。

```python
# 事前インストール: pip install sella
import numpy as np
from sella import Sella
from ase.constraints import FixAtoms
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

# TS 候補の取得
# 手法 A / B: 最高エネルギー image を抽出
energies = [img.get_potential_energy() for img in images]
ts = images[int(np.argmax(energies))].copy()
# 手法 C: transition_state_indices から取得（「手法 C の解析」参照）
# ts = rs[0][rp_result.transition_state_indices[0]].ase_atoms.copy()

# 表面系では下層固定などの拘束を経路探索時と同じ条件で設定する
# ts.set_constraint(FixAtoms(indices=[...]))

ts.calc = ASECalculator(
    Estimator(model_version="v8.0.0", calc_mode="PBE")
)

# Sella はデフォルトで一次の鞍点に向かって最適化する
# 収束条件は ReactionString の fmax_ts と同水準の 0.01 eV/A を推奨
# （精密化が目的のため、NEB バンドの収束判定 0.05 より厳しくする）
ts_opt = Sella(ts)
ts_opt.run(fmax=0.01)

print(f"TS energy:  {ts.get_potential_energy():.3f} eV")
print(f"max force:  {np.abs(ts.get_forces()).max():.4f} eV/A")
```

### 虚振動解析による検証

厳密な遷移状態確認には、得られた TS 構造に対して振動解析を行い、虚の振動数が 1 つだけ存在することを確認します。計算コストが非常に高いため、まずはエネルギープロファイルの形状で判断し、必要に応じて実施してください。

```python
from matlantis_features.features.common.vibration import VibrationFeature

vib = VibrationFeature(delta=0.01, estimator_fn=estimator_fn)
vib_result = vib(ts_atoms)
# 虚振動数が1つだけであることを確認
```

## ベストプラクティス

### 端点の最適化は必須

始状態と終状態は必ず事前に構造最適化を完了させてください。最低でも `fmax < 0.01` eV/A、理想的には `fmax < 0.001` eV/A まで収束させます（公式チュートリアルでは 0.01〜0.005、ReactionStringFeature の公式例では `fmax_eq=0.001` が使われています）。端点の最適化が緩いと経路が不安定になり、非物理的な結果やバリアの過大・過小評価の原因になります。なお、NEB バンド自体の収束判定（0.05 程度）と端点の収束判定は別物であることに注意してください。

### Atom Mapping の検証

始状態と終状態で原子のインデックス順序が完全に一致している必要があります。Atoms[0] が始状態で「炭素 0」なら、終状態でも「炭素 0」でなければなりません。ここがずれていると、原子が物理的にありえない移動を行い、エネルギーが発散します。

### PBC 補正の徹底

周期境界条件をまたぐ原子移動がある場合、`modify_trajectory_by_mic` による補正を必ず行ってください。「原子が突然消えた」「結合が異常に伸びた」ように見える場合は、PBC の補正ミスが原因です。

### 補間法の使い分け

| 補間法 | 特徴 | 推奨用途 |
| :--- | :--- | :--- |
| **Linear** | 単純な直線補間。高速だが初期経路が粗くなりやすい | 端点が近く単純移動の場合 |
| **IDPP** | Image Dependent Pair Potential。原子間距離を考慮した補間 | 初期経路の粗さが強い場合 |

### 最適化器の選択

NEB 本計算では `FIRE` を使うと、バンド全体の力が大きい初期段階で LBFGS より安定するケースがあります。

### NEB の 2 段階実行

NEB バンドの収束判定は `fmax=0.05` eV/A 以下が推奨です（小さすぎると収束に時間がかかります）。1 回目は緩い条件（`fmax=0.2` 程度）で実行して経路の妥当性を確認し、妥当であれば 2 回目で厳しい条件に締めると効率的です。

### RestScan の経路をそのまま活性化障壁として解析に使わない

RestScan の経路は最適な反応経路上を通らず、TS も通らないことが多いため、スキャンプロファイルの最高点は活性化障壁の概算にとどめてください。活性化障壁が必要な場合は ReactionString などで経路と TS を最適化します。

### 終了判定の二軸化

`NEBTools.get_barrier()` だけでなく `get_fmax()` も同時評価し、「見かけ上バリアが出たが未収束」の誤判定を防いでください。

### 虚振動解析による TS 確認

厳密な遷移状態確認には虚振動解析（虚の振動数が 1 つだけ存在すること）を実施してください（「遷移状態の精密化と検証」参照）。

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| エネルギーが発散する | Atom Mapping のずれ（始状態と終状態で原子インデックスの対応が不一致） | 始状態と終状態の原子順序を確認し、必ず同じ原子が同じインデックスに対応するように修正する |
| 原子が突然消える、結合が異常に伸びる | PBC の補正ミス | `modify_trajectory_by_mic` を必ず適用する。`wrap()` を補間前に使わない |
| `MultiCalculatorUseDetected` | 1 つの Estimator を複数の Calculator で共有した | 各 image に個別の Estimator + Calculator を生成する |
| NEB が収束しない | 初期経路が不適切、または steps 不足 | IDPP 補間を試す、image 数を増やす、steps を増やす。RestScan で初期経路を作る |
| Timeout で計算が中断される | ReactionStringFeature のデフォルトタイムアウト超過 | `timeout` 引数で時間を延長する、中間構造を増やして初期パスを改善する |
| 「見かけ上バリアが出たが未収束」 | `get_barrier()` のみで判定し `get_fmax()` を確認していない | `NEBTools.get_fmax()` も同時に評価し、fmax < 0.05 を確認する |
| 非物理的な経路が得られる | 端点の最適化が不十分 | 端点を最低でも `fmax < 0.01`（理想的には 0.001）まで最適化してから NEB を実行する |
| `transition_state_indices` に `None` が含まれる | TS 未計算（`fmax_ts=inf` 指定、タイムアウト、障壁が許容値未満） | 仕様どおりの動作。TS が必要なセグメントは `fmax_ts` / `fmax_rd` を有限値にして再実行する |
| Sella が TS に収束しない、別の鞍点に行く | TS 候補が鞍点から遠い | NEB をより収束させてから最高エネルギー image を抽出する。image 数を増やして TS 近傍の解像度を上げる |
| `ConcurrentUseDetected` | Calculator を複数スレッドで共有した | 各スレッドで新しく Calculator を生成する |

## 関連ガイド

- **構造最適化** (optimization/SKILL.md): 端点の最適化手順の詳細
- **吸着** (adsorption/SKILL.md): 吸着構造の構築・吸着エネルギー計算・吸着サイト探索（吸着構造を端点とした NEB の前準備）
- **物性計算** (property-analysis/SKILL.md): 振動解析による TS 検証
- **可視化** (visualization/SKILL.md): 反応経路の構造確認
- **外部データ連携** (external-data/SKILL.md): 反応物・生成物の構造取得
