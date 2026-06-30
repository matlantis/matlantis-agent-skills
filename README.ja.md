# Matlantis Agent Skills

Claude Code / Codex 向けの Matlantis / PFP 原子シミュレーションスキルを提供するプラグインマーケットプレイスです。

## Overview

Matlantis Agent Skills は、[Matlantis](https://matlantis.com/ja/) プラットフォーム上での原子シミュレーションを Claude Code と Codex がサポートするためのスキル集です。PFP (Preferred Potential) ニューラルネットワークポテンシャルを活用した構造最適化、分子動力学、物性計算などの Python コード生成を、16 のスキルモジュールがガイドします。

各スキルは、ASE (Atomic Simulation Environment) を基盤とした実装パターン、ベストプラクティス、よくあるエラーと対処法を体系的にまとめています。

## Skills Catalog

| カテゴリ | スキル | 概要 |
|---------|--------|------|
| **基盤** | [matlantis](#matlantis) | プラットフォーム全体概要、PFP計算モード選択、ライブラリ構成 |
| **基盤** | [setup](#setup) | Calculator初期化、calc_mode設定、構造ファイルI/O |
| **基盤** | [ase-basics](#ase-basics) | ASE基礎操作、Atomsオブジェクト、分子・結晶の構築 |
| **構造構築** | [modeling](#modeling) | 表面スラブ生成、SMILES→3D変換、置換・ドーピング、液体・アモルファス |
| **構造構築** | [external-data](#external-data) | PubChem・Materials Project等からの構造取得、ファイル形式変換 |
| **計算** | [optimization](#optimization) | 固定セル・可変セル構造緩和、対称性保存、多段階最適化 |
| **計算** | [dynamics](#dynamics) | NVE/NVT/NPT分子動力学、メタダイナミクス、輸送特性 |
| **計算** | [reaction](#reaction) | NEB反応経路探索、遷移状態探索、活性化エネルギー |
| **計算** | [adsorption](#adsorption) | 吸着構造構築、吸着エネルギー計算、吸着サイト探索 |
| **計算** | [property-analysis](#property-analysis) | 弾性定数、フォノン、QHA、拡散係数、粘度、熱伝導率 |
| **ワークフロー** | [integrated-workflows](#integrated-workflows) | 構造スクリーニング、吸着キャンペーン、欠陥解析パイプライン |
| **ワークフロー** | [background-job](#background-job) | 長時間計算のバックグラウンド実行・チェックポイント管理 |
| **高速化** | [light-pfp](#light-pfp) | Light-PFPカスタムモデル(15-40倍高速)の学習・推論 |
| **ユーティリティ** | [pfcc-extras](#pfcc-extras) | 可視化・構造処理・ジョブ実行制御のユーティリティ群 |
| **ユーティリティ** | [visualization](#visualization) | 原子構造・トラジェクトリの可視化、nglview、POV-Ray |
| **インフラ** | [storage](#storage) | Group Drive永続ストレージ、チーム間データ共有 |
| **インフラ** | [ssh](#ssh) | SSH接続設定、VSCode Remote-SSH、WinSCP連携 |

## Getting Started

### 前提条件

- Claude Code で使う場合は **Claude Code** CLI がインストールされていること
  - プラグイン機能に対応したバージョンが必要です
  - インストール: https://claude.ai/code
- Codex plugin として使う場合は **Codex** CLI がインストールされていること
  - インストール: https://developers.openai.com/codex/cli
- **GitHub Copilot** — 公式のプラグインサポートはありません。手動での設置方法は下記の [GitHub Copilot (参考)](#github-copilot-参考) セクションを参照してください

### インストール手順

#### Claude Code

Claude Code の REPL 上で以下のコマンドを実行します。

```
# 1. マーケットプレイスを登録する
/plugin marketplace add https://github.com/matlantis/matlantis-agent-skills
```

```
# 2. プラグインマネージャーを開く
/plugin

# 3. 「Marketplaces」タブから matlantis プラグインを有効化する
```

#### Codex

シェル上で以下のコマンドを実行します。

```bash
# 1. マーケットプレイスを登録する
codex plugin marketplace add matlantis/matlantis-agent-skills

# 2. マーケットプレイスから plugin をインストールする
codex plugin add matlantis@matlantis-agent-skills
```

インストールまたは有効化後、Matlantis 関連のコードや質問に対して、エージェントが自動的に適切なスキルを参照するようになります。

#### GitHub Copilot (参考)

> **注意:** GitHub Copilot はこのプラグインマーケットプレイスの公式サポート対象外です。以下はあくまで参考のための手動設置方法であり、Copilot のバージョンやエディタ設定によって動作が異なる場合があります。

このリポジトリをクローンし、`skills` ディレクトリをプロジェクトの `.copilot/` 以下にコピーします。

```bash
git clone https://github.com/matlantis/matlantis-agent-skills.git
cp -r matlantis-agent-skills/plugins/matlantis/skills <your-project>/.copilot/
```

## スキル一覧と使用例

### matlantis

プラットフォーム全体の概要スキル。PFP の計算モード（PBE, R2SCAN, WB97XD 等）の選択指針、ライブラリ階層、重要な制約事項をカバーします。

```
「銅の表面エネルギーを計算したい。どの calc_mode を使うべき？」
```

### setup

Calculator の初期化と計算環境セットアップ。Estimator/ASECalculator の構成パターン、構造ファイルの入出力、バッチ実行設定を扱います。

```
「PFP の Calculator を R2SCAN モードで初期化するコードを書いて」
```

### ase-basics

ASE の基礎操作。Atoms オブジェクトの作成・編集、分子・結晶の構築、周期境界条件、ファイル I/O を扱います。

```
「CIF ファイルから結晶構造を読み込んで 2x2x2 のスーパーセルを作りたい」
```

### modeling

構造モデリング。表面スラブ生成（ASE builders / pymatgen）、SMILES からの 3D 分子生成（RDKit）、原子置換・ドーピング、液体・アモルファス・高エントロピー合金の構築を扱います。

```
「FCC銅の(111)表面スラブを4層で作成して、真空層を15Å入れたい」
```

### external-data

外部データベースからの構造取得。PubChem（有機分子）、Materials Project（無機結晶）、AFLOW、COD との連携と、ファイル形式変換を扱います。

```
「Materials Project から BaTiO3 の結晶構造を取得して最適化したい」
```

### optimization

構造最適化（構造緩和）。固定セル最適化（原子位置のみ）と可変セル最適化（格子ベクトル含む）、LBFGS/BFGS/FIRE オプティマイザ、対称性保存最適化を扱います。

```
「可変セル最適化で結晶構造を緩和したい。対称性は保存して」
```

### dynamics

分子動力学（MD）シミュレーション。NVE/NVT/NPT アンサンブル、Langevin サーモスタット、PLUMED メタダイナミクス、トラジェクトリ解析、輸送特性（拡散・粘度・熱伝導率）の算出を扱います。

```
「300K の NVT アンサンブルで 100ps の MD シミュレーションを実行したい」
```

### reaction

反応経路探索・遷移状態探索。NEB（Nudged Elastic Band）、CI-NEB、ReactionStringFeature を用いた反応障壁計算を扱います。

```
「表面上での CO 解離反応の活性化エネルギーを NEB で求めたい」
```

### adsorption

表面への分子吸着。add_adsorbate による吸着構造の構築、D3 分散補正を用いた吸着エネルギー計算、吸着サイト探索（Optuna によるブラックボックス最適化を推奨）を扱います。

```
「Pt(111) 上の CO の最安定吸着サイトを探索して吸着エネルギーを計算したい」
```

### property-analysis

物性計算と MD 後処理。弾性定数・弾性率、フォノンバンド構造・状態密度、準調和近似（QHA）、拡散係数、粘度、熱伝導率、動径分布関数（RDF）を扱います。

```
「シリコンのフォノンバンド構造と状態密度を計算したい」
```

### integrated-workflows

統合ワークフロー（パイプライン）テンプレート。構造スクリーニング（多構造バッチ最適化＋エネルギーランキング）、吸着キャンペーン（サイト探索＋ランキング）、欠陥解析パイプラインを提供します。

```
「10種類の合金組成を一括スクリーニングして安定性をランキングしたい」
```

### background-job

長時間計算のバックグラウンド実行。ノートブックセッションから切り離した計算実行、チェックポイント/リスタート、ジョブ監視を扱います。

```
「大規模 MD をバックグラウンドで実行して、途中経過をモニタリングしたい」
```

### light-pfp

Light-PFP カスタムモデル。特定材料向けの MTP（Moment Tensor Potential）を学習し、PFP の 15-40 倍高速な推論を実現するワークフローを扱います。

```
「銅系の材料に特化した Light-PFP モデルを作成したい」
```

### pfcc-extras

pfcc_extras ユーティリティツール群。可視化レイヤー（show_gui, view_ngl）、構造処理レイヤー（部分占有率、分子ラッピング）、MD イベントスケジューリング、ジョブ実行制御を扱います。

```
「複数の構造最適化ジョブをリソース制限付きで並列実行したい」
```

### visualization

原子構造・トラジェクトリの可視化。pfcc_extras.show_gui による簡易表示、nglview によるインタラクティブ可視化、トラジェクトリアニメーション、POV-Ray 高品質レンダリングを扱います。

```
「MD トラジェクトリをアニメーション表示して、特定フレームを高品質で描画したい」
```

### storage

Group Drive 永続ストレージ。クラウド上のオブジェクトストレージへのファイルアップロード・ダウンロード、チーム間データ共有を扱います。

```
「計算結果を Group Drive に保存してチームメンバーと共有したい」
```

### ssh

SSH 接続設定。WebSocket ベースの SSH トンネリング、VSCode Remote-SSH 連携、WinSCP によるファイル転送を扱います。

```
「VSCode から Matlantis 環境に Remote-SSH で接続したい」
```

## ライセンス

Copyright 2026 Matlantis Corporation.

本プロジェクトは [Apache License 2.0](LICENSE) のもとで公開されています。

## サードパーティライセンス

このリポジトリには、Creative Commons Attribution 4.0 International (CC BY 4.0) ライセンスの下で提供されているMatlantisのドキュメントが含まれています。
Copyright 2021-present, Preferred Networks, Inc. and ENEOS Corporation.
