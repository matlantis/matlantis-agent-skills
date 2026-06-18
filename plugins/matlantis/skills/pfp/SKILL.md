---
name: pfp
description: >
  PFP (Preferred Potential) モデルの計算モード (calc_mode) 選択、対応元素、適用範囲を扱うスキルです。
  PFP, Preferred Potential, calc_mode選択, EstimatorCalcMode,
  PBE, PBE_U, PBE_PLUS_D3, PBE_U_PLUS_D3, R2SCAN, R2SCAN_PLUS_D3, WB97XD,
  旧モード名 CRYSTAL_U0, CRYSTAL, CRYSTAL_U0_PLUS_D3, CRYSTAL_PLUS_D3, MOLECULE,
  Hubbard U, D3分散補正,
  対応元素, supported_elements,
  PFP適用範囲外, 励起状態, 帯電系, 電子物性, HOMO, LUMO, バンドギャップ, 粗視化, 磁気秩序,
  model_version (calc_mode との組み合わせの範囲で)
  に関するコード生成・モード選択判断時に使用してください。
---

# PFP モデル・計算モード (calc_mode)

## 概要

PFP (Preferred Potential) は Matlantis のコアとなる、グラフニューラルネットワークベースの汎用原子間ポテンシャルです。
PFP は複数の汎関数に基づく計算モード (calc_mode) を提供しており、対象系の性質や用途に応じて適切なモードを選択する必要があります。
計算モードは結果の精度と物理的妥当性を大きく左右する、最も重要な設定の 1 つです。

## calc_mode の選択フローチャート

```
計算対象は？
+-- 遷移金属酸化物・強相関系 -> PBE_U (分散力が重要なら PBE_U_PLUS_D3)
+-- 分子結晶・液体・吸着系
|   +-- 実験値との定量比較が重要 -> R2SCAN_PLUS_D3
|   +-- 定性的・スクリーニング用途 -> PBE_PLUS_D3
+-- バルク結晶の形成エネルギー
|   +-- 実験値と比較 -> R2SCAN
|   +-- Materials Project と比較 -> PBE (または PBE_U)
+-- 有機分子系
|   +-- 単分子系、9元素のみ (H,C,N,O,F,P,S,Cl,Br) かつDFT結果の再現用途-> WB97XD
|   +-- それ以外 -> R2SCAN_PLUS_D3
+-- 迷ったとき -> PBE　(分散力が重要なら PBE_PLUS_D3、遷移金属酸化物を含むなら PBE_U_PLUS_D3)
```

## モード一覧

| モード | 汎関数 | 対応元素 | 主な用途 |
|--------|--------|----------|----------|
| `PBE` (旧: `CRYSTAL_U0`)| GGA-PBE | 96元素（v3.0.0 ~ v6.0.0では72元素, v2.0.0では55元素）| バルク系で汎用（v8.0.0以前でデフォルト） |
| `PBE_U`  (旧: `CRYSTAL`)| GGA-PBE + Hubbard U | 96元素（v3.0.0 ~ v6.0.0では72元素, v2.0.0では55元素） | 遷移金属酸化物・強相関系 |
| `PBE_PLUS_D3` (旧: `CRYSTAL_U0_PLUS_D3`) | GGA-PBE + D3分散力補正 | 96元素（v3.0.0 ~ v6.0.0では72元素, v2.0.0では55元素） | 分子結晶、液体（定性）、複合系 |
| `PBE_U_PLUS_D3`  (旧: `CRYSTAL_PLUS_D3`)| GGA-PBE + U + D3 | 96元素（v3.0.0 ~ v6.0.0では72元素, v2.0.0では55元素） | 遷移金属酸化物+有機配位子 |
| `R2SCAN` | meta-GGA R2SCAN | 96元素（v8.0.0では72元素, v7.0.0以前では使用不可） | バルク系で汎用・高精度（v9.0.0でデフォルト） |
| `R2SCAN_PLUS_D3` | meta-GGA R2SCAN + D3 | 96元素（v8.0.0では72元素, v7.0.0以前では使用不可） | 汎用・高精度（v8.0.0では定量評価） |
| `WB97XD`  (旧: `MOLECULE`)| hybrid omegaB97X-D | 9元素のみ | 有機分子単分子（非推奨、DFT結果の再現用途のみ） |

## 対応元素

- PBE/PBE_U/R2SCAN: 最大96元素 (H〜Cm)、バージョン依存
- PBE_UのHubbard U補正の対象: 特定の遷移金属(V, Cr, Mn, Fe, Co, Ni, Cu, Mo, Wの9元素)
- WB97XD: 9元素のみ (H, C, N, O, F, P, S, Cl, Br)
- 実行時確認: `estimator.supported_elements()`

## PFP の適用範囲外

PFP は以下を扱えません:
- 励起状態
- 帯電系
- 電子物性 (HOMO/LUMO, バンドギャップ)
- 粗視化モデル
- 磁気秩序 (Hubbard U 以上の精度)
- Uパラメーター数値の変更（特定のUパラメーターで学習データを収集しているため）
