# DFT-D3補正に関するアップデート ​

## 概要 ​

  * v7.0.0以降のPFPで提供されるD3補正のパラメータを、以下のように変更しました 
    * cutoff_smoothingをnoneからpolyに変更
    * カットオフ距離を21 Åから14 Åに変更
  * パラメータ変更前後での有機物の体積と吸着エネルギーの変化を調べ、精度への悪影響がないことを確認しました
  * 上記の変更によってPFP＋D3補正を用いた場合に、より原子数の多い系を扱えるようになることを確認しました
  * v6.0.0以前に関しては、後方互換性のためにD3補正のパラメータは変更されません。



## 背景 ​

以前のPFPではD3補正のカットオフ距離をtorch-dftdのオリジナルの値である95 Bohrよりも短い、40 Bohr(21 Å)に設定していました。その後PFPのパフォーマンスの向上により、PFPではなくD3補正が扱える原子数のボトルネックとなる状況が出てきました。そのため推論精度に悪影響を与えない範囲で、より多くの原子数の系を扱えるようにすることを目的に、適切なカットオフ距離を設定しました。

また以前のバージョンでは内部的な不具合のため、カットオフ距離付近での平滑化に関するパラメータであるcutoff_smoothingがnoneに設定されており、適切な平滑化がされない状況でした。今回この不具合も同時に修正し、cutoff_smoothingをpolyに設定しました。なお、v6.0.0以前のバージョンでは、後方互換性を維持するために、cutoff_smoothingがnoneに設定されています。

## 体積の再現性 ​

ここではD3補正のcutoff_smoothingおよびカットオフ距離の変更前後で、PFP Validationで議論している体積の再現性が低下しないことを確認します。 検証に用いる構造としては従来のPFP Validation 2024.03と同様、Crystallography Open Database (COD) [1]の構造を使用しました。 DFTによる構造最適化後の構造を初期構造とし、PFP v5.0.0で格子定数を含む構造最適化を行い、D3補正のパラメータ変更前を基準にした場合の体積の再現性を各カットオフ距離で確認しました。 なおD3補正はPFPで求めたエネルギーに加算するモデルであり、PFP計算とは独立に計算されるため、以降での再現結果についてPFPのバージョンによる違いは無視できると考えています。

`ase.Atoms`である`atoms`に対して、以下のコードを実行しています。
    
    
    from ase.optimize import LBFGS
    from ase.constraints import UnitCellFilter
    
    ecf = UnitCellFilter(atoms)
    opt = LBFGS(ecf)
    opt.run(fmax=0.03, steps=2000)

PFPでは`PBE_PLUS_D3`を指定しており、DFT計算では`PBE`での条件に加えて`IVDW=12`を指定することでDFT-D3補正を有効にしています。

図中のVref は変更前のD3パラメータである、カットオフ距離21 Åおよびcutoff_smoothing=noneを用いた場合の測定結果です。 D3カットオフ距離が14 Åの場合は、cumulative ratioが0.9となる90パーセンタイルの位置はパラメータ変更前と比べて変化しませんでした。 しかしそれより短いカットオフ距離では、90パーセンタイルの位置がrelative errorの大きい側にずれることが確認できました。 そのため有機物の体積に関する調査からは、D3カットオフ距離を14 Åまで減らしても従来と比べて精度に悪影響がないことを確認できました。

## 吸着エネルギーの再現性 ​

ここではD3カットオフ距離を14 Åに変更することによって、吸着エネルギーの精度に悪影響がないことを確認します。 検証には吸着エネルギーに関するベンチマーク論文[2]に記載されている吸着構造のうち、解離反応が起きない以下の系を用いました。

金属| 吸着分子  
---|---  
Ni| CO  
Pt| CO  
Pd| CO  
Pd| CO  
Rh| CO  
Ir| CO  
Cu| CO  
Ru| CO  
Co| CO  
Pt| NO  
Pd| NO  
Pd| NO  
Ni| O  
Ni| O  
Pt| O  
Rh| O  
Pt| H  
Ni| H  
Ni| H  
Rh| H  
Pd| H  
Pt| I  
Pt| CH3OH  
Pt| CH3  
Pt| CH4  
Pt| C2H6  
Pt| C3H8  
Pt| C4H10  
Pt| C6H6  
Cu| C6H6  
  
これらの吸着構造についてPFP v5.0.0およびD3補正を用いて吸着エネルギーを求め、D3補正のcutoff_smoothingおよびカットオフ距離の変更前後での値を比較しました。 ここで吸着分子間の相互作用が結果に影響しないように、吸着分子間の距離が元のD3カットオフ距離の21 Åよりも広くなるように設定しました。

上図から、カットオフ距離の変更前後を比較した場合のMAEは0.0024 eVであることがわかります。 これは吸着エネルギーで議論される典型的なエネルギースケールである、0.01 eVに比べて無視できる大きさと言えます。

そのためcutoff_smoothingとカットオフ距離の変更によって、吸着エネルギーの予測精度は悪影響を受けないことを確認できました。

## パラメータ変更による推論速度および扱える系のサイズの向上 ​

ここではパラメータ変更により、PFP+D3の推論速度および扱える原子数が増加することを確認します。

Ni3Alについて、原子数を増やしながらPFP+D3での力の推論速度を測定しました。

図中のreferenceは変更前のD3パラメータである、カットオフ距離21 Åおよびcutoff_smoothing=noneを用いた場合の測定結果です。 上図から、カットオフ距離の変更前では約5000原子以上の系に対する推論に失敗するのに対し、14 Åに変更した場合はその約３倍の原子数の系を扱えることが分かりました。 加えて、速度面でも元のカットオフ距離の場合と比較して高速になっています。 またD3以外の計算との兼ね合いで、扱える原子数はカットオフ距離14 Åで頭打ちとなり、それよりカットオフ距離を短くしてもわずかな速度面の向上を除いてメリットが薄いと言えます。

精度および扱える原子数の大きさに関する以上の検証をもとに、PFP v7.0.0以降ではD3カットオフ距離を14 ÅとしてD3補正を提供します。

## 参考文献 ​

* * *

  1. “Crystallography Open Database” <https://www.crystallography.net/cod/> ↩︎

  2. Wellendorff, Jess, Trent L. Silbaugh, Delfina Garcia-Pintos, Jens K. Nørskov, Thomas Bligaard, Felix Studt, and Charles T. Campbell. “A Benchmark Database for Adsorption Bond Energies to Transition Metal Surfaces and Comparison to Selected DFT Functionals.” Surface Science, Reactivity Concepts at Surfaces: Coupling Theory with Experiment, 640 (October 1, 2015): 36–44. <https://doi.org/10.1016/j.susc.2015.03.023>. ↩︎


