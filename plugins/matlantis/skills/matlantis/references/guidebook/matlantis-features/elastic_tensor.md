# ElasticTensorFeature ​

## はじめに ​

弾性係数テンソル(Elastic tensor)は固体において重要な物理量です。 [ElasticTensorFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature) は、微小ひずみに対する応力を計算し、弾性係数テンソルを求めます。 体積弾性率、ヤング率、ポアソン比などさまざまな力学特性が弾性係数テンソルから計算することができます。

## 計算手法 ​

弾性係数テンソルは以下のように計算します。

  * まず、入力構造の原子位置とセル形状の両方を最適化します。
  * 格子ベクトルにひずみを加えます。
  * 格子形状を固定しながら原子位置のみを緩和し、変形に伴う応力を計算します。
  * ひずみと応力の関係を線形フィッティングすることで、6×6の弾性係数テンソルを求めます。



## 弾性係数テンソルの対称性 ​

弾性係数テンソルは6x6の対称行列なので、21個の独立な行列要素を持ちます。 しかし、結晶格子の対称性を考慮するとより多くの拘束条件がかけられ、独立な要素数を絞り込むことができます。 以下に、7つの [結晶系](https://ja.wikipedia.org/wiki/%E7%B5%90%E6%99%B6%E7%B3%BB) での既約テンソルを示します。[1,2]

### 立方晶系 (cubic lattice) ​

3個の独立な行列要素を持ちます。

### 六方晶系 (hexagonal lattice) ​

5個の独立な行列要素を持ちます。

### 正方晶系 (tetragonal lattice) ​

7個の独立要素を持ちます。

### 三方晶系 (rhombohedral (or trigonal) lattice) ​

7個の独立要素を持ちます。

### 直方晶系 (orthorhombic lattice) ​

9個の独立要素を持ちます。

### 単斜晶系 (monoclinic lattice) ​

13個の独立な行列要素を持ちます。

### 三斜晶系 (triclinic lattice) ​

21個の独立な行列要素を持ちます。

## Tips ​

  * `ElasticTensorFeature` による弾性係数テンソルの計算を行う前に、原子位置とセル形状の両方を構造最適化しておく必要があります。
  * `ElasticTensorFeature` は `check_symmetry=True` にした場合、格子の対象性を考慮して既約弾性係数テンソルを出力します。



## Reference ​

[1] <https://arxiv.org/abs/1812.03367>

[2] <http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf>
