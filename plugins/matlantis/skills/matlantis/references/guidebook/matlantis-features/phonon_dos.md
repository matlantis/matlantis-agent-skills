# PostPhononDOSFeature ​

## はじめに ​

[PostPhononDOSFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.phonon.dos.PostPhononDOSFeature.html#matlantis_features.features.phonon.dos.PostPhononDOSFeature) はフォノンの状態密度 (density of states, DOS) を求めるための機能です。

## 計算手法 ​

結晶のフォノンは多様な振動モードを持ちます。DOS とは、微小な振動数区間 (ω から ω+dω) の中にあるフォノンモードの数を表す量です。 波数 q→ に対するフォノン振動数の計算方法については `PostPhononBandFeature` のページをご参照ください。

### k点メッシュ ​

フォノンDOSを計算するには、第一ブリルアンゾーン内の q→ に対する一様メッシュをとり、各メッシュ点(k-point)でのフォノン振動数を求める必要があります。 以下では、 `PostPhononDOSFeature` でサポートされているk点メッシュの2種類のとり方、Monkhorst-Pack mesh および Gamma-centered mesh を紹介します。

#### (1) Monkhorst-Pack ​

Monkhorst-Pack のk点メッシュを指定した場合、各軸のk点数を指定すると自動的にメッシュを生成します。 Monkhorst-Pack は以下の式に従ってk点を取ります。

k→=n1+12N1b→1+n2+12N2b→2+n3+12N3b→3ni=0,1,2...,Ni

ここで、 b→i は逆格子ベクトル、 Ni は b→i 方向のk点数を表します。

#### (2) Gamma-centered ​

Gamma-centered メッシュは Monkhorst-Pack mesh を (12,12,12) だけシフトし、ガンマ点 q→=(0,0,0) が含まれるようにしたものです。

### フォノン状態密度 ​

生成したk点メッシュ上でフォノン振動数を計算します。 得られたフォノン振動数の分布から、DOSは以下のように計算されます:

D(ω)=1N∑m,q→δ(ω−ωm,q→)

ここで、 N は全フォノンモードの総数、 ωm,q→ は 波数 q→ における第 m 番目の振動数を表します。

`PostPhononDOSFeature` では、DOSを計算する際の振動数区間をパラメーター `freq_min`, `freq_max`, `freq_bin` (区間の最小値、最大値、微小区間の幅) で指定します。 各微小区間に属するフォノン振動数の数をk点メッシュ上でカウントし、DOSが計算されます。DOSをプロットすると図のようになります。

上図はフォノンDOSのプロット例です。k点メッシュをより細かくとれば、さらに滑らかなプロットにすることができます。

#### partial DOS ​

partial (projected) DOS とは、DOSに対するある元素、または原子からの寄与を表す量です。 ある原子 j からの寄与は次のように計算することができます:

Dj(ω)=1N∑m,q→|em,q→j|2δ(ω−ωm,q→)

|em,q→| フォノンモード (m,q→) の固有ベクトル(単位長さに正規化済み)を表します。 |em,q→j| は |em,q→| の中で原子 j に対応する要素を取り出したベクトルです。 ある元素 e に対するpartial DOS は、その元素である原子 j∈e の Dj(ω) を合計することで得ることができます。

`PostPhononDOSFeature` では、 `partial=True` を指定することで各元素に対するpartial DOSを計算します。

#### フォノン分散との関係 ​

フォノン分散 (`PostPhononBandFeature` 参照) はフォノンDOSと密接に関係しています。 一般に、フォノンバンドが比較的平らになっている振動数区間においてフォノンDOSが極大となります。

## Examples ​

ランチャーページから開ける `Matlantis-features: Phonon` および `Matlantis-features: Phonon advanced usage` をご参照ください。
