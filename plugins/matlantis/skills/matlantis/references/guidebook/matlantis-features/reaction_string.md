# ReactionStringFeature ​

## はじめに ​

[ReactionStringFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature) は反応経路と遷移状態 (TS) を求めるのに使います。実際の使用例についてはMatlantis Examplesの `Reaction string method` をご参照ください。

## 計算手法 ​

`ReactionStringFeature` のアルゴリズムは `CI-NEB` に近いアルゴリズムですが、TSを高精度かつ自動的に計算するため、String法をベースにいくつかの改善を行いました。

主な利点は以下のようになっています。

  * CI-NEBと異なり、イメージ間距離を指定します。イメージ数はイメージ間距離を用いて自動的に決定され、アルゴリズムの途中でも変化します。NEBの失敗の典型的な原因として、イメージ間隔が広すぎることがありますが、これが起こりにくくなります。
  * 反応経路最適化とTS最適化が自動で進行します。
  * Dimer法やSellaに近い精度でTSの最適化を行います。 `fmax_ts` と `fmax_rd` を指定することにより各TSの収束条件を細かく決めることができます。
  * 障壁が複数あるような反応経路も適切に取り扱うことが出来ます。全ての障壁のTSを求める、一番高い障壁のTSだけ求める、興味のある障壁のTSのみを求めるなど、柔軟に計算条件を設定できます。
  * kink(反応経路のねじれ)を解消する機能があるため、kinkに由来した計算失敗が起こりにくいです。初期値が余程悪くない限り、多くの場合に最後まで失敗せずに計算が進行します。



### String法について ​

String法では、現在の経路を補間した線の上に `dx_rp` 以下の間隔で並ぶようにイメージを配置します。力の経路垂直方向成分を計算し、その方向に配置したイメージを移動させることで経路全体を最適化します。

`ReactionStringFeature` のアルゴリズムはNEBとString法を組み合わせたアルゴリズムです。イメージ間隔が大きくなりすぎないように反応経路を最適化します。

### TSの精度について ​

CI-NEBでTSの精度を上げにくい原因は以下の2つがあります。

  * CI-NEBでTSの収束条件を厳しく設定するためには、イメージの数を大きく設定しなければならず、それに伴い計算量が大きくなってしまう。
  * CI-NEBではTSの周りのイメージ間隔が反応経路の他のイメージ間隔と同じくらい大きいため、反応経路方向の正確な推定が困難である。



そこで、本実装ではTSの精度を上げるためTS付近のイメージ間隔を非常に小さく (dx_ts) 取るようになっています。 これによってTSの収束条件の厳密化および反応経路方向の正確な推定が可能になり、高精度なTSを得ることができます。

## 変数の説明 ​

引数については [ReactionStringFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature) 、返り値については [ReactionStringFeatureResult](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult) をご参照ください。

また、引数では共通するsuffixやprefixが使われています。以下にこれらの定義を示します。

#### 変数のsuffix ​

suffix| Long form| 意味  
---|---|---  
rp| reaction path| EQおよびTSを除く、反応経路を構成する全ての点  
ts| transition state| 反応経路中でエネルギーが極大値となる点(TS)  
rd| rate determining step| 一連の反応経路の中で最もエネルギーの高いTS  
is| initial state| 最初の構造  
fs| final state| 最後の構造  
eq| equilibrium structure| エネルギーが極小となる点(EQ)  
  
#### 変数のprefix ​

prefix| 意味  
---|---  
fmax| 収束条件を示します。各原子にかかるForceの最大値がfmax以下であるかどうかで収束を判定します。単位はeV/Åです。  
dx| イメージ間の距離を示します。単位はÅです。  
optimize| ISやFSを最適化するかどうかを決定します。最適化の経路も反応経路に追加されます*。デフォルトでTrueです。  
logfile| ログを書き込む先を指定します。ASEのoptimizer等で採用されている形式に準拠しています。’-‘で標準出力に出力します。Noneで出力しません。文字列を指定した場合、ファイルパスと解釈されます。  
  
* 最適化途中の構造からいくつかサンプルし、イメージとして反応の前後に追加します。構造最適化による経路の不連続なジャンプを防ぎます。

## ReactionStringFeature objectの引数の型 ​

`ReactionStringFeature` objectの引数は `List[List[MatlantisAtoms | AseAtoms]]` の型を持っています。 これは返り値と同じ型を意図した型となっており、複数の障壁からなる反応経路の概形を入力として取ることが出来るようになっています。 入力を収束した反応経路に近くすることで計算時間を大きく削減することが出来ますが、ISとFSの線形補完から最適化を開始したい場合、 `[[IS, FS]]` という形でISとFSのみ指定することも出来ます。

## ログ出力について ​

logfile_rp 引数にファイル名を指定すると、そこに反応経路とTSの最適化に関するログが記されます。 ログの大部分はOptimizerのログですが、時折”2-3:4/5”のような行が出ることがあります。 このログ(“2-3:4/5”)は以下を意味しています。

  * 経路全体に関して2周目の最適化を行っている。
  * 現在、経路は5つの山に分割されている
  * 現在最適化中の山はISから数えて4番目に位置し、全体の中で3番目にエネルギーが高い山である。



## うまくいかない場合 ​

### 計算がいつまで経っても終了しない ​

まずは経路が長すぎないかを確認してみてください。経路が長いとimage数が増えて時間がかかります。 本当に欲しい経路よりも長い経路を与えてしまっている場合、初期値を改善することで計算時間が短くなることがあります。

また、fmax_tsが厳しく(小さく)設定されていると、TSを収束させるのに時間がかかります。fmax_tsを大きく設定したり、interesting引数*を使って不要なTSにかける時間を削減してみてください。

* 興味のある山かそうでない山かを幾何的情報などから判定する関数を定義し、interesting引数に指定することができます。興味がない山について最適化を行わなくなり、高速化を図ることができます。使い方についてはMatlantis Examplesの `Reaction string method` および [API reference](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature) をご覧ください。

### 解が思った経路にならない ​

反応経路の初期値の作り方を工夫してみてください。例えばase.constraints.FixInternalを用いてscan計算を行うと直線補間より良い初期値を作成することが出来ることがあります。

### TSが二次以上の鞍点になってしまった ​

高次の鞍点付近ではforceが小さくなるため誤って収束条件を満たしてしまうことがあります。これは FIRE 等のlocal optimizerが鞍点に収束してしまうケースがあるのと同様の現象です。 local optimizerにおいて `fmax` を小さくすると鞍点で止まりにくくなるのと同様に、 `fmax_ts` や `fmax_rd` を小さくすることが高次の鞍点での停止の防止に役立ちます。
