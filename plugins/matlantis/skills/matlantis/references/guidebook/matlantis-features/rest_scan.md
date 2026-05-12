# RestScanFeature ​

## はじめに ​

RestScanFeature はGaussianのscan機能に似た機能で、構造を変化させてトラジェクトリを得るための機能です。単純な分子や物質を逐次的に変化させてより複雑な構造を得たり、反応経路の初期値を得るのに使用することができます。 実際の使用例についてはMatlantis Exampleの `RestScanFeature` をご参照ください。

## 計算手法 ​

`RestScanFeature` のアルゴリズムは Gaussian の `Scan` に近いアルゴリズムですが、複数の変化を同時にもしくは逐次に適用できるようにいくつかの改善を行いました。 主な利点は以下のようになっています。

  * GaussianのScanと異なり、構造変化は不等式の形で与えます。例えば「0-1の距離がdÅ以下になるまで縮める」といった指示ができます。
  * 複数の指示を合わせて使うことができます。例えば、「0-1の距離がd1Å以下になるまで縮めながら3-0-1-2からなる二面角を120度回転させ、その後で1-2の距離をd2Å以上になるまで回す」といった指示ができます。
  * 指示が不等式の形で表されているため、変形の目的値が曖昧でもうまく動作します。厳密に何Åの距離にするといった指示は必要ありません。



計算手法

系の座標を x 、系のエネルギーを E(x) とします。   
x の全ての要素による偏微分を ∇ で表します。   
この時、系に働く力は以下のように記述されます。   


F(x)=−∇E(x)

i 番目のCollective Variable (CV) を Ci(x) とし、例えばその全ての値が Cei 以上になるまで変化させたいとします。   
この時、scanに用いる力 Fscan は理想的には以下のようになります。   


Fscan(x)=F(x)+∑i∈{i|Ci(x)<Cei}max(ci−F(x)⋅∇C(x)|∇C(x)|,0)∇Ci(x)|∇Ci(x)|

ここで、 ci は以下の条件を満たすように決められます。

ci=fmaximaxnorm(∇C(x)|∇C(x)|)

ここで maxnorm とは原子ごとにノルムをとり、その値を全ての原子間に関しては max をとる演算です。この演算はaseの一般的なoptimizerの収束条件に用いられるものと同じです。

値を Cei 以下にしたいものに関しては不等号の向きや項全体の符号が逆になります。   


この Fscan を用いて構造を変化させていきます。

dxds=Fscan(x)

また、理想的にはscanの式はこのような式になりますが、値が急激に変わると微分方程式が不安定になるため Cei にまつわる不等式に関しては少しぼかしています。 このぼかす範囲をexceedとしています。   


ScanDistanceやScanDihedralのsuper classであるScanRestraintはdirection, destination, exceed, fmaxというパラメーターを持っています。directionはCVを増やすか減らすかを決めるパラメーターです。destinationは目的値を決めるパラメーターです。exceedは上記のexceedです。地形によってはdestinationをexceed程度超えることがあります。fmaxは上記外力の大きさを決めるパラメーターです。

Scanを同時または逐次的に実行するため、Scanの引数は `List[List[ScanRestraint]]` の形で与えられます。例えば以下の例は57-74の距離をrより小さくし、その後で54-79の距離を小さくしながら74-79の距離を大きくすることを指示しています。

default
    
    
    r = covalent_radii + 0.1
    scans = [
        [ScanDistance(indices=(57, 74), direction=-1, destination=r[57] + r[74])],
        [
            ScanDistance(indices=(54, 79), direction=-1, destination=r[54] + r[79]),
            ScanDistance(indices=(74, 79), direction=1, destination=(r[74] + r[79]) * 1.5),
        ]
    ]

## ユーザー定義CV ​

組み込みで `ScanDistance` と `ScanDihedral` を用意していますが、 `ScanRestraint` を継承して簡単に新しいCVを定義することができます。 詳しくは `ScanDistance` の実装を参考にしていただきたいですが、最低限 `get_colvar` と `get_colvar_forces` を実装すると新しいCVを動かすことができます。 また、 `index_shuffle` は一般的な aseの constraintと同じ流儀で記述して構いません。
