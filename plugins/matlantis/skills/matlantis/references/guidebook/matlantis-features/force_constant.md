# ForceConstantFeature ​

## はじめに ​

[ForceConstantFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.phonon.force_constant.ForceConstantFeature.html#matlantis_features.features.phonon.force_constant.ForceConstantFeature) は結晶の力定数を求めるのに用います。 力定数の値から、フォノンの分散や状態密度、いくつかの熱的な特性などを求めることができます。

## 計算手法 ​

ForceConstantFeature では [VibrationFeature](./vibration.html) と同じく、結晶の力定数を調和近似により計算します。 結晶の周期性を取り入れるため、ユニットセルを3方向に繰り返したしたスーパーセルを作ります。(たとえば  倍など) このスーパーセルのポテンシャルエネルギーは次のように近似されます。

ここで、  ,  はx, y, zを、  ,  は原子のインデックスを、 そして  ,  はユニットセルのインデックスを表します。  は原子の平衡位置です。  は  番目の原子の  番目のユニットセルでの  方向の座標を表します。

上式のポテンシャルエネルギーの2階微分  は力定数と呼ばれ、 これを  で表わします。 `matlantis-features` では、力の微分を有限差分法から計算します。

ここで、  は、 セル  内の原子  がその平衡位置から  方向に距離  だけ動いたときの、 セル  内の原子  に作用する力の  方向成分を示しています。

結晶の力定数は  をユニットセルあたりの原子数、  をユニットセルの総数とすると、  のテンソルになりますが、 結晶の並進対称性から  だけ考慮すればよいことがわかります。 たとえば、ユニットセルを  倍した場合、 得られる力定数のサイズは  となります。 最初のスライス、  はセル (0,0,0) の自分自身との相互作用を表します。 第2スライス、  はセル (0,0,0) とセル (0,0,1) との相互作用を表します。

力定数と実験から得られる物性値を直接比較することはできませんが、Post-calculation Featuresを用いて力定数をさらに解析することで、フォノンのバンドや状態密度および一部の熱力学特性を得ることができます。

## Tips ​

  * `ForceConstantFeature` の入力として使用できるのは安定構造のみです。 原子位置をあらかじめ最適化してからご使用ください。

  * 入力構造にはprimitive cellの使用が推奨されます。 conventional cellでも動作しますが、結果から得られるフォノンのバンド構造が煩雑になってしまいます。

primitive_matrix 引数を指定することで手動または自動的に入力構造をprimitive cellに変換することができます。

default
        
        # Manually set the basis vector of the primitive cell
        force_constant = fc(
          atoms_Si,
          primitive_matrix=np.array(
              [
                  [0.00, 0.25, 0.25],
                  [0.25, 0.00, 0.25],
                  [0.50, 0.50, 0.00],
              ]
          ),
          prec=1e-3,
        )
        
        # automatically detect the primitive cell
        force_constant = fc(atoms_Si, primitive_matrix="auto", prec=1e-3)

この例の前半では、primitive cellの基本格子ベクトルを

として指定しています。 ここで、 , ,  は入力構造の基本格子ベクトルです。 注意点として、 手動で設定する場合には `primitive_matrix` を正確に入力するようにしてください。 もともとのcellにおいて折りたたまれる対象となる原子は、折りたたみ後の座標がprimitive cell内の原子のどれかと重なっている(距離が `prec` 以下である)必要があります。そうでない場合、例外 `ValueError` が生じます。

一方、 `primitive_matrix="auto"` を指定した場合、 `spglib` により対称性から自動的にprimitive cellを判定します。

  * パラメータ `supercell` を指定することで、ユニットセルを繰り返してスーパーセルを作成することができます。 スーパーセルは、力定数における無視できないセル間相互作用をすべて考慮できるように、十分に大きくする必要があります。 とるべきスーパーセルの大きさは材料に強く依存するため決定するのは困難ですが、 大きさを変化させながら収束性の確認を行なうことが望ましいです。

  * パラメータ `delta` により有限差分における変位の大きさを指定することができます。 PFPの出力に含まれる微細なノイズの影響を抑えるため、比較的大きな値（0.1程度）に設定することを推奨します。




## Post-calculation Features ​

  * PostPhononBandFeature
  * PostPhononDOSFeature
  * PostPhononThermochemistryFeature
  * PostPhononModeFeature



## 使用例 ​

ランチャーページから起動できる `Matlantis-features: Phonon` をご覧ください。
