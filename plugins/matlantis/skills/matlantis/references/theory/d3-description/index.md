# DFT-D3補正の概要 ​

## DFT-D3補正について ​

Matlantisで使われるニューラルネットワークポテンシャル(NNP)であるPFPは、密度汎関数法(DFT)によって作られたデータセットのエネルギーを再現するように学習しています。DFTは幅広い物質に対して適用され、また成功を収めていますが、それでもいくつかの現象はDFTの近似を超えた効果と関係していることが知られています。基底状態を取り扱うDFTの枠組みで扱うことが比較的難しく、かつ一部のシミュレーションで無視できない効果をもたらすものとして分散力(Van der Waals力)があげられます。分散力とは、電子状態が励起することで生じる双極子がもたらす相互作用によって生じる力です。

この分散力の効果を扱う手法として、広く使われているものにDFTD3[1, 2]があります。DFTD3補正は、明示的に分散力からの寄与分のエネルギーを計算し、DFT計算で求めたエネルギーに加算するモデルとなっています。幸いなことに、DFT計算本体と独立に計算できること、またDFT計算と比較して要求される計算コストが非常に軽いことから、D3補正はPFPに対しても有効な補正となっています。

## D3補正の実装とtorch-dftdについて ​

DFTD3のオリジナル実装は著者のグループが公開しています。[3] 多くのDFT計算手法でこのD3補正が利用できるようになっています。

ここでは、オリジナルとは別の実装であるtorch-dftd [4] について紹介します。torch-dftdは、DFTと比較して圧倒的に高速に動作するNNPの推論と組み合わせて使うことを考えて、高速化を志向したDFTD3の実装となっています。計算手法や結果自体はオリジナルの実装と同一になるようにしつつ、相互作用をPyTorchで実装することで、GPUなどの専用デバイスを利用した並列化の恩恵が受けられるようになっています。

torch-dftdについてのより詳しい解説および計算内容についての情報は、[解説記事](https://tech.preferred.jp/ja/blog/dft-dispersion-pytorch-oss/)および[実装](https://github.com/pfnet-research/torch-dftd)を参照してください。[4, 5]

## PFPにおけるD3補正 ​

PFPでは、上記のtorch-dftdをバックエンドとしたD3補正を利用できるようにしました。PFPにおけるD3補正はオプション扱いとなっており、計算モードとして`PBE_PLUS_D3`または`R2SCAN_PLUS_D3`を選ぶとPFPの出力に、PBE汎関数用またはr2SCAN汎関数用のD3補正の値を加算した結果を返します。r2SCAN汎関数用のD3補正パラメータについての情報は、[文献](https://pubs.aip.org/aip/jcp/article/154/6/061101/199772/r2SCAN-D4-Dispersion-corrected-meta-generalized)を参照してください。[6]

通常の系ではD3補正の計算時間は支配的ではないと考えられますが、それでも計算時間の多少の増加が予想されます。

なお、WB97XDモードには学習対象のDFT計算にすでに分散力補正が入っているため、D3補正を追加する計算モードは提供していません。

PFPで提供されるD3補正のパラメータは、以下のようになっています。

Parameter| Setting  
---|---  
xc| pbe or r2scan  
damping| bj  
old| False  
abc| False  
cutoff| 14 Å  
cutoff_smoothing| poly  
precision| FP32  
  
## 元素の対応状況 ​

PFPにおけるD3補正は、HからPuまでの94元素に対応しています。

## 参考文献 ​

[1] Stefan Grimme, Jens Antony, Stephan Ehrlich, Helge Krieg, J. Chem. Phys. 132, 154104 (2010).

[2] Stefan Grimme, Stephan Ehrlich, Lars Goerigk, Comp. Chem. 32, 7, 1456-1465 (2011).

[3] <https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/get_dft-d3>

[4] <https://tech.preferred.jp/ja/blog/dft-dispersion-pytorch-oss/>

[5] <https://github.com/pfnet-research/torch-dftd>

[6] Sebastian Ehlert, et al., J. Chem. Phys. 154, 061101 (2021).
