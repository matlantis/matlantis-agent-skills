# PFP v8.0.0の検証 ​

## 概要 ​

  * PFP v8.0.0で提供されたR2SCANモードに関する検証を行いました。 
    * 結晶の形成エネルギーの対実験値での再現性がPBEモードよりも良くなりました。
    * 分子系ベンチマーク(Wiggle150)において、PBE_PLUS_D3モードおよびいくつかのOSSのMLIPと比較して、R2SCAN_PLUS_D3モードがより高い再現性を示しました。
    * 分子系ベンチマーク(GMTKN55)において、PBE汎関数を用いたDFT結果よりも良い性能を示しました。
    * R2SCANモードでの液体密度の実験値の再現性が、PBEモードよりも改善されました。
    * R2SCANモードでの水の再現性においては、先行研究のDFT計算結果とよく一致することを示しました。分散力補正を加えたR2SCAN_PLUS_D3モードでは液体の水の密度を過大評価する傾向がありますが、同様の現象がr2SCAN汎関数の元となったSCAN汎関数を用いたDFT計算でも報告されています。
  * PBEモードに関しては、従来と同様の検証に加えてMDR Phonon Benchmark(熱物性の再現性の検証)を追加しました。 
    * 水素単結晶に関する外れ値があるものの、全般的にPFP v7.0.0比で精度が改善していることを確認しました。



前回のPFP v7.0.0の検証結果については、[PFP Validation 2024.09](/api/docs/ja/pfp-validation/2024.09/index.html) を参照してください。

## R2SCANモードの検証 ​

ここではPFP v8.0.0で新たに提供されたR2SCANモードの検証を行います。なお、PFP v8.0.0ではランタノイド、アクチノイド、クラスター、錯体、分子吸着構造などが r2SCAN の学習データセットに含まれていないため、PFP v8.0.0のR2SCANモードではこれらの物質の物性予測において再現性が限定される可能性があります。

### 結晶の形成エネルギー ​

ここでは結晶構造の標準形成エンタルピーの実験値のデータセットを用いた検証を行います。対象のデータセットとして、Matminerに掲載されているexpt_formation_enthalpy_kingsburyデータセット [1]を用います。このデータセットはKubaschewski [2]、NIST JANAF [3]、Kim et al [4] [5]に記載されているデータから構成されています。

データセット内のそれぞれの系に対応するMaterials Project [6]の構造に対して、PFP v8.0.0のR2SCAN、および、PBEあるいはPBE_Uモードを用いて構造最適化を行った後、形成エネルギーを算出しました。 ここでは検証対象をPFP v8.0.0のR2SCANモードの対応元素のみを含む738構造に限定しています。 PBEモードに関しては、補正を考慮しない場合と考慮する場合をそれぞれ検証しました。 補正を考慮する場合は、PBEおよびPBE_Uモードを用いて対象の系に応じて適切にHubbard U補正の有無を切り替え、また形成エネルギーの算出時にMaterials Projectで用いているものと同様の実験値に近づけるための補正 [7]を行いました。 また、参照値としてMaterials Project [6:1]から取得したR2SCAN汎関数を用いたDFT計算の結果を用いて、同様に形成エネルギーを求めました。 各モードでの推論結果およびMaterials Projectの結果を実験値と比較した図を以下に示します。 なお、温度、圧力、零点振動を考慮していない形成エネルギーと、標準形成エンタルピーとを比較していることに注意してください。

R2SCAN モード| PBE+PBE_Uモード(補正を考慮)| PBEモードのみ| 参考：Materials Project (R2SCAN)[6:2]  
---|---|---|---  
| | |   
  
対実験値でのMAEの値を比較すると、PBEモードとPBE_Uモードを組み合わせて補正を加えた場合よりも、補正なしのR2SCANモードの方が良い再現性を示しています。 また、Materials ProjectのR2SCAN汎関数を用いたDFT結果とのMAEの値は0.080 eV/atom であり、PFP v8.0.0のR2SCANモードのMAEはそれと同程度の再現性を示すことが分かりました。

### Wiggle150: 分子系ベンチマーク ​

Wiggle150 [8]は、アデノシン、ベンジルペニシリン、エファビレンツという3つの有機分子についての、合計150個の非常にひずんだ立体配座からなるベンチマークセットです。 これらの構造は、結合長、結合角、二面角が通常の安定構造から大きく外れており、エネルギーも非常に高いという特徴があります。 このような特徴は、特に化学反応の経路探索などで重要となります。

このベンチマークセット中には参照データとして、各構造についてDLPNO-CCSD(T)で計算されたエネルギーが含まれています。DLPNO-CCSD(T)は、化学的精度（1 kcal/mol）を達成できるとされる高精度な手法です。この値との比較は、実験値との比較と近しいものになります。各構造についてPFP v8.0.0の各モードでエネルギーを計算し、基底状態からのエネルギー差を求め、DLPNO-CCSD(T)による参照データと比較しました。特にv8.0.0のR2SCAN_PLUS_D3モードとPBE_PLUS_D3モードによる計算結果を下に示します。

R2SCAN_PLUS_D3| PBE_PLUS_D3  
---|---  
|   
  
DLPNO-CCSD(T)での値に対するMAEを見ると、R2SCAN_PLUS_D3モードでの相対エネルギーの再現性がPBE_PLUS_D3モードを上回ることがわかります。

同様に、他のPFPのバージョンや計算モード、及びいくつかのOSSのMLIP [9] [10] [11] [12] [13]を用いて相対エネルギーを計算し、DLPNO-CCSD(T)での値に対するMAEを求めた結果を以下に示します。

この結果からWiggle150ベンチマークセットにおいて、OSSを含む他のモデルと比較しても、PFP v8.0.0 のR2SCAN_PLUS_D3モードがより高い再現性を示すことが分かります。

### GMTKN55: 分子系ベンチマーク ​

GMTKN55 [14]は分子系を対象としたデータセットで、さまざまなサイズの系における基本的な特性、反応エネルギー、反応障壁や相互作用などについての各種サブデータセットから構成されています。 今回はこのデータセットを用いて、PFP v8.0.0とDFTの分子系に関する性能を比較します。

現在のPFPでは荷電した系は扱えないため、そのような系はデータセットから除外しています。また単原子分子についてもエネルギーの基準点の相違のため対象外としました。

具体的な検証対象および用いた38個のサブベンチマークデータセットは以下の通りです。

| 検証対象| サブベンチマークデータセット  
---|---|---  
basic + small| 基本的な特性と小規模な系の反応エネルギー| W4-11, ALKBDE10, YBDE18, AL2X6, HEAVYSB11, NBPRC, ALK8, G2RC, BH76RC, FH51, TAUT15, DC13  
iso. + large| 大規模な系の反応エネルギーと異性化反応| MB16-43, DARC, RSE43, BSR36, CDIE20, ISO34, ISOL24, C60ISO  
barriers| 反応障壁| BH76, BHPERI, BHDIV10, INV24, BHROT27, PX13, WCPT18  
intermolecular NCIs| 分子間非共有結合相互作用| RG18, ADIM6, S22, S66, HEAVY28, WATER27, CARBHB12, PNICO23, HAL59, IL16  
intramolecular NCIs| 分子内非共有結合相互作用| IDISP, ICONF, ACONF, Amino20x4, PCONF21, MCONF, SCONF, BUT14DIOL  
all NCIs| 分子間および分子内相互作用| intermolecular NCIsとintramolecular NCIsのサブベンチマークデータセット  
All| GMTKN55のうちの対象38個| 上記のすべてのデータセット  
  
上記ベンチマークを用いた検証結果を以下の表にまとめました。 評価指標には加重平均絶対偏差(WTMAD-2 [kcal/mol])を用います。これはより低い値の方が良い再現性となる指標です。 DFT結果についてはGMTKN55のwebサイト [14:1]に記載された値を用いて、今回対象としたデータセットに関してWTMAD-2を計算した結果です。

model_version/functional| calc_mode| basic + small| iso. + large| barriers| intermol. NCIs| intramol. NCIs| all NCIs| All  
---|---|---|---|---|---|---|---|---  
PFP v8.0.0| R2SCAN_PLUS_D3| 5.47| 7.87| 15.38| 10.60| 9.82| 10.21| **9.28**  
|  R2SCAN| 5.48| 7.87| 15.53| 13.52| 11.32| 12.42| 10.51  
| PBE_PLUS_D3| 7.24| 11.19| 18.72| 13.26| 10.73| 11.99| 11.54  
| PBE| 7.21| 15.28| 19.12| 19.69| 19.15| 19.42| 15.43  
| WB97XD| 8.66| 9.91| 14.83| 20.12| 18.47| 19.29| 14.23  
PFP v7.0.0| PBE_PLUS_D3| 7.47| 13.19| 17.99| 12.24| 10.67| 11.45| 11.61  
PBE-D3(BJ) [14:2]| N/A| 6.33| 11.49| 17.80| 11.23| 9.92| 10.57| 10.62  
PBE [14:3]| N/A| 6.46| 15.76| 16.10| 17.72| 19.65| 18.69| 14.58  
SCAN-D3(BJ) [14:4]| N/A| 5.12| 6.92| 14.30| 9.20| 6.85| 8.02| 7.94  
SCAN [14:5]| N/A| 5.05| 8.13| 13.84| 11.73| 8.29| 10.00| 8.91  
wB97X-D3(0) [14:6]| N/A| 2.82| 7.95| 4.69| 5.03| 4.64| 4.83| 4.82  
  
この表から、PFP v8.0.0のR2SCANモードが分子系の各特性について、PBEモードと比較して良い再現性を示すことが分かります。 またD3補正を加えたR2SCAN_PLUS_D3モードを用いることで、さらに再現性が向上することが確認できます。 特にGMTKN55の対象データセット全体に関する指標(All)を見ると、PFP v8.0.0のR2SCAN_PLUS_D3モードの性能は、PBE汎関数およびD3補正を用いたDFT計算の精度を上回っていることが分かります。

(補足) PFPのWB97XDモードはwB97X-D3(0)汎関数でのDFT結果に比べて数段劣った再現性を示しています。これはWB97XDモードの学習に用いたデータセットの制限のためであると考えられます。

### R2SCANモードでの液体の密度の再現性 ​

このベンチマークではR2SCANモードでの液体の密度の再現性を検証します。 PFPの計算モードにはPFPのR2SCAN_PLUS_D3モードを使用しました。 その他の詳細な検証方法は後述のPBEに関する検証の中の、液体の密度の再現性を参照してください。

この検証により得られた液体の密度の再現性を以下の図に示します。R2SCANモードの方がPBEモードでの結果よりもMAEの値が改善していることが分かります。

R2SCAN_PLUS_D3| PBE_PLUS_D3  
---|---  
|   
  
一方で、水の密度に関してはR2SCAN_PLUS_D3モードでは実験値を大きく過大評価しています。水の密度に関しては、次節でさらに詳しく検証します。

各物質に対する詳細な値については有機液体密度のベンチマークの詳細な値を参照してください。

### R2SCANモードでの水の動径分布関数、密度の再現性 ​

液体構造の水を対象として、R2SCANモードの水の動径分布関数および密度の再現性をDFT計算結果、PBEモードおよび実験結果と比較します。

#### 計算条件 ​

MDシミュレーションの条件は以下の通りです。

  * 分子数: 64分子の水
  * 初期構造: 一辺12.42 Åの立方体シミュレーションボックスに分子を配置し、体積密度1 g/cm³のバルク水構造を作成
  * 動径分布関数の計算: 
    * 350 KのNVTアンサンブル
    * 時間刻み: 0.5 fs
    * 初期平衡化: 1 fs刻みで1 psのMD
    * 30 psのMDトラジェクトリから平均動径分布関数を計算
    * 動径分布関数の計算はpymatgenを使用
  * 密度計算: 
    * 300 KのNVTアンサンブルで1 psのMDによる初期平衡化
    * 300 KのNPTアンサンブルで30 psのMD
    * 時間刻みは動径分布関数の計算と同様
    * 最後の20 psのMDトラジェクトリから平均密度を計算



#### 動径分布関数 ​

核量子効果(Nuclear Quantum Effects, NQE)が水分子の運動に大きな影響を与えることが知られており、水素核を古典的に扱った場合、動径分布関数のピーク強度が増大することが報告されています[15]。そこで、今回の検証では、水におけるNQEを模擬するため、比較対象のDFT計算と同様に、温度を約50K上昇させて350 KのNVTアンサンブルで行っており、室温(298K)における実験値(X-Ray Diffraction)と比較しています。R2SCANモードとDFT計算結果[16]、PBEモードおよび実験値[17]との比較を下図に示します。

R2SCANモードは、DFT計算結果と良い一致を示しています。PBEモードと比較すると、ピーク強度などの特徴が改善しており、実験値との一致も良くなっています。水素結合した隣接原子間の相関は、gOO(r)の第1ピークとgOH(r)の第2ピークに反映されており、R2SCANモードは水素結合を記述する能力がPBEモードよりも向上していると考えられます。また、R2SCANモードとR2SCAN_PLUS_D3モードの動径分布関数の違いは無視できる程度でした。これは、最近報告されたr2SCAN汎関数の元となったSCAN汎関数と、これにvdW補正を加えたSCAN+rVV10汎関数との動径分布関数の比較結果[18]とも整合しており、vdW相互作用の取り扱いが局所的な構造に大きな影響を与えないことを示唆しています。

#### 密度 ​

300Kにおける水の密度について、PFPとDFT計算結果[18:1]との比較を示します。

Model| Liquid Water density (g/cm³)  
---|---  
PFP-R2SCAN| 1.01±0.02  
PFP-R2SCAN_PLUS_D3| 1.11±0.02  
PFP-PBE| 0.90±0.01  
PFP-PBE_PLUS_D3| 0.99±0.02  
DFT-SCAN[18:2]| 1.05  
DFT-SCAN+rVV10[18:3]| 1.16  
  
R2SCANモードでは、vdW補正を考慮しない場合に実験値の密度をよく再現している一方で、vdW補正を考慮した場合(R2SCAN_PLUS_D3)では、密度を大幅に過大評価していることが分かります。先行研究のSCAN, SCAN+rVV10のDFT計算結果と同じ傾向を示していることから、r2SCAN汎関数の特徴としてvdW補正が水分子間相互作用に対して過度に強く働いており、R2SCAN_PLUS_D3モードではその特徴を再現していると考えています。

### COD構造によるベンチマーク ​

PFPの現実世界に存在する多様な物質に対する再現性を確認するために、Crystallography Open Database (COD) [19]を用いた検証を行いました。なお、ここで用いるCOD構造を用いたベンチマークデータセットはPFPの学習には利用していません。

#### ベンチマークデータセット ​

CODから、PFP v8.0.0のr2SCANモードでサポートする70元素のみで構成される結晶構造のうち、単位胞の体積が2200Å³以下のものを取得しました。ただし、各サイトの占有率が0.99未満のものや、対称性の指定に問題があるもの等は除外しています。さらに、学習データセットに含まれる結晶構造を除外しています。ただし、これによって除外されなかった結晶構造についても、学習用データセットに同じ結晶構造が含まれている場合もあります。例えば、Materials Project [6:3] に同じ結晶構造が登録されており、この構造を学習用データセットでも利用している場合があります。これらの結晶構造に対して、構造最適化とサイト微小変位を加えた構造を用意し、DFT計算を行いベンチマークデータセットを生成しました。DFT計算の計算条件は`R2SCAN`モードに相当するものです。`R2SCAN` モードの計算条件に関しては、[About PFP](/api/docs/ja/pfp-description/2025.07/index.html)を参照してください。

##### 構造最適化 ​

比較対象とするDFT条件での構造を得るために、格子定数（長さと角度）を含めた構造最適化をVASPで行いました。構造最適化の収束条件は、forceが0.02eV/Å以下としました。データセットには、構造最適化後の構造のみが含まれます。

##### 微小変位 ​

安定構造周辺の局所的なエネルギー曲面の再現性を確認するため、構造最適化後に微小変位を加えました。詳しくは、PFP論文[20]のSupporting Information, NOTE 10にある“site position displacement”を参照ください。ただし、変位の大きさは論文中よりも小さくしています。

#### エネルギーと力の再現性 ​

微小変位を加えた有機結晶、錯体結晶、無機結晶に対するエネルギー・力、および、微小変位を加える前後でのエネルギー・力の変化をDFT結果と比較しました。なお、有機結晶、錯体結晶、無機結晶の分類については後述のPBEモードに関する検証の中の、エネルギーと力の再現性と同様です。

y-yプロットではエネルギーは1原子あたりで比較し、力は各原子にかかる力のxyz成分を独立に比較しています。ヒストグラムでは、微小変位を加える前後でのエネルギーと力の変化をDFTとPFPとで比較し、その誤差をプロットしています。

| エネルギー| 力| 微小変位を加える前後でのエネルギーの変化| 微小変位を加える前後での力の変化  
---|---|---|---|---  
有機結晶| | | |   
錯体結晶| | | |   
無機結晶| | | |   
  
後述のPBEモードに関するエネルギーと力の再現性と比較すると、無機結晶についてはPBEモードとMAE比で同程度である一方で、有機結晶と錯体結晶についてはR2SCANモードで悪化していることが分かります。この問題については将来のバージョンでの改善を予定しています。

### 無秩序構造に関する力とエネルギーの再現性 ​

ここでは無秩序構造に関する力とエネルギーの再現性を検証します。今回は簡易的に、PFP v8.0.0の学習データセットの一部を用いて検証を行います。そのため、ここでの検証は限定的なものとなります。今後、改善のために検証用のデータセットの作成を検討しています。

ランダムに選んだ初期構造に対し高温でMDシミュレーションを実行し、無秩序構造を作成しました。詳しくは、PFP論文[20:1]のSupporting Information, NOTE 10にある“disordered”を参照ください。ただし、温度や元素の選択対象は論文とは異なる場合もあり、特に元素に関してはR2SCANモードがサポートしている70元素を対象にしています。このデータセットに対し、PBEモードおよびR2SCANモードでの力とエネルギーの再現性を検証し、両者を比較しました。

y-yプロットではエネルギーは1原子あたりで比較し、力は各原子にかかる力のxyz成分を独立に比較しています。 ヒストグラムでは1原子あたりのエネルギーの誤差の絶対値、および、各原子にかかる力の誤差のL2ノルムを示しています。

| R2SCANモード| PBEモード| ヒストグラム  
---|---|---|---  
エネルギー| | |   
力| | |   
  
PFP v8.0.0のR2SCANモードのDFT結果の再現性は、PBEモードのDFT結果の再現性と同程度であり、学習が適切に行われていることが分かります。なお、ここでの検証は学習データセットを用いた簡易的なものであり、あくまでも限定的な検証であることに注意してください。

## PBEモードの検証 ​

### COD構造によるベンチマーク ​

PFPの現実世界に存在する多様な物質に対する再現性を確認するため、Crystallography Open Database (COD) [19:1]を利用して再現性の検証を行いました。なお、ここで用いるCOD構造を用いたベンチマークデータセットはPFPの学習には利用していません。

#### ベンチマークデータセット ​

CODから、PFP v7.0.0以降でサポートする96元素のみで構成される結晶構造のうち、単位胞の体積が2200Å³以下のものを取得しました。これはv3.0.0から使用していた72元素対応のベンチマークデータセットを拡張したものです。ただし、各サイトの占有率が0.99未満のものや、対称性の指定に問題があるもの等は除外しています。これらの結晶構造に対して、構造最適化とサイト微小変位を加えた構造を用意し、DFT計算を行いベンチマークデータセットを生成しました。DFT計算の計算条件は`PBE` モードに相当するものですが、構造最適化を安定させるためSCFの収束条件(EDIFF)を10-7とより厳しくしています。`PBE` モードの計算条件に関しては、[About PFP](/api/docs/ja/pfp-description/2025.07/index.html)を参照してください。

PFP v5.0.0以降では、CODに含まれる結晶構造の一部を学習データセットとして用いています。このため、本検証にあたってはこれらの結晶構造を除外しています。ただし、これによって除外されなかった結晶構造についても、学習用データセットに同じ結晶構造が含まれている場合もあります。例えば、Materials Project [6:4] に同じ結晶構造が登録されており、この構造を学習用データセットでも利用している場合があります。

今回の検証ではベンチマーク実施速度の向上のため、CODベンチマークデータセットのサイズを縮小しています。これによるベンチマーク結果への影響がないことを確認しました。

##### 構造最適化 ​

比較対象とするDFT条件での構造を得るために、格子定数（長さと角度）を含めた構造最適化をVASPで行いました。構造最適化の収束条件は、forceが0.03eV/Å以下としました。 データセットには、構造最適化後の構造のみが含まれます。

##### 微小変位 ​

安定構造周辺の局所的なエネルギー曲面の再現性を確認するため、構造最適化後に微小変位を加えました。詳しくは、PFP論文[20:2]のSupporting Information, NOTE 10にある“site position displacement”を参照ください。ただし、変位の大きさは論文中よりも小さくしています。

#### エネルギーと力の再現性 ​

微小変位を加えた構造を対象として、エネルギーと力の再現性をPFP v8.0.0とPFP v7.0.0で比較しました。比較しやすくするため、プロット範囲を固定しており、範囲外にも比較的少数のデータ点が存在することに注意してください。

##### 有機結晶 ​

H, C, N, O, P, S, F, Cl, Br, Iのみで構成される構造を有機結晶とみなして、比較をしました。 y-yプロットではエネルギーは1原子あたりで比較し、力は各原子にかかる力のxyz成分を独立に比較しています。 ヒストグラムでは1原子あたりのエネルギーの誤差の絶対値、および、各原子にかかる力の誤差のL2ノルムを示しています。

| v8.0.0| v7.0.0| histogram  
---|---|---|---  
エネルギー| | |   
力| | |   
  
PFP v8.0.0では有機結晶に関して、エネルギー・力ともにMAEが改善していることがわかります。

微小変位を加える前後でのエネルギーと力の変化をDFTとPFPとで比較し、その誤差をヒストグラムとしてプロットしています。

エネルギー| 力  
---|---  
|   
  
微小変位を加える前後のエネルギーと力の変化は、v8.0.0とv7.0.0で同等です。

##### 錯体結晶 ​

以下の条件をすべて満たすものを錯体結晶とみなしました。

  * H, C, N, O, P, S, F, Cl, Br, Iのうちの2元素以上を含むこと
  * H, C, N, O, P, S, F, Cl, Br, I以外の元素を1元素以上含むこと
  * H, C, N, O, P, S, F, Cl, Br, Iで原子数の8割以上を占めること



有機結晶と同様に、微小変位を加えた錯体結晶に対するエネルギー・力、および、錯体結晶に微小変位を加える前後でのエネルギー・力の変化を比較しました。

| v8.0.0| v7.0.0| histogram  
---|---|---|---  
エネルギー| | |   
力| | |   
エネルギー| 力  
---|---  
|   
  
有機結晶に対する結果と同様、v8.0.0では錯体結晶に対するエネルギーと力の再現性が改善しています。

##### 無機結晶 ​

有機結晶と錯体結晶のいずれにも分類されない構造を、無機結晶とみなし、同様に比較しました。

| v8.0.0| v7.0.0| histogram  
---|---|---|---  
エネルギー| | |   
力| | |   
エネルギー| 力  
---|---  
|   
  
有機および錯体結晶に対する結果と同様、v8.0.0では無機結晶に対するエネルギーと力の再現性が改善しています。

#### 体積と格子定数の再現性 ​

DFTによる構造最適化後の構造を初期構造とし、PFPで格子定数を含む構造最適化を行い、体積の再現性をv8.0.0、v7.0.0で比較しました。 PFP v3.0.0以降にサポートされていた72元素のみを含む構造に対する結果です。

`ase.Atoms`である`atoms`に対して、以下のコードを実行しています。
    
    
    from ase.optimize import LBFGS
    from ase.constraints import UnitCellFilter
    
    ecf = UnitCellFilter(atoms)
    opt = LBFGS(ecf)
    opt.run(fmax=fmax, steps=2000)
    

fmaxの値として有機結晶と錯体結晶では0.03、無機結晶では0.01を用いています。 有機結晶と錯体結晶ではDFT-D3補正を用いられることが多く、本検証にあたってはDFT-D3補正を用いた場合の体積を有機結晶と錯体結晶に対して比較しています。 PFPでは`PBE_PLUS_D3`を指定しており、DFT計算では`PBE`での条件に加えて`IVDW=12`を指定することでDFT-D3補正を有効にしています。 無機結晶に対してはDFT-D3補正を用いていません。

| 相対体積| 相対誤差  
---|---|---  
有機結晶| |   
錯体結晶| |   
無機結晶| |   
  
v8.0.0では、v7.0.0と比較して無機結晶に関する体積の再現性が低下していました。この結果だけでは必ずしもPFP自体の性能低下と断定することはできず、その原因については調査中です。 有機結晶と錯体結晶に関しては、v7.0.0と同等の再現性が得られています。

格子定数についての検証結果が以下です。

| 長さ| 角度  
---|---|---  
有機結晶| |   
錯体結晶| |   
無機結晶| |   
  
v8.0.0では、v7.0.0と比較して無機結晶に関する格子定数の再現性がやや低下していました。原因については調査中です。 その他に関しては、v7.0.0と同等の再現性が得られています。

#### 高密度構造における安定性検証 ​

ここでは実際よりも高密度な構造に対するポテンシャルエネルギー曲面(Potential Energy Surface, PES)の安定性を検証します。 このような構造は通常は存在しませんが、超高圧下でのシミュレーションや、結晶構造探索などにおける機械的な構造生成で出現することがあります。

高密度構造でのPESの不安定により、これらの構造に対する構造最適化が収束しない、あるいは、本来の局所安定構造から大きく離れた構造に収束してしまう、といった問題が生じます。 過度に高密度な構造に対して非常に低いエネルギーが返されると、たとえばより安定でエネルギーが低い構造を探索するというタスクである、結晶構造探索を破綻させてしまいます。

##### COD構造 ​

Crystallography Open Database (COD) [19:2]に含まれる、PFP v3.0.0以降でサポートする72元素のみで構成される結晶構造を対象としました。 これらの構造をDFTを用いて構造最適化を行った後、各構造の格子長を変えて、等方圧縮におけるPESの安定を検証しました。 構造とDFTによる構造最適化の詳細は、前述のベンチマークデータセットと構造最適化の節を参照してください。 PFPでは、構造最適化に使用されたDFT条件に相当する`PBE`を指定しました。

0.985=0.9039倍に圧縮した構造から開始し、ステップごとに格子長を0.98倍にすることを繰り返しました。 エネルギーが前のステップと比較して低下またはnanとなった場合、PESが不安定になったとみなしました。 ヒストグラムは、不安定となったステップの一つ前での相対格子長の分布を示しています。

有機結晶| 錯体結晶| 無機結晶  
---|---|---  
| |   
  
グラフの左側にヒストグラムが集中しているほど、高密度でのPESの安定性が高いと言えます。 有機結晶、錯体結晶、無機結晶の分類はエネルギーと力の再現性と同様に行っています。 特に有機結晶と錯体結晶について、v8.0.0はv7.0.0よりも高い安定性を示しています。無機結晶に関してはv7.0.0と同等です。

有機結晶の場合、v7.0.0では格子長を0.64倍まで圧縮すると2%の構造に対してPESが不安定になり始ます。 一方で、v8.0.0では2%の構造に対してPESが不安定になり始めるのは、格子長を0.56倍まで圧縮したときです。 無機結晶にはより多様な元素が含まれるため、ヒストグラムがより広がりを持っています。

### 液体の密度の再現性 ​

25種類の分子を対象として液体構造の密度の実験値との比較を行います。

#### 計算条件 ​

計算に使用した25種類の分子のリストは以下の通りです。

ethylene glycol, acetic acid, benzene, n-octane, dimethyl sulfoxide, water, pyridine, toluene, diethyl ether, isopropyl alcohol, chloroform, N,N-dimethylformamide, acetonitrile, methane sulfonic acid, trifluoroethanol, 1-methylnaphthalene, γ-Butyrolactone, xylene, propylene carbonate, tetrahydrofuran, n-hexane, ethyl acetate, perfluorooctane, hexamethylphosphoric triamide, dioxane

液体の構造のため、有限温度での分子動力学(MD)計算を行います。それぞれの分子について、以下の手順でシミュレーションを行いました。

  1. 系の原子数が400以上になるように分子を複製
  2. 300 KのNVTアンサンブルで1 psのMD
  3. 300 K、1気圧のNPTアンサンブルで400 psのMD
  4. NPTアンサンブルのtrajectoryのうち最後の20 psについて平均密度を計算



今回はすべての分子について上記の計算条件で計算しました。

細かい計算条件は以下の通りです。

  * PFPの計算モードは`PBE_U_PLUS_D3` (D3補正あり)
  * MDのtimestepは全体を通して1 fs
  * NPTアンサンブルでの体積変化は等方変形のみ許可



#### 密度の比較 ​

計算された密度と実験値との比較を以下に示します。実験値はPubChem[21]およびSigma-Aldrich[22]から取得しました。参考までに、PFP v7.0.0で同じ計算を行ったものを同様に示しています。

PFP v8.0.0| PFP v7.0.0  
---|---  
|   
  
v8.0.0ではv7.0.0に比べてMAEが改善しております。特に過去のバージョンで課題となっていた、perfluorooctaneに関する再現性が向上していることがわかります。 perfluorooctaneは比較的長鎖のパーフルオロアルカンです。

各物質に対する詳細な値については有機液体密度のベンチマークの詳細な値を参照してください。

### 相図の再現性によるベンチマーク ​

ここでは相図の再現性を確認するため、絶対零度における組成をふった相図(状態図)を描画し、DFT計算結果との比較を行いました。

#### ベンチマークデータセット ​

結晶構造のデータベースとして、Materials Project [6:5]に登録されている結晶構造を使いました。DFTの計算条件の違いを吸収するため、PFPデータセットの計算条件(`PBE`に相当するもの)により再計算を行いました。再計算では格子定数を含めた構造最適化を行っており、収束条件はforceが0.005 eV/Å以下としました。なお、今回はベンチマークの問題設定の都合上、これらの初期構造の多くはPFPの学習用データセットにも含まれています。

#### 構造最適化計算 ​

これらの結晶構造を初期構造として、PFPの`PBE`モードで構造最適化計算を行い、PFPにおける緩和構造を得ました。PFPの構造最適化計算にはASEのLBFGS法を使い、`fmax=0.01, alpha=300, steps=4000`としています。PFPの構造最適化計算においても、格子定数を含めた構造最適化を行っています。

最終的に得られた結晶構造は95347個です。

#### 相図の再現性 ​

##### Energy above hull ​

相図は、元素の組成と形成エネルギーの空間に様々な結晶構造を配置した図であり、特に安定に存在しうる結晶構造群からなる凸包(convex hull)として表現されます。用語の区別として、ここでは構造最適化計算で得られた局所安定な構造を緩和構造、相図を描いたときに凸包を構成する構造を安定構造と呼ぶことにします。ここでは相図の再現性の指標として、相図を描いたときの安定構造以外の結晶構造のエネルギー(energy above hull)の予測精度に着目しました。これは、相図における凸包が正しく描けていることに加え、類似の組成を持つ結晶構造間のエネルギーの差が正しく表現できているかどうかを測る指標になります。

具体的な手順は以下のようになります。まず結晶構造のデータに対し、1, 2, 3元系の相図として考えられるすべての元素の組み合わせにおいて相図を描くことを考えます。これをDFT計算で得られたエネルギーとPFPにおける構造最適化計算後のエネルギーの2パターンについてそれぞれ行います。得られたすべての相図に対して、その中のすべての構造についてenergy above hullを計算します。DFT計算の相図でenergy above hullが0より大きい構造(=安定構造以外の構造)を抽出し、DFTとPFPの2つの相図でenergy above hullを比較します。

得られた結果を以下に示します。比較のためにPFP v7.0.0の結果も示しています。x, y軸がそれぞれlogスケールになっていることに注意してください。PFP v8.0.0はv7.0.0と比較して、結晶構造のエネルギーの予測精度が改善しています。

PFP v8.0.0| PFP v7.0.0  
---|---  
|   
  
PFPとDFTとでのenergy above hull (Ehull)の差をプロットしたものが以下の図になります。 それぞれ、DFTでのenergy above hullが 1 eV/atom以下、0.01 eV/atom以下の構造のみに絞ってプロットしています。 両者ともv8.0.0ではv7.0.0と比較してMAEが改善しています。 また、0.01 eV/atom以下の構造に絞った場合の方がMAEが小さく、PFPではより安定な構造はより高い精度でenergy above hullを計算することができることがわかります。

Ehull ≤ 1 eV/atom| Ehull ≤ 0.01 eV/atom  
---|---  
|   
  
### MDR Phonon benchmark ​

ここでは、PFPにおける熱物性の再現性を確認するため、MDR Phonon benchmark[23]のワークフローに沿ってDFT計算との比較検証を行います。MDR Phonon benchmarkは、約10,000の化合物半導体の第一原理フォノン計算結果を検証データセットとして、汎用機械学習ポテンシャルの調和フォノン特性の予測能力を評価するベンチマークです。

#### データセット ​

PFPの`PBE`モードとの定量比較を行うため、MDR Phonon database[24]をPBE汎関数で再計算したもの[23:1]を対象のデータセットとして使用します。

#### 計算条件 ​

DFTで構造緩和した構造を初期構造として、PFPの `PBE` モードで構造緩和および力の計算を行います。構造緩和計算はASEのFIRE, FrechetCellFilterを使用し、力の収束基準を0.005 eV/Åに設定しています。熱物性（エントロピー、ヘルムホルツ自由エネルギー、熱容量）は、 Phonopy を使用して、有限変位を伴うスーパーセル法により求めた力定数から計算します。力定数計算時の変位距離が結果に与える影響を調べるため、変位距離(distance)を[0.01Å, 0.03Å, 0.05Å, 0.07Å, 0.1Å, 0.2Å]と変えて評価します。

#### 変位距離と熱物性の関係 ​

下図は、300Kでのエントロピー、ヘルムホルツ自由エネルギー、熱容量に関して、変位距離を変化させた際のPFPとDFT計算 (distance=0.01Å) との絶対平均誤差(MAE)を示したものです。なお、計算時間短縮のため、データセットからランダムに500点を抽出して評価を行いました。変位距離が大きい領域(distance>0.1Å)では非調和項の影響が無視できなくなり、誤差が増加することが確認されました。変位距離を小さくすることで調和近似の精度が向上する傾向が見られる一方で、変位距離が小さい領域(distance<0.03Å)において逆に誤差が増大する現象が観測されました。これは、局所安定点付近のポテンシャルエネルギー曲面の滑らかさに起因する問題であり、微小変位に対する力の計算精度に課題があることを示唆しています。このことから、PFP v8.0.0においては最適な変位距離は0.05~0.1Åの範囲であり、この条件下でフォノン物性の最も良好な再現性が得られることが期待できます。微小変位領域における精度改善については今後の課題として継続的に取り組んでいます。

#### 熱物性の再現性 ​

フォノンの平均振動数、フォノンの最大振動数、および熱物性に関して、変位距離0.05ÅにおけるPFP v8.0.0とDFTとのMAEを以下に示します。比較のためにPFP v7.0.0の結果も示しています。

Model| Δ Average Frequency [K]| Δ Max Frequency [K]| Δ Entropy [J/K/mol]| Δ Free Energy [kJ/mol]| Δ Heat Capacity [J/K/mol]  
---|---|---|---|---|---  
v8.0.0| 4.33| 24.50| 9.46| 3.07| 2.49  
v7.0.0| 5.25| 16.42| 10.60| 3.34| 2.90  
v8.0.0 (w/o mp-23907)| 4.27| 12.83| 9.46| 3.02| 2.49  
v7.0.0 (w/o mp-23907)| 5.16| 16.33| 10.60| 3.34| 2.90  
  
(注) v7.0.0では構造緩和が収束しなかった構造([mp-722890](https://next-gen.materialsproject.org/materials/mp-722890))および熱物性の計算が発散した構造([mp-37906](https://legacy.materialsproject.org/materials/mp-37906/), [mp-961706](https://next-gen.materialsproject.org/materials/mp-961706/))を評価から除外しています。また、v7.0.0 (w/o mp-23907) および v8.0.0 (w/o mp-23907)は、水素単結晶([mp-23907](https://next-gen.materialsproject.org/materials/mp-23907))を外れ値として評価から除外した結果です。

v8.0.0では、v7.0.0と比べて熱物性の精度が全体的に向上している一方で、フォノンの最大振動数については精度が悪化しています。この原因は、水素の単結晶([mp-23907](https://next-gen.materialsproject.org/materials/mp-23907))のフォノン最大振動数がDFTの結果と大きく乖離しているためでした。具体的には、PFP v8.0.0では117,971 K、DFT計算では1,724 Kと約68倍の差が生じていました。この構造を外れ値として評価から除外した場合、表および下図に示すように、フォノンの最大振動数を含む各種指標が改善し、PFP v8.0.0の予測精度がv7.0.0に比べ向上していることが確認できます。

PFP v8.0.0| PFP v7.0.0  
---|---  
|   
  
水素単結晶のフォノンの最大振動数の外れ値は、他の構造・熱物性と比較し実用上大きな問題とはなりにくい一方で、次バージョン以降で取り組むべき課題の一つと考えています。

## Appendix ​

### 有機液体密度のベンチマークの詳細な値 ​

PFP v8.0.0での有機液体密度のベンチマークにおいて、計算された密度と実験値との比較を以下に示します。実験値はPubChem[21:1]およびSigma-Aldrich[22:1]から取得しました。

物質名| R2SCAN_PLUS_D3モード(g/cm3)| PBE_PLUS_D3モード(g/cm3)| 実験値(g/cm3)[21:2] [22:2]  
---|---|---|---  
ethylene glycol| 1.166409| 1.090114| 1.115  
acetic acid| 1.111795| 0.995436| 1.051  
benzene| 0.831481| 0.832955| 0.879  
n-octane| 0.728450| 0.685701| 0.703  
dimethyl sulfoxide| 1.098966| 1.072171| 1.10  
water| 1.131187| 0.993455| 0.9950  
pyridine| 0.935546| 0.932742| 0.98272  
toluene| 0.851612| 0.876004| 0.867  
diethyl ether| 0.769595| 0.704671| 0.714  
isopropyl alcohol| 0.842338| 0.807849| 0.785  
chloroform| 1.415056| 1.342654| 1.4832  
N,N-dimethylformamide| 0.960495| 0.910953| 0.9445  
acetonitrile| 0.810502| 0.777934| 0.787  
methanesulfonic acid| 1.564527| 1.443646| 1.4812  
trifluoroethanol| 1.445350| 1.201619| 1.387  
1-methylnaphthalene| 0.974886| 1.019425| 1.0202  
γ-Butyrolactone| 1.143243| 1.079012| 1.1286  
xylene| 0.858634| 0.843388| 0.861  
propylene carbonate| 1.179535| 1.099879| 1.2047  
tetrahydrofuran| 0.913074| 0.878096| 0.888  
n-hexane| 0.732410| 0.677622| 0.659  
ethyl acetate| 0.948614| 0.848718| 0.9003  
perfluorooctane| 1.705958| 1.404300| 1.766  
hexamethylphosphoric triamide| 1.075276| 1.022574| 1.024  
dioxane| 1.045178| 0.964017| 1.0337  
  
## 参考文献 ​

* * *

  1. "matminer" <https://hackingmaterials.lbl.gov/matminer/dataset_summary.html#expt-formation-enthalpy-kingsbury> ↩︎

  2. O. Kubaschewski, C. Alcock, and P. Spencer, Materials Thermochemistry, 6th ed. (Pergamon Press, 1993). ↩︎

  3. M. W. Chase, Nist-janaf thermochemical tables (1998). ↩︎

  4. G. Kim, S. V. Meschel, P. Nash, and W. Chen, Scientific Data 4, 10.1038/sdata.2017.162 (2017). ↩︎

  5. G. Kim, S. Meschel, P. Nash, and W. Chen, Experimental formation enthalpies for intermetallic phases and other inorganic compounds (2017). ↩︎

  6. “Materials Project” <https://materialsproject.org/> ↩︎ ↩︎ ↩︎ ↩︎ ↩︎ ↩︎

  7. "Materials Project: Anion and GGA/GGA+U Mixing" <https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing> ↩︎

  8. Brew, R. R.; Nelson, I. A.; Binayeva, M.; Nayak, A. S.; Simmons, W. J.; Gair, J. J.; Wagen, C. C. Wiggle150: Benchmarking Density Functionals and Neural Network Potentials on Highly Strained Conformers. J. Chem. Theory Comput. 2025, 21 (8), 3922–3929. <https://doi.org/10.1021/acs.jctc.5c00015>. ↩︎

  9. Anstine, D. M.; Zubatyuk, R.; Isayev, O. AIMNet2: A Neural Network Potential to Meet Your Neutral, Charged, Organic, and Elemental-Organic Needs. Chem. Sci. 2025, 16 (23), 10228–10244. <https://doi.org/10.1039/D4SC08572H>. ↩︎

  10. "AceFF Machine Learning Potentials" <https://huggingface.co/collections/Acellera/aceff-machine-learning-potentials-6810a5f0ab054efab65cab08> ↩︎

  11. Devereux, C.; Smith, J. S.; Huddleston, K. K.; Barros, K.; Zubatyuk, R.; Isayev, O.; Roitberg, A. E. Extending the Applicability of the ANI Deep Learning Molecular Potential to Sulfur and Halogens. J. Chem. Theory Comput. 2020, 16 (7), 4192–4202. <https://doi.org/10.1021/acs.jctc.0c00121>. ↩︎

  12. Neumann, M.; Gin, J.; Rhodes, B.; Bennett, S.; Li, Z.; Choubisa, H.; Hussey, A.; Godwin, J. Orb: A Fast, Scalable Neural Network Potential. arXiv October 29, 2024. <https://doi.org/10.48550/arXiv.2410.22570>. ↩︎

  13. Batatia, I.; Kovács, D. P.; Simm, G. N. C.; Ortner, C.; Csányi, G. MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields. arXiv January 26, 2023. <https://doi.org/10.48550/arXiv.2206.07697>. ↩︎

  14. "The GMTKN55 database" <https://goerigk.chemistry.unimelb.edu.au/research/the-gmtkn55-database/> ↩︎ ↩︎ ↩︎ ↩︎ ↩︎ ↩︎ ↩︎

  15. W. Chen, F. Ambrosio, G. Miceli, and A. Pasquarello, Ab initio Electronic Structure of Liquid Water (2016), <https://doi.org/10.1103/PhysRevLett.117.186401> ↩︎

  16. F. Belleflamme, and Jürg Hutter, Radicals in aqueous solution: assessment of density-corrected SCAN functional (2023), <https://doi.org/10.1039/D3CP02517A> ↩︎

  17. A. K. Soper, The Radial Distribution Functions of Water as Derived from Radiation Total Scattering Experiments: Is There Anything We Can Say for Sure? (2013), <https://doi.org/10.1155/2013/279463> ↩︎

  18. J. Wiktor, F. Ambrosio, A. Pasquarello, Note: Assessment of the SCAN+rVV10 functional for the structure of liquid water (2017), <https://doi.org/10.1063/1.5006146> ↩︎ ↩︎ ↩︎ ↩︎

  19. “Crystallography Open Database” <https://www.crystallography.net/cod/> ↩︎ ↩︎ ↩︎

  20. So Takamoto, Chikashi Shinagawa, Daisuke Motoki, Kosuke Nakago, Wenwen Li, Iori Kurata, Taku Watanabe, Yoshihiro Yayama, Hiroki Iriguchi, Yusuke Asano, Tasuku Onodera, Takafumi Ishii, Takao Kudo, Hideki Ono, Ryohto Sawada, Ryuichiro Ishitani, Marc Ong, Taiki Yamaguchi, Toshiki Kataoka, Akihide Hayashi, Nontawat Charoenphakdee, and Takeshi Ibuka, “Towards universal neural network potential for material discovery applicable to arbitrary combination of 45 elements,” Nature Communications 13, 2991 (2022). <https://doi.org/10.1038/s41467-022-30687-9> ↩︎ ↩︎ ↩︎

  21. "PubChem" <https://pubchem.ncbi.nlm.nih.gov/> ↩︎ ↩︎ ↩︎

  22. "Sigma-Aldrich" <https://www.sigmaaldrich.com/> ↩︎ ↩︎ ↩︎

  23. A. Loew, D. Sun, HC. Wang, et al. Universal machine learning interatomic potentials are ready for phonons. npj Comput Mater 11, 178 (2025). <https://doi.org/10.1038/s41524-025-01650-1> ↩︎ ↩︎

  24. "MDR phonon calculation database" <https://mdr.nims.go.jp/collections/d7aab932-8512-4b9a-b93d-b61f6e5e7019> ↩︎



