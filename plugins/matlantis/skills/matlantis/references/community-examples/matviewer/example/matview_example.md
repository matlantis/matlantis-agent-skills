Copyright ENEOS, Corp. and  Preferred Computational Chemistry as contributors to Matlantis contrib project  
:2023/10/30 ENEOS Ibuka  
### Matviewer  
 Matlantis で分子、結晶構造の取得、加工、OPT,NEB,VIB などの解析をGUIで行うツールです。  
 インストールをInstall.ipynbで一度実行すれば他のフォルダでも実行出来ます。  Python3.8,3.9に対応しています

機能説明  

・全体  
　　・初期左側　[Setting] が表示されCalc,面サイズ,色を設定できます。全画面モードも対応しています  
　　・MatviewerではTraj番号#　とTraj内のimage番号##で表示、編集、計算を管理しています  
　　・複数の構造を管理するためにはAdd# Add##　ボタンで更新領域Noを確保します。  
　　　Del#,Del##で表示されている領域を削除できます。  
　　・複数の##を持つTraj#ではplot ボタンを押すことでepot ,force を表示できます。  
　　・[Save]  現在の登録すべてが保存されます。　次回起動するときはこれが読み込まれます。(zip化も同時にされ過去復元も可能）  
　　・下のチェックボックスcalcで現在のenergy と max forceが表示されます。  
　　　setボタンで登録した値との差分dEも表示されるのでNEBや吸着で便利です  
  
・ADD MENU  
　　[FILE] 既存のAtoms または　traj ファイルを読込みます。複数構造に対応しています  
　　[SILES or Draw]    JSMEによる分子描画またはSMILESからAtoms作成できます。  
　　　　Conformer ボタンで同一SMILESでの他の構造を探索できます。  
　　[Crystal]   Materials Project のデータ（mp.gz)から検索読込できます。serach後　view ボタンclick で表示読込できます。  
  
・EDIT MENU  
　　[Move] 指定したindexの原子を移動、回転できます。　また原子距離を指定した移動もできます。　undo 対応しています  
  　　　　index指定方法はpotistion index また conectでは指定番号とつながった分子の設定ができます（-1で最終原子指定可能）  
  　　　　mult で接続判定値を制御できます  
　　[Utiledit]  原子置換、削除が可能です。また番号のソートやwrap も可能です。   
　　[Repeat]    atoms.repeat の機能です。　pbc=False の分子の場合自動で設定されますが、allowanceでサイズを設定できます   
　　[Cell]    Cell の情報表示と設定が可能です  
　　[Surface]  pymatgenのSlabGeneratorの機能です。 x,yは小さいスラブができるのでその後Repeat([2,2,1])等行うとよいです  
　　[Attach]  Atoms0のa0 とAtoms1 のa1 の距離を指定して合体したAtoms を作成します。  
　　　　他の原子同士の干渉ができるだけ無いように自動回転・移動した構造が作成されます。  
　　　　NEBの初終構造作成や　モノマー→ダイマー→オリゴマー作成, スラブへの分子追加等で利用できます。  
　　　　OPTを1以上の数設定すると距離固定のoptをすることができます。単純にA＋B(cellはAを踏襲）することもできます  
　　[Liquid] 分子と数および密度を指定し液体構造を作成します。PFN社が開発したliquid_generator Matlantis contrib を利用しています &lt;br/&gt; 
　　[Traj] Traj番号# image番号##を指定したコピー追加やスワップができます。 

  ・CALC MENU  
　　[Setting]    Calcの設定　や画面サイズ（右側）、色（アプリ全体　NGL）設定できます。　全画面モードも対応しています。         
　　[OPT]    LGBFS での構造最適化を行います。cellopt も選択可能。結果は./output/optにも保存.constraint はFixAtomsのみに対応 &lt;br/&gt; 
　　[VIB]  　振動解析を行います。 view で指定した振動モードの動きも確認できます。温度も設定できます。  
　　　　　前回実施したatoms および./output/vib/内のtrajではVIB情報をREAD可能です。  
　　[Thermo]  VIB設定されているatoms でThermo、HARMONIC条件でのエネルギー計算が可能です。  
　　[MD_NVT]  NVTができます。数千イタレーション以内での実行のみ動作確認しています。  


```python
import matviewer 
```


```python
matv = matviewer.Matviewer()
```


```python
matv
```


```python
#表示しているatoms
matv.atoms
```


```python
#latoms_list　ADD# ##で登録したAtoms のリストのリスト
matv.latoms_list[-1][0]
```


```python
#別作業からのAtomsの追加
import ase
o2 = ase.Atoms("O2",[[0,0,0],[0,0,1]])
matv + o2 
#matv.add_atoms(o2) 　でも可能
```


```python
#別作業からのAtoms listの追加
h2 = ase.Atoms("H2",[[0,0,0],[0,0,0.75]])
matv.add_atoms([h2,o2])
#matv + [o2,h2] 出も可能
```


```python
#ADD-Crystal で検索後
display("results",matv.crystalsearch.cands.head())
display("mp-db all",matv.crystalsearch.df.head())
```


```python
#ADD-MOF で検索後
display("results",matv.mofsearch.cands.head())
display("mof-db all",matv.mofsearch.df.head())
```

更新履歴

2023/8/1  :v0.0.1 初回リリース  
2023/8/16 :v0.0.2  
・[Setting] Calcの設定　や画面サイズ,全画面モード,色対応、LightMode可視対応  
・Installer autodE およびcython もない人が入るように 追加  
・Liquidgenerator を別ファイルに（torch  import メモリが大きく必要なときのみ実行）   
2023/8/18 :v0.0.3  
・OPTにfixbond 追加、　画面各種微修正、nglview メモリ削除し、動作重くなる問題回避  
2023/8/24 :v0.0.4:  
・画面修正　trajあるもののみスライダー表示するようにし　そこにplotボタンも追加  
・OPT,NEB,グラフ表示機能追加　、特に速さを体感する上でも価値あり。またTrajのplotも並列計算で高速化   
・plotly 一回目遅い問題回避のため、treadバックグラウンド実行  
・AttachのPBC対応 触媒表面吸着でも利用できるように  
2023/9/12 :v0.1.0  
 ・OOM対応最大値設定 Atoms:20000原子:Traj数500:ファイル表示500,500MB Liq:1000原子,vib500原子  
 ・分子　SMILES conformer 対応  
 ・vib 保存・読込機能　（Atoms の atoms.info にdeltaやH等書き込み  
 ・thermo 対応 IDEAL と HARMO （CRYSTALは未）　グラフも表示  
 ・EDITにULITEDIT:wrap機能,EDIT CELL設定機能,Traj:Traj編集機能,Attach:A＋B（Aセル保持）を追加  
 ・複数のイメージを返す処理は　Ｔraj# を追加し、それ以外は　# ## を上書きするように設定   
2023/10/30 :v0.1.1  
  ・Python3.9に対応（Pandas pikcle 互換)  
  ・Crystal IDやFormulaでも検索可能とした  
  ・MOFS追加  
2024/11/25 : v0.1.2  
  ・ASE 3.23に対応 (FixSymmetry, NEB等)  
  ・GUIの背景色を変更  

既知の問題
  ・メモリ１G以上程度消費します。エラー処理は不十分。自己責任にてご利用ください
