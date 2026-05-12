Copyright Preferred Computational Chemistry, Inc as contributors to Matlantis contrib project

# Structure search with CrySPY &amp; PFP

結晶構造探索を行うツールであるCrySPYをMatlantis上のPFP力場を用いて動かす方法。

## 事前セットアップ

以下のライブラリを用意します。

**CrySPY**
 - https://github.com/Tomoki-YAMASHITA/CrySPY

執筆時点の最新版はv1.1.1: https://github.com/Tomoki-YAMASHITA/CrySPY/releases/tag/v1.1.1

インストールに関する注意事項
- CrySPY v1.0.0以降のversionではPython 3.8以上が要求されるためご注意ください。
- Matlantis Featuresが既にインストールされている場合、双方のライブラリで要求されるNumPy, SciPyのバージョン不整合に起因したErrorが表示される可能性があります。この場合、CrySPYが要求するバージョンを再インストールした上でCrySPYのインストールが完了します。

**LAMMPS: matlantis-lammpsを利用**

 - Experimental packageからコピーし、ipynbに従ってインストール。


```python
# Download CrySPY v1.1.1
!pip install csp-cryspy
```

 - pyxtal
   - docs: https://pyxtal.readthedocs.io/en/latest/
   - github: https://github.com/qzhu2017/PyXtal

CrySPYが内部で利用をしている、空間群の対称性を考慮した操作を行うためのライブラリ。

以下のコマンドでインストールできます。


```python
!pip install pyxtal
```

## 公式のLAMMPS Exampleの実行

CrySPYの実行は、python programを実行していく形となり、Terminal上で行うことを推奨します。

以下のようなコードをTerminal上で実行することで、CrySPY packageに事前に含まれているLAMMPS連携Exampleを実行することができます。

```
# python=3.8のActivate
source ~/.py38/bin/activate

# LAMMPS Exampleのディレクトリへ移動 (ディレクトリ位置はユーザーごとに異なるので注意)
cd CrySPY/example/LAMMPS_bash_Si16_lj_RS

# CrySPYの実行
cryspy
```

※ 2回目以降の実行で、初期化してから実行したい場合は、出力として作成される data, work ディレクトリおよび cryspy.out, cryspy.stat, lock_cryspy のファイルのうち存在するものを削除してから実行してください。

上記のように、実行したい設定ファイル郡が格納されているディレクトリ以下に `cd` コマンドで移動後、その場所から`cryspy` を実行することでツールを動かすことができます。

`cryspy` はJobの進捗状況に応じて、

 - 1. 構造候補生成
 - 2.1 次の構造候補の計算に移る
 - 2.2 今行っている構造の計算が終了したかを確認
 - 2. すべての計算が終了したか判定

といったCheckを繰り返すため、 **定期的に `cryspy` を呼び出す必要があります。**



例えば、上記のLAMMPS exampleの場合、正常に動作していると以下のようなログが出力されます。

**1回目の`cryspy`:**

```
2023/07/04 09:20:42
CrySPY 1.1.1
Start cryspy.py
Number of MPI processes: 1

Read input file, cryspy.in
Save input data in cryspy.stat

# --------- Generate initial structures
# ------ mindist
Si - Si 1.11
Structure ID      0 was generated. Space group: 203 --&gt; 227 Fd-3m
Structure ID      1 was generated. Space group: 176 --&gt; 176 P6_3/m
Structure ID      2 was generated. Space group:  59 --&gt;  59 Pmmn
Structure ID      3 was generated. Space group:   4 --&gt;   4 P2_1
Structure ID      4 was generated. Space group: 159 --&gt; 159 P31c

Elapsed time for structure generation: 0:00:00.560896
```

→候補構造を5つ生成しています。実際に生成される構造はRandomとなるため実行毎に異なります。


**2回目の`cryspy`:**

```
2023/07/04 09:21:09
CrySPY 1.1.1
Restart cryspy.py
Number of MPI processes: 1



# ---------- job status
ID      0: submit job, Stage 1
ID      1: submit job, Stage 1
```

ID=0, 1 の構造についてのみ、実際にエネルギー計算(構造最適化含む)が開始されます。

**3回目の`cryspy`:**

```
# ---------- job status
ID      0: Stage 1 Done!
    submitted job, ID      0 Stage 2
ID      1: Stage 1 Done!
    submitted job, ID      1 Stage 2
```

2回目の実行から、ある程度の時間を空けてから(構造最適化が終了したあとに)実行を行うと、

このような形で、実際にすべての計算が終了するまで、何度も`cryspy` を実行し確認する必要があります。

### watch コマンドによる実行


この操作は watch コマンドを利用することで、自動化することができます。

以下のコマンドは20秒ごとに、 `cryspy` を実行し続ける例です。

```
watch -n 20 cryspy
```


これをTerminalで実行し、待ち続けていると、CrySPYの構造探索が進み、最終的に以下のような出力が得られます。

```
2023/07/04 09:28:03
CrySPY 1.1.1
Restart cryspy.py
Number of MPI processes: 1



# ---------- job status
Done all structures!
```

## PFPを用いたExampleの実行


上記のLAMMPS Exampleでは、Si 結晶の構造探索をLJ potentialで行っていました。&lt;br/&gt;
ここでは、力場をPFPに変更して構造探索を行ってみましょう。

`LAMMPS_bash_Si16_lj_RS_PFP` ディレクトリにPFPを使うように設定ファイルを書き換えたファイルが存在します。
具体的には`calc_in/in.si_1` および `calc_in/in.si_2` がLAMMPSの入力設定ファイルになりますが、これらを以下のようにPFP力場を使うように書き換えています。

```
pair_style pfp_api v4.0.0
pair_coeff * * species Si
```

Terminalを起動し、 `LAMMPS_bash_Si16_lj_RS_PFP` ディレクトリに移動後、以下のコマンドを実行してみてください。

```
# python=3.8のActivate
source ~/.py38/bin/activate

# PFP LAMMPS Exampleのディレクトリへ移動
cd LAMMPS_bash_Si16_lj_RS_PFP

# CrySPYの実行 (watchコマンドを使って定期実行するVersion)
watch -n 20 cryspy
```

(実行には、`matlantis-lammps` のインストールが必要です。)

※ 2回目以降の実行で、初期化してから実行したい場合は、出力として作成される data, work ディレクトリおよび cryspy.out, cryspy.stat, lock_cryspy のファイルのうち存在するものを削除してから実行してください。

PFPを用いたExampleは上記のLJ potentialを使用した場合と比べ、構造最適化などに時間がかかるため、１つ１つのJobが終わるまでに時間がかかり、"ID      0: still queueing or running" のようなログがある程度長く表示され続けます (数分から10分程度)。
時間を待つと、全てのJobが終了します。


### PFPが正常に実行されているか確認する方法

まず注意点として、 **`cryspy` の実行は、Terminal上で行う必要があります。**&lt;br/&gt;
Jupyter LabのNotebook上のセルから `!cryspy` と実行した場合には、"ID      0: still queueing or running" が表示されるが、PFPがBackgroundで実行されずにスタックしてしまうことが確認されています。

Terminal上で、"ID      1: still queueing or running" が表示され続けていて、正常に実行が進んでいるかわからないときは、`work/000001/log.struc`, `work/000001/out.si`, `work/000001/stage1_out.si` などを開き計算ログに進捗があるかを確認してください。("000001" は実行Jobの番号に応じて変更してください。)
