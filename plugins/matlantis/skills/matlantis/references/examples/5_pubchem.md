# PubChemの利用例
## アクセス方法
PubChemのトップページにアクセスすると、検索ウインドウがあるページが表示されます。  
https://pubchem.ncbi.nlm.nih.gov/  
Webページからの検索、構造のダウンロード等に登録は必要ありません。

&lt;img src="./assets/tutorial_pubchem/toppage.png" alt="toppage" /&gt;

## 検索方法
ここでは、トルエン(toluene)を検索したいケースを考えます。  

### 化合物名で検索
PubChemでは、化合物名や化学式などで検索することができます。さっそく、`toluene`と検索ウインドウに入力して、検索してみましょう。  
https://pubchem.ncbi.nlm.nih.gov/#query=toluene

検索結果一覧が表示されます。今回（2021/03時点）は、検索ウインドウのすぐ下に`COMPOUND BEST MATCH`として、目的とするトルエンが表示され、その下に検索結果一覧が表示されました。  
検索対象によっては、大量の検索結果が表示され、目的とする化合物を見つけにくいことがあります。その場合は、検索結果一覧の上にある`Filters`をクリックし、分子量などで絞り込むことができます。

&lt;img src="assets/tutorial_pubchem/search_by_name.gif" alt="search_by_name" /&gt;

### 化学式で検索
別の方法で検索することもできます。たとえば、`C6H5CH3`のように化学式で検索することができます。  
https://pubchem.ncbi.nlm.nih.gov/#query=C6H5CH3

&lt;img src="assets/tutorial_pubchem/search_by_formula.gif" alt="search_by_formula" /&gt;

### 構造式で検索
また、トップページ( https://pubchem.ncbi.nlm.nih.gov/ )からであれば、検索ウインドウ下の`Draw Structure`をクリックすると、構造式を書くことができます。  
構造式を書いた後、`Search for This Structure`をクリックすると、検索結果が表示され、検索ウインドウには`C1=CC=CC=C1C`と入力されています。
これは、SMILES表記と呼ばれ、文字列で化学構造を表現したものです。  
https://pubchem.ncbi.nlm.nih.gov/#query=C1%3DCC%3DCC%3DC1C

&lt;img src="assets/tutorial_pubchem/search_by_draw_structure.gif" alt="search_by_structure" /&gt;

さて、目的とするトルエンのページを開いて、次に進みましょう。  
https://pubchem.ncbi.nlm.nih.gov/compound/1140

## 構造のダウンロード

化合物のページの最初には、分子式(Molecular Formula)や分子量(Molecular Weight)などとともに、PubChem CIDが表示されています。PubChem CIDは、各化合物に紐付けられているユニークな番号です。  
少し下に進むと、"1 Structures"があり、ここから構造ファイルをダウンロードすることができます。  
https://pubchem.ncbi.nlm.nih.gov/compound/1140#section=Structures

今回は、2次元構造("1.1 2D Structure")に加えて、3次元構造("1.2 3D Conformer")があります。  
3次元構造を初期構造とした場合の方が、分子のoptなどにかかる時間を削減することができます。
ただし、3次元構造はPubChemによって生成されたものであり、妥当でない構造の場合もあるので注意してください。

構造が表示されている枠の外側、右上をの`Download`をクリックすると、複数の選択肢があらわれます。
`SDF`の横の`SAVE`をクリックすると、拡張子sdfのファイルがダウンロードされます。 

ファイルをダウンロードして、次に進みましょう。
今回のファイル名は、`Conformer3D_CID_1140.sdf`でした。

&lt;img src="assets/tutorial_pubchem/download.gif" alt="download" /&gt;

## Material Notebookでの利用例
ダウンロードしたファイルをNotebookにアップロードすることができます。
Material Notebookのページ左上の`Upload files`から、あるいは単純にファイルを左のフォルダ欄にドラッグ&amp;ドロップすることでファイルがアップロードできます。

アップロードしたSDFファイルは、pythonのASEライブラリで扱うことができます。

&lt;img src="assets/tutorial_pubchem/upload.gif" alt="upload" /&gt;


```python
# If you have installed 'pfp-api-client', you can skip.
!pip install pfp-api-client
```


```python
import ase.io
from ase.visualize import view

atoms = ase.io.read("Conformer3D_CID_1140.sdf")

v = view(atoms, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


ここから先は、PFPと組合せて任意の計算が出来ます。

ここでは、構造最適化を行ってみましょう。


```python
from ase.io import Trajectory
from ase.optimize import BFGS

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
calculator = ASECalculator(estimator)
atoms.calc = calculator

opt = BFGS(atoms, trajectory="pubchem_sample_opt.traj")
opt.run()

traj = Trajectory("pubchem_sample_opt.traj")
v = view(traj, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```
