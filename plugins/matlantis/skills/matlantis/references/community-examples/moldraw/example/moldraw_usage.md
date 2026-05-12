Copyright ENEOS, Corp., Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project

# moldraw

JSMEを用いて分子を描画しASE　atomsを返すUIツールです。  #eneos ibuka 2022/12/8

v0.1.0 2023/10/25 　Pandasバージョンアップ対応 パッケージ化対応,　引数からcalculator 設定可能に変更　3.8,3.9で動作確認

## セットアップ
 Install.ipynb を開き実行してください。一度実行すれば他のフォルダでも実行できます。
 python3.8,3.9で動作確認しています。利用バージョンそれぞれでのinstallが必要です。


```python
import moldraw
```


    



```python
m = moldraw.Moldraw()
```


    VBox(children=(HBox(children=(VBox(children=(VBox(children=(Output(layout=Layout(height='290px', width='440px'…


TIPS:　 ▼▲の一番下のメニューから　SMILES→図構造ができる。左回転ボタンでUNDOができる。

登録した分子はdf経由でアクセス可能。


```python
display(m.df)
# atomsへアクセス
m.df.iloc[0].atoms
```


&lt;div&gt;
&lt;style scoped&gt;
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
&lt;/style&gt;
&lt;table border="1" class="dataframe"&gt;
  &lt;thead&gt;
    &lt;tr style="text-align: right;"&gt;
      &lt;th&gt;&lt;/th&gt;
      &lt;th&gt;SMILES&lt;/th&gt;
      &lt;th&gt;formula&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;maxforce&lt;/th&gt;
      &lt;th&gt;opt&lt;/th&gt;
      &lt;th&gt;conf_no&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;c1ccccc1&lt;/td&gt;
      &lt;td&gt;C6H6&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-0.2995, -1.35...&lt;/td&gt;
      &lt;td&gt;-60.97546&lt;/td&gt;
      &lt;td&gt;1.115038&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;NaN&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;c1ccccc1&lt;/td&gt;
      &lt;td&gt;C6H6&lt;/td&gt;
      &lt;td&gt;(Atom('C', [0.7943, -1.149...&lt;/td&gt;
      &lt;td&gt;-60.96755&lt;/td&gt;
      &lt;td&gt;1.540152&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;0.0&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;





    Atoms(symbols='C6H6', pbc=False)



イテレータでfor loopを用いてもアクセス可能


```python
for i, atoms in enumerate(m):
    print(i, atoms)
```

    0 Atoms(symbols='C6H6', pbc=False, calculator=ASECalculator(...))
    1 Atoms(symbols='C6H6', pbc=False)
    2 Atoms(symbols='C6H6', pbc=False)


## SMILESからのConformaer search

`get_conformersdf` 関数を使うことで、SMILES指定をしての直接コンフォーマ探索も可能


```python
from moldraw.moldraw import get_conformersdf

confo = get_conformersdf("CCCCC", top=3, opt=True)
confo 
```

    6 conformers found (n_confs:100 n_atoms:17,rmsd_threshold:0.22)
    start opt top 3 conformers
     CCCCC 1 / 3 opt 36 steps done
     CCCCC 2 / 3 opt 46 steps done
     CCCCC 3 / 3 opt 32 steps done





&lt;div&gt;
&lt;style scoped&gt;
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
&lt;/style&gt;
&lt;table border="1" class="dataframe"&gt;
  &lt;thead&gt;
    &lt;tr style="text-align: right;"&gt;
      &lt;th&gt;&lt;/th&gt;
      &lt;th&gt;SMILES&lt;/th&gt;
      &lt;th&gt;formula&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;maxforce&lt;/th&gt;
      &lt;th&gt;opt&lt;/th&gt;
      &lt;th&gt;conf_no&lt;/th&gt;
      &lt;th&gt;conf_no_opt&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.573962872795...&lt;/td&gt;
      &lt;td&gt;-69.826895&lt;/td&gt;
      &lt;td&gt;0.040378&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;0&lt;/td&gt;
      &lt;td&gt;0&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.39375624465...&lt;/td&gt;
      &lt;td&gt;-69.786662&lt;/td&gt;
      &lt;td&gt;0.049567&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;1&lt;/td&gt;
      &lt;td&gt;1&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.372965376100...&lt;/td&gt;
      &lt;td&gt;-69.779179&lt;/td&gt;
      &lt;td&gt;0.045156&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;2&lt;/td&gt;
      &lt;td&gt;2&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [1.7747, -0.788...&lt;/td&gt;
      &lt;td&gt;-69.373104&lt;/td&gt;
      &lt;td&gt;2.622761&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;3&lt;/td&gt;
      &lt;td&gt;3&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-1.8406, -0.26...&lt;/td&gt;
      &lt;td&gt;-69.348328&lt;/td&gt;
      &lt;td&gt;2.297797&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;4&lt;/td&gt;
      &lt;td&gt;4&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;5&lt;/th&gt;
      &lt;td&gt;CCCCC&lt;/td&gt;
      &lt;td&gt;C5H12&lt;/td&gt;
      &lt;td&gt;(Atom('C', [1.5662, 0.0888...&lt;/td&gt;
      &lt;td&gt;-69.146562&lt;/td&gt;
      &lt;td&gt;2.633826&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;5&lt;/td&gt;
      &lt;td&gt;5&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;



コンフォーマーは直接関数でも実行可能。大分子では遅い。


```python
%%time
confo2 = get_conformersdf("CCCCCCC", opt=False)
confo2
```

    25 conformers found (n_confs:100 n_atoms:23,rmsd_threshold:0.28)
    CPU times: user 1.99 s, sys: 536 ms, total: 2.53 s
    Wall time: 4.52 s





&lt;div&gt;
&lt;style scoped&gt;
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
&lt;/style&gt;
&lt;table border="1" class="dataframe"&gt;
  &lt;thead&gt;
    &lt;tr style="text-align: right;"&gt;
      &lt;th&gt;&lt;/th&gt;
      &lt;th&gt;SMILES&lt;/th&gt;
      &lt;th&gt;formula&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;maxforce&lt;/th&gt;
      &lt;th&gt;opt&lt;/th&gt;
      &lt;th&gt;conf_no&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [3.1262, 0.1417...&lt;/td&gt;
      &lt;td&gt;-95.186902&lt;/td&gt;
      &lt;td&gt;2.073040&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;0&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.8209, -0.42...&lt;/td&gt;
      &lt;td&gt;-95.166704&lt;/td&gt;
      &lt;td&gt;1.845676&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;1&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.781, 1.2439...&lt;/td&gt;
      &lt;td&gt;-95.147992&lt;/td&gt;
      &lt;td&gt;1.424633&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;2&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [3.569, 0.1564,...&lt;/td&gt;
      &lt;td&gt;-95.142429&lt;/td&gt;
      &lt;td&gt;1.780890&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;3&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.6159, -0.15...&lt;/td&gt;
      &lt;td&gt;-95.130356&lt;/td&gt;
      &lt;td&gt;1.795929&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;4&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;5&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.731, 1.0508,...&lt;/td&gt;
      &lt;td&gt;-95.124852&lt;/td&gt;
      &lt;td&gt;1.902729&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;5&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;6&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.8554, 0.007...&lt;/td&gt;
      &lt;td&gt;-95.105848&lt;/td&gt;
      &lt;td&gt;1.855395&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;6&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;7&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-3.5235, -0.42...&lt;/td&gt;
      &lt;td&gt;-95.101927&lt;/td&gt;
      &lt;td&gt;1.861386&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;7&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;8&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.8141, 1.0146...&lt;/td&gt;
      &lt;td&gt;-95.073047&lt;/td&gt;
      &lt;td&gt;1.621446&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;8&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;9&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.9278, 0.448...&lt;/td&gt;
      &lt;td&gt;-95.069335&lt;/td&gt;
      &lt;td&gt;2.097057&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;9&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;10&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [3.7123, -0.064...&lt;/td&gt;
      &lt;td&gt;-95.067109&lt;/td&gt;
      &lt;td&gt;1.901971&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;10&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;11&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.1378, 1.0463...&lt;/td&gt;
      &lt;td&gt;-95.067095&lt;/td&gt;
      &lt;td&gt;1.967215&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;11&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;12&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [3.4119, 0.1057...&lt;/td&gt;
      &lt;td&gt;-95.017548&lt;/td&gt;
      &lt;td&gt;1.423544&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;12&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;13&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [3.4199, -0.111...&lt;/td&gt;
      &lt;td&gt;-95.011299&lt;/td&gt;
      &lt;td&gt;1.934434&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;13&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;14&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.7396, -0.22...&lt;/td&gt;
      &lt;td&gt;-94.984482&lt;/td&gt;
      &lt;td&gt;1.808125&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;14&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;15&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.9296, -1.134...&lt;/td&gt;
      &lt;td&gt;-94.944889&lt;/td&gt;
      &lt;td&gt;2.153117&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;15&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;16&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.6008, -0.40...&lt;/td&gt;
      &lt;td&gt;-94.943989&lt;/td&gt;
      &lt;td&gt;1.951802&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;16&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;17&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.6142, 0.667...&lt;/td&gt;
      &lt;td&gt;-94.932814&lt;/td&gt;
      &lt;td&gt;2.050968&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;17&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;18&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.7207, -0.19...&lt;/td&gt;
      &lt;td&gt;-94.925487&lt;/td&gt;
      &lt;td&gt;1.604887&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;18&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;19&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.1879, 1.309...&lt;/td&gt;
      &lt;td&gt;-94.910040&lt;/td&gt;
      &lt;td&gt;2.521500&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;19&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;20&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-1.9826, 1.417...&lt;/td&gt;
      &lt;td&gt;-94.754874&lt;/td&gt;
      &lt;td&gt;3.535223&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;20&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;21&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.4726, 0.6282...&lt;/td&gt;
      &lt;td&gt;-94.669979&lt;/td&gt;
      &lt;td&gt;3.534907&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;21&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;22&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-1.729, -0.812...&lt;/td&gt;
      &lt;td&gt;-94.591812&lt;/td&gt;
      &lt;td&gt;3.556740&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;22&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;23&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-2.3894, 0.250...&lt;/td&gt;
      &lt;td&gt;-94.576384&lt;/td&gt;
      &lt;td&gt;4.205362&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;23&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;24&lt;/th&gt;
      &lt;td&gt;CCCCCCC&lt;/td&gt;
      &lt;td&gt;C7H16&lt;/td&gt;
      &lt;td&gt;(Atom('C', [2.0665, -0.021...&lt;/td&gt;
      &lt;td&gt;-94.553474&lt;/td&gt;
      &lt;td&gt;4.005065&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;24&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;



## DataFrameに対するグラフ描画ツール

DataFrameに対して汎用的に使えるグラフ用のUIを提供するクラス`Dfplot`。

引数にDataFrameを入れることで、x軸・y軸を指定しての描画が可能です。


```python
from moldraw import Dfplot
p = Dfplot(confo2)
```


    VBox(children=(HBox(children=(Select(description='col1:', layout=Layout(border='0px', height='27px', width='20…

