Copyright ENEOS, Corp., Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project

# DFTデータベース（MaterialsProject等）から構造を探索するGUI

Matlantis の利用バージョン表示機能もついています #eneos ibuka 2022/10/1

`crystalsearch` moduleは内部でMaterials Projectのデータベースを使用しています。


Database is from  
A. Jain*, S.P. Ong*, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson (*=equal contributions)  
The Materials Project: A materials genome approach to accelerating materials innovation  
APL Materials, 2013, 1(1), 011002.  
[doi:10.1063/1.4812323](http://dx.doi.org/10.1063/1.4812323)  
[[bibtex]](https://materialsproject.org/static/docs/jain_ong2013.349ca3156250.bib)  
Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)  


```python
import  crystalsearch as cs
```


    



```python
m = cs.Crystalsearch()
```


&lt;style&gt; .widget-gridbox {background-color: #222266;padding: 3px} .widget-toggle-button {padding: 1px; font-size: 12px;border-radius: 6px;}&lt;/style&gt;



    VBox(children=(GridspecLayout(children=(ToggleButton(value=False, description='H', layout=Layout(grid_area='wi…


以下[Fe,O,C]を選択し、[hullview] ボタン実施後の動作です。 検索結果は `インスタンス名.cands` で取得できます。


```python
cands = m.cands
```


```python
cands
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
      &lt;th&gt;material_id&lt;/th&gt;
      &lt;th&gt;pretty_formula&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;e_above_hull&lt;/th&gt;
      &lt;th&gt;is_hubbard&lt;/th&gt;
      &lt;th&gt;elements&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;mp-13&lt;/td&gt;
      &lt;td&gt;Fe&lt;/td&gt;
      &lt;td&gt;-8.470009&lt;/td&gt;
      &lt;td&gt;0.000000&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Fe]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;mp-47&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-36.261055&lt;/td&gt;
      &lt;td&gt;0.161517&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-1.1144813925056951e-06, 1.4509361...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;mp-48&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-36.881187&lt;/td&gt;
      &lt;td&gt;0.006484&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [0.0, 0.0, 6.5137786865234375], ind...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;mp-66&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-18.180737&lt;/td&gt;
      &lt;td&gt;0.136413&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [1.263497233390808, 0.7294805049896...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;mp-136&lt;/td&gt;
      &lt;td&gt;Fe&lt;/td&gt;
      &lt;td&gt;-16.744431&lt;/td&gt;
      &lt;td&gt;0.097793&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Fe]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [-1.0390743909738376e-06, 1.423605...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;173&lt;/th&gt;
      &lt;td&gt;mvc-12204&lt;/td&gt;
      &lt;td&gt;FeO2&lt;/td&gt;
      &lt;td&gt;-218.275894&lt;/td&gt;
      &lt;td&gt;0.271338&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [6.014296531677246, 3.451520204544...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;174&lt;/th&gt;
      &lt;td&gt;mvc-12205&lt;/td&gt;
      &lt;td&gt;Fe9O13&lt;/td&gt;
      &lt;td&gt;-142.633408&lt;/td&gt;
      &lt;td&gt;0.263188&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [5.969217777252197, 1.786785602569...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;175&lt;/th&gt;
      &lt;td&gt;mvc-12905&lt;/td&gt;
      &lt;td&gt;FeO2&lt;/td&gt;
      &lt;td&gt;-73.856781&lt;/td&gt;
      &lt;td&gt;0.179825&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [-1.005743384361267, 1.72052621841...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;176&lt;/th&gt;
      &lt;td&gt;mvc-13234&lt;/td&gt;
      &lt;td&gt;FeO2&lt;/td&gt;
      &lt;td&gt;-218.290543&lt;/td&gt;
      &lt;td&gt;0.270931&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [4.6083173751831055, 0.87333792448...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;177&lt;/th&gt;
      &lt;td&gt;mvc-15135&lt;/td&gt;
      &lt;td&gt;FeO2&lt;/td&gt;
      &lt;td&gt;-217.965958&lt;/td&gt;
      &lt;td&gt;0.279947&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [4.428077697753906, 2.556548357009...&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;178 rows × 7 columns&lt;/p&gt;
&lt;/div&gt;



pfpとMaterials Projectのエネルギーを比較します


```python
pfpdf = m.get_pfp_df()
pfpdf.head()
```

    calc 178 cands by pfp v6.0.0
    [WARNING] 111 / 178 atoms contain Hubbard U correction, its U parameter might be different between Materials project &amp; PFP. Please check carefully.



      0%|          | 0/178 [00:00&lt;?, ?it/s]


    end



    
![png](output_7_3.png)
    





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
      &lt;th&gt;material_id&lt;/th&gt;
      &lt;th&gt;pretty_formula&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;e_above_hull&lt;/th&gt;
      &lt;th&gt;is_hubbard&lt;/th&gt;
      &lt;th&gt;elements&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
      &lt;th&gt;pfpene&lt;/th&gt;
      &lt;th&gt;pfpfmax&lt;/th&gt;
      &lt;th&gt;lenatoms&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;mp-13&lt;/td&gt;
      &lt;td&gt;Fe&lt;/td&gt;
      &lt;td&gt;-8.470009&lt;/td&gt;
      &lt;td&gt;0.000000&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Fe]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
      &lt;td&gt;-8.255982&lt;/td&gt;
      &lt;td&gt;0.00000&lt;/td&gt;
      &lt;td&gt;1&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;mp-47&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-36.261055&lt;/td&gt;
      &lt;td&gt;0.161517&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [-1.1144813925056951e-06, 1.4509361...&lt;/td&gt;
      &lt;td&gt;-36.337563&lt;/td&gt;
      &lt;td&gt;0.00091&lt;/td&gt;
      &lt;td&gt;4&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;mp-48&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-36.881187&lt;/td&gt;
      &lt;td&gt;0.006484&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [0.0, 0.0, 6.5137786865234375], ind...&lt;/td&gt;
      &lt;td&gt;-36.896410&lt;/td&gt;
      &lt;td&gt;0.00008&lt;/td&gt;
      &lt;td&gt;4&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;mp-66&lt;/td&gt;
      &lt;td&gt;C&lt;/td&gt;
      &lt;td&gt;-18.180737&lt;/td&gt;
      &lt;td&gt;0.136413&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[C]&lt;/td&gt;
      &lt;td&gt;(Atom('C', [1.263497233390808, 0.7294805049896...&lt;/td&gt;
      &lt;td&gt;-18.185245&lt;/td&gt;
      &lt;td&gt;0.00000&lt;/td&gt;
      &lt;td&gt;2&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;mp-136&lt;/td&gt;
      &lt;td&gt;Fe&lt;/td&gt;
      &lt;td&gt;-16.744431&lt;/td&gt;
      &lt;td&gt;0.097793&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Fe]&lt;/td&gt;
      &lt;td&gt;(Atom('Fe', [-1.0390743909738376e-06, 1.423605...&lt;/td&gt;
      &lt;td&gt;-16.305622&lt;/td&gt;
      &lt;td&gt;0.00001&lt;/td&gt;
      &lt;td&gt;2&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;



pfp の結果でphasediagram　を描きます


```python
pd = cs.get_phasediagram(pfpdf , energy = "pfpene")
pfppdfig = cs.get_phasediagram_fig(pd, show_unstable=1 )
print("pfp")
display (pfppdfig)
print("Original Materials Project")
#1回上記で表示し　2回目の表示でなぜか綺麗（線青）になる
display (m.pdfig) 
```

    pfp



    FigureWidget({
        'data': [{'a': [0, 0, None, 16.0, 0, None, 0, 1.0, None, 0, 1.0, None, 16.0,
                        1.0, None, 0, 0, None, 0, 16.0, None],
                  'b': [4.0, 0, None, 0, 0, None, 4.0, 0, None, 4.0, 0, None, 0, 0,
                        None, 4.0, 4.0, None, 4.0, 0, None],
                  'c': [8.0, 2.0, None, 24.0, 2.0, None, 8.0, 0, None, 0, 0, None,
                        24.0, 0, None, 8.0, 0, None, 8.0, 24.0, None],
                  'hoverinfo': 'none',
                  'line': {'color': 'black', 'width': 1.5},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': 'ecf10d41-1608-4ef8-a888-bf49fe6bfdbe'},
                 {'a': [0, 0, 1.0],
                  'b': [4.0, 4.0, 0],
                  'c': [0, 8.0, 0],
                  'fill': 'toself',
                  'fillcolor': '#2E91E5',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'C—CO&lt;sub&gt;2&lt;/sub&gt;—Fe',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': 'ce3e5b7c-7132-46dc-8a28-532a0906fbc8'},
                 {'a': [0, 1.0, 16.0],
                  'b': [4.0, 0, 0],
                  'c': [8.0, 0, 24.0],
                  'fill': 'toself',
                  'fillcolor': '#E15F99',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'CO&lt;sub&gt;2&lt;/sub&gt;—Fe—Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '02fa485a-54cc-4edc-86da-97c9c89c42ba'},
                 {'a': [0, 16.0, 0],
                  'b': [4.0, 0, 0],
                  'c': [8.0, 24.0, 2.0],
                  'fill': 'toself',
                  'fillcolor': '#1CA71C',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'CO&lt;sub&gt;2&lt;/sub&gt;—Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;—O&lt;sub&gt;2&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '1d0baa73-cf30-4087-b7bb-053ae2b9f916'},
                 {'a': [0, 0, 16.0, 1.0, 0],
                  'b': [4.0, 0, 0, 0, 4.0],
                  'c': [8.0, 2.0, 24.0, 0, 0],
                  'cliponaxis': False,
                  'hoverinfo': 'text',
                  'hoverlabel': {'font': {'size': 14}},
                  'hovertext': [CO&lt;sub&gt;2&lt;/sub&gt; (18) &lt;br&gt; -1.314 eV/atom (Stable),
                                O&lt;sub&gt;2&lt;/sub&gt; (105) &lt;br&gt; 0.0 eV/atom (Stable),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (33) &lt;br&gt; -0.335 eV/atom
                                (Stable), Fe (0) &lt;br&gt; 0.0 eV/atom (Stable), C (103)
                                &lt;br&gt; 0.0 eV/atom (Stable)],
                  'marker': {'color': 'green', 'line': {'color': 'black', 'width': 2.0}, 'size': 4, 'symbol': 'circle'},
                  'mode': 'markers',
                  'name': 'Stable',
                  'showlegend': True,
                  'type': 'scatterternary',
                  'uid': 'e2c50637-077d-4cd8-9bc1-64528f94e53b'},
                 {'a': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        2.0, 1.0, 28.0, 100.0, 2.0, 4.0, 4.0, 2.0, 4.0, 2.0, 2.0, 12.0,
                        4.0, 12.0, 40.0, 4.0, 4.0, 1.0, 2.0, 4.0, 2.0, 1.0, 2.0, 4.0,
                        8.0, 8.0, 1.0, 4.0, 12.0, 4.0, 12.0, 12.0, 2.0, 4.0, 4.0, 64.0,
                        4.0, 8.0, 16.0, 64.0, 16.0, 8.0, 16.0, 16.0, 4.0, 16.0, 4.0,
                        8.0, 32.0, 32.0, 32.0, 8.0, 4.0, 8.0, 8.0, 4.0, 8.0, 8.0, 6.0,
                        4.0, 6.0, 6.0, 6.0, 12.0, 12.0, 6.0, 24.0, 24.0, 12.0, 6.0,
                        12.0, 12.0, 12.0, 3.0, 12.0, 12.0, 12.0, 24.0, 24.0, 24.0,
                        24.0, 12.0, 12.0, 12.0, 12.0, 6.0, 12.0, 12.0, 12.0, 12.0,
                        12.0, 8.0, 4.0, 8.0, 5.0, 10.0, 20.0, 5.0, 5.0, 10.0, 20.0,
                        7.0, 14.0, 21.0, 14.0, 28.0, 16.0, 8.0, 16.0, 9.0, 9.0, 20.0,
                        10.0, 11.0, 12.0, 13.0, 13.0, 14.0, 15.0, 15.0, 17.0, 21.0,
                        21.0, 23.0, 23.0, 23.0, 23.0, 25.0, 32.0, 38.0, 41.0, 43.0],
                  'b': [0, 0, 0, 0, 0, 0, 0, 0, 4.0, 4.0, 2.0, 2.0, 4.0, 4.0, 2.0,
                        4.0, 8.0, 10.0, 14.0, 12.0, 8.0, 12.0, 16.0, 2.0, 2.0, 4.0,
                        2.0, 2.0, 8.0, 8.0, 71.0, 16.0, 2.0, 4.0, 2.0, 8.0, 2.0, 4.0,
                        4.0, 8.0, 48.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 8.0, 8.0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 4.0,
                        4.0, 4.0, 8.0, 8.0, 6.0, 6.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 4.0, 4.0, 4.0,
                        4.0, 4.0, 0, 1.0, 2.0, 0, 0, 0, 0, 0, 4.0, 8.0, 0, 0, 0, 6.0,
                        12.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 6.0, 0, 0, 0, 0, 0],
                  'c': [8.0, 2.0, 2.0, 2.0, 4.0, 8.0, 8.0, 4.0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0, 8.0,
                        4.0, 16.0, 4.0, 8.0, 8.0, 16.0, 36.0, 0, 0, 0, 0, 0, 0, 0, 2.0,
                        4.0, 2.0, 2.0, 12.0, 4.0, 12.0, 40.0, 4.0, 4.0, 2.0, 4.0, 8.0,
                        4.0, 2.0, 4.0, 8.0, 16.0, 16.0, 2.0, 8.0, 24.0, 8.0, 24.0,
                        24.0, 6.0, 24.0, 24.0, 96.0, 6.0, 12.0, 24.0, 96.0, 24.0, 12.0,
                        24.0, 24.0, 6.0, 24.0, 6.0, 12.0, 48.0, 48.0, 48.0, 12.0, 0, 0,
                        20.0, 14.0, 28.0, 28.0, 21.0, 18.0, 8.0, 8.0, 8.0, 16.0, 16.0,
                        8.0, 32.0, 32.0, 16.0, 8.0, 16.0, 16.0, 16.0, 4.0, 16.0, 16.0,
                        16.0, 32.0, 32.0, 32.0, 32.0, 16.0, 16.0, 16.0, 16.0, 0, 0, 0,
                        0, 0, 0, 10.0, 0, 0, 7.0, 14.0, 32.0, 8.0, 8.0, 0, 0, 8.0,
                        16.0, 27.0, 0, 0, 18.0, 9.0, 34.0, 10.0, 13.0, 22.0, 11.0,
                        12.0, 13.0, 14.0, 15.0, 15.0, 16.0, 16.0, 18.0, 23.0, 32.0,
                        25.0, 32.0, 32.0, 0, 32.0, 35.0, 39.0, 56.0, 64.0],
                  'cliponaxis': False,
                  'hoverinfo': 'text',
                  'hoverlabel': {'font': {'size': 14}},
                  'hovertext': [O&lt;sub&gt;2&lt;/sub&gt; (11) &lt;br&gt; 0.007 eV/atom (+0.007
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (42) &lt;br&gt; 0.113 eV/atom (+0.113
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (43) &lt;br&gt; 0.019 eV/atom (+0.019
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (46) &lt;br&gt; 0.016 eV/atom (+0.016
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (108) &lt;br&gt; 0.108 eV/atom
                                (+0.108 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (115) &lt;br&gt; 0.03 eV/atom
                                (+0.03 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (116) &lt;br&gt; 0.008 eV/atom
                                (+0.008 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (126) &lt;br&gt; 0.063
                                eV/atom (+0.063 eV/atom), C (1) &lt;br&gt; 0.141 eV/atom
                                (+0.141 eV/atom), C (2) &lt;br&gt; 0.001 eV/atom (+0.001
                                eV/atom), C (3) &lt;br&gt; 0.133 eV/atom (+0.133 eV/atom), C
                                (6) &lt;br&gt; 0.003 eV/atom (+0.003 eV/atom), C (34) &lt;br&gt;
                                0.003 eV/atom (+0.003 eV/atom), C (35) &lt;br&gt; 0.003
                                eV/atom (+0.003 eV/atom), C (36) &lt;br&gt; 0.005 eV/atom
                                (+0.005 eV/atom), C (37) &lt;br&gt; 0.004 eV/atom (+0.004
                                eV/atom), C (38) &lt;br&gt; 0.007 eV/atom (+0.007 eV/atom), C
                                (39) &lt;br&gt; 0.136 eV/atom (+0.136 eV/atom), C (40) &lt;br&gt;
                                0.135 eV/atom (+0.135 eV/atom), C (41) &lt;br&gt; 0.004
                                eV/atom (+0.004 eV/atom), C (44) &lt;br&gt; 0.137 eV/atom
                                (+0.137 eV/atom), C (45) &lt;br&gt; 0.136 eV/atom (+0.136
                                eV/atom), C (49) &lt;br&gt; 0.135 eV/atom (+0.135 eV/atom), C
                                (50) &lt;br&gt; 0.008 eV/atom (+0.008 eV/atom), C (99) &lt;br&gt;
                                0.003 eV/atom (+0.003 eV/atom), C (100) &lt;br&gt; 0.003
                                eV/atom (+0.003 eV/atom), C (101) &lt;br&gt; 0.003 eV/atom
                                (+0.003 eV/atom), C (106) &lt;br&gt; 0.003 eV/atom (+0.003
                                eV/atom), C (113) &lt;br&gt; 0.237 eV/atom (+0.237 eV/atom),
                                C (114) &lt;br&gt; 0.33 eV/atom (+0.33 eV/atom), C (118) &lt;br&gt;
                                0.104 eV/atom (+0.104 eV/atom), C (141) &lt;br&gt; 0.323
                                eV/atom (+0.323 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (10) &lt;br&gt;
                                -1.308 eV/atom (+0.007 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (32)
                                &lt;br&gt; -1.3 eV/atom (+0.015 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (51)
                                &lt;br&gt; -1.304 eV/atom (+0.01 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt;
                                (102) &lt;br&gt; -1.013 eV/atom (+0.301 eV/atom),
                                CO&lt;sub&gt;2&lt;/sub&gt; (110) &lt;br&gt; -1.282 eV/atom (+0.033
                                eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (120) &lt;br&gt; -1.298 eV/atom
                                (+0.016 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (134) &lt;br&gt; -1.292
                                eV/atom (+0.022 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (142) &lt;br&gt;
                                -1.304 eV/atom (+0.01 eV/atom),
                                C&lt;sub&gt;4&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (104) &lt;br&gt; -0.616 eV/atom
                                (+0.229 eV/atom), Fe (4) &lt;br&gt; 0.103 eV/atom (+0.103
                                eV/atom), Fe (5) &lt;br&gt; 0.173 eV/atom (+0.173 eV/atom),
                                Fe (145) &lt;br&gt; 0.199 eV/atom (+0.199 eV/atom), Fe (160)
                                &lt;br&gt; 0.247 eV/atom (+0.247 eV/atom), Fe (161) &lt;br&gt;
                                0.113 eV/atom (+0.113 eV/atom), Fe (162) &lt;br&gt; 0.022
                                eV/atom (+0.022 eV/atom), Fe (163) &lt;br&gt; 0.108 eV/atom
                                (+0.108 eV/atom), FeO (68) &lt;br&gt; 0.065 eV/atom (+0.344
                                eV/atom), FeO (79) &lt;br&gt; 0.108 eV/atom (+0.387 eV/atom),
                                FeO (80) &lt;br&gt; 0.157 eV/atom (+0.436 eV/atom), FeO (81)
                                &lt;br&gt; 0.156 eV/atom (+0.436 eV/atom), FeO (96) &lt;br&gt;
                                0.157 eV/atom (+0.436 eV/atom), FeO (123) &lt;br&gt; 0.055
                                eV/atom (+0.335 eV/atom), FeO (124) &lt;br&gt; 0.164 eV/atom
                                (+0.443 eV/atom), FeO (158) &lt;br&gt; 0.297 eV/atom (+0.577
                                eV/atom), FeO (165) &lt;br&gt; 0.055 eV/atom (+0.335
                                eV/atom), FeO (166) &lt;br&gt; 0.056 eV/atom (+0.335
                                eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (20) &lt;br&gt; 0.044 eV/atom
                                (+0.323 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (21) &lt;br&gt; -0.064
                                eV/atom (+0.215 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (67) &lt;br&gt;
                                0.045 eV/atom (+0.324 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (97)
                                &lt;br&gt; -0.082 eV/atom (+0.198 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt;
                                (107) &lt;br&gt; 0.024 eV/atom (+0.303 eV/atom),
                                FeO&lt;sub&gt;2&lt;/sub&gt; (119) &lt;br&gt; 0.051 eV/atom (+0.331
                                eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (122) &lt;br&gt; -0.028 eV/atom
                                (+0.251 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (131) &lt;br&gt; -0.051
                                eV/atom (+0.228 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (149) &lt;br&gt;
                                -0.046 eV/atom (+0.234 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (152)
                                &lt;br&gt; 0.044 eV/atom (+0.323 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt;
                                (170) &lt;br&gt; 0.05 eV/atom (+0.329 eV/atom),
                                FeO&lt;sub&gt;2&lt;/sub&gt; (173) &lt;br&gt; 0.042 eV/atom (+0.321
                                eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (175) &lt;br&gt; -0.051 eV/atom
                                (+0.228 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (176) &lt;br&gt; 0.04
                                eV/atom (+0.319 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (177) &lt;br&gt;
                                0.048 eV/atom (+0.327 eV/atom), FeCO&lt;sub&gt;3&lt;/sub&gt; (15)
                                &lt;br&gt; -0.89 eV/atom (+0.01 eV/atom),
                                Fe(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;2&lt;/sub&gt; (78) &lt;br&gt; -0.877 eV/atom
                                (+0.092 eV/atom), Fe(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;2&lt;/sub&gt; (167)
                                &lt;br&gt; -0.878 eV/atom (+0.091 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (7) &lt;br&gt; -0.292 eV/atom
                                (+0.043 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (17) &lt;br&gt;
                                -0.321 eV/atom (+0.014 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (25) &lt;br&gt; -0.28 eV/atom
                                (+0.055 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (31) &lt;br&gt;
                                -0.298 eV/atom (+0.037 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (54) &lt;br&gt; -0.276 eV/atom
                                (+0.059 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (59) &lt;br&gt;
                                -0.256 eV/atom (+0.079 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (70) &lt;br&gt; -0.279 eV/atom
                                (+0.056 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (74) &lt;br&gt;
                                -0.331 eV/atom (+0.004 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (95) &lt;br&gt; -0.288 eV/atom
                                (+0.047 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (111)
                                &lt;br&gt; -0.176 eV/atom (+0.159 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (125) &lt;br&gt; -0.295 eV/atom
                                (+0.04 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (133) &lt;br&gt;
                                -0.246 eV/atom (+0.089 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (148) &lt;br&gt; -0.233 eV/atom
                                (+0.102 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (156)
                                &lt;br&gt; -0.119 eV/atom (+0.217 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (157) &lt;br&gt; -0.136 eV/atom
                                (+0.199 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (159)
                                &lt;br&gt; -0.119 eV/atom (+0.216 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (171) &lt;br&gt; -0.142 eV/atom
                                (+0.193 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;C (8) &lt;br&gt; 0.041
                                eV/atom (+0.041 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;C (121) &lt;br&gt;
                                0.174 eV/atom (+0.174 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;CO&lt;sub&gt;5&lt;/sub&gt; (147) &lt;br&gt; -0.525 eV/atom
                                (+0.177 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (77) &lt;br&gt;
                                -0.822 eV/atom (+0.047 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (84) &lt;br&gt;
                                -0.811 eV/atom (+0.058 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (86) &lt;br&gt;
                                -0.804 eV/atom (+0.065 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (88) &lt;br&gt;
                                -0.823 eV/atom (+0.046 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;3&lt;/sub&gt; (85) &lt;br&gt;
                                -0.927 eV/atom (+0.038 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (14) &lt;br&gt; -0.19 eV/atom
                                (+0.13 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (16) &lt;br&gt;
                                -0.232 eV/atom (+0.087 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (22) &lt;br&gt; -0.239 eV/atom
                                (+0.08 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (24) &lt;br&gt;
                                -0.132 eV/atom (+0.187 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (26) &lt;br&gt; -0.233 eV/atom
                                (+0.086 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (47) &lt;br&gt;
                                -0.115 eV/atom (+0.204 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (53) &lt;br&gt; -0.24 eV/atom
                                (+0.079 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (56) &lt;br&gt;
                                -0.152 eV/atom (+0.167 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (69) &lt;br&gt; -0.184 eV/atom
                                (+0.135 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (72) &lt;br&gt;
                                -0.188 eV/atom (+0.131 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (73) &lt;br&gt; -0.234 eV/atom
                                (+0.085 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (75) &lt;br&gt;
                                -0.233 eV/atom (+0.086 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (76) &lt;br&gt; -0.237 eV/atom
                                (+0.082 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (127)
                                &lt;br&gt; -0.16 eV/atom (+0.16 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (128) &lt;br&gt; -0.236 eV/atom
                                (+0.084 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (129)
                                &lt;br&gt; -0.163 eV/atom (+0.156 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (130) &lt;br&gt; -0.25 eV/atom
                                (+0.069 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (132)
                                &lt;br&gt; -0.245 eV/atom (+0.074 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (135) &lt;br&gt; -0.245 eV/atom
                                (+0.075 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (136)
                                &lt;br&gt; -0.243 eV/atom (+0.076 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (139) &lt;br&gt; -0.241 eV/atom
                                (+0.079 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (143)
                                &lt;br&gt; -0.189 eV/atom (+0.13 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (164) &lt;br&gt; -0.242 eV/atom
                                (+0.077 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (168)
                                &lt;br&gt; -0.239 eV/atom (+0.08 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (169) &lt;br&gt; -0.24 eV/atom
                                (+0.079 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (12) &lt;br&gt; 0.046
                                eV/atom (+0.046 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (27) &lt;br&gt;
                                0.042 eV/atom (+0.042 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (48)
                                &lt;br&gt; 0.114 eV/atom (+0.114 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C
                                (138) &lt;br&gt; 0.238 eV/atom (+0.238 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;C (150) &lt;br&gt; 0.144 eV/atom (+0.144
                                eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (151) &lt;br&gt; 0.15 eV/atom
                                (+0.15 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;O&lt;sub&gt;5&lt;/sub&gt; (140) &lt;br&gt;
                                -0.126 eV/atom (+0.185 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;C (109)
                                &lt;br&gt; 0.275 eV/atom (+0.275 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;C
                                (112) &lt;br&gt; 0.111 eV/atom (+0.111 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (117) &lt;br&gt; -0.182 eV/atom
                                (+0.144 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (172)
                                &lt;br&gt; -0.133 eV/atom (+0.193 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (28) &lt;br&gt; -0.215 eV/atom
                                (+0.107 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (154)
                                &lt;br&gt; -0.192 eV/atom (+0.13 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (155) &lt;br&gt; -0.2 eV/atom
                                (+0.122 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt; (9) &lt;br&gt;
                                0.044 eV/atom (+0.044 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt; (52) &lt;br&gt; 0.046 eV/atom
                                (+0.046 eV/atom), Fe&lt;sub&gt;7&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (23) &lt;br&gt;
                                -0.046 eV/atom (+0.252 eV/atom),
                                Fe&lt;sub&gt;7&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (71) &lt;br&gt; -0.047 eV/atom
                                (+0.251 eV/atom), Fe&lt;sub&gt;7&lt;/sub&gt;O&lt;sub&gt;9&lt;/sub&gt; (63) &lt;br&gt;
                                -0.037 eV/atom (+0.277 eV/atom),
                                Fe&lt;sub&gt;7&lt;/sub&gt;C&lt;sub&gt;3&lt;/sub&gt; (13) &lt;br&gt; 0.063 eV/atom
                                (+0.063 eV/atom), Fe&lt;sub&gt;7&lt;/sub&gt;C&lt;sub&gt;3&lt;/sub&gt; (19) &lt;br&gt;
                                0.046 eV/atom (+0.046 eV/atom),
                                Fe&lt;sub&gt;8&lt;/sub&gt;O&lt;sub&gt;9&lt;/sub&gt; (64) &lt;br&gt; -0.042 eV/atom
                                (+0.254 eV/atom), Fe&lt;sub&gt;8&lt;/sub&gt;O&lt;sub&gt;9&lt;/sub&gt; (90) &lt;br&gt;
                                -0.03 eV/atom (+0.266 eV/atom),
                                Fe&lt;sub&gt;8&lt;/sub&gt;O&lt;sub&gt;17&lt;/sub&gt; (137) &lt;br&gt; 0.062 eV/atom
                                (+0.33 eV/atom), Fe&lt;sub&gt;9&lt;/sub&gt;O&lt;sub&gt;10&lt;/sub&gt; (87) &lt;br&gt;
                                -0.033 eV/atom (+0.261 eV/atom),
                                Fe&lt;sub&gt;9&lt;/sub&gt;O&lt;sub&gt;13&lt;/sub&gt; (174) &lt;br&gt; -0.124 eV/atom
                                (+0.206 eV/atom), Fe&lt;sub&gt;10&lt;/sub&gt;O&lt;sub&gt;11&lt;/sub&gt; (62)
                                &lt;br&gt; -0.019 eV/atom (+0.273 eV/atom),
                                Fe&lt;sub&gt;10&lt;/sub&gt;O&lt;sub&gt;11&lt;/sub&gt; (89) &lt;br&gt; -0.01 eV/atom
                                (+0.283 eV/atom), Fe&lt;sub&gt;11&lt;/sub&gt;O&lt;sub&gt;12&lt;/sub&gt; (57)
                                &lt;br&gt; -0.01 eV/atom (+0.281 eV/atom),
                                Fe&lt;sub&gt;12&lt;/sub&gt;O&lt;sub&gt;13&lt;/sub&gt; (91) &lt;br&gt; -0.004 eV/atom
                                (+0.286 eV/atom), Fe&lt;sub&gt;13&lt;/sub&gt;O&lt;sub&gt;14&lt;/sub&gt; (55)
                                &lt;br&gt; 0.003 eV/atom (+0.293 eV/atom),
                                Fe&lt;sub&gt;13&lt;/sub&gt;O&lt;sub&gt;15&lt;/sub&gt; (82) &lt;br&gt; -0.047 eV/atom
                                (+0.253 eV/atom), Fe&lt;sub&gt;14&lt;/sub&gt;O&lt;sub&gt;15&lt;/sub&gt; (92)
                                &lt;br&gt; 0.011 eV/atom (+0.3 eV/atom),
                                Fe&lt;sub&gt;15&lt;/sub&gt;O&lt;sub&gt;16&lt;/sub&gt; (60) &lt;br&gt; 0.043 eV/atom
                                (+0.332 eV/atom), Fe&lt;sub&gt;15&lt;/sub&gt;O&lt;sub&gt;16&lt;/sub&gt; (153)
                                &lt;br&gt; 0.046 eV/atom (+0.334 eV/atom),
                                Fe&lt;sub&gt;17&lt;/sub&gt;O&lt;sub&gt;18&lt;/sub&gt; (58) &lt;br&gt; 0.016 eV/atom
                                (+0.303 eV/atom), Fe&lt;sub&gt;21&lt;/sub&gt;O&lt;sub&gt;23&lt;/sub&gt; (66)
                                &lt;br&gt; -0.011 eV/atom (+0.281 eV/atom),
                                Fe&lt;sub&gt;21&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (29) &lt;br&gt; -0.223 eV/atom
                                (+0.109 eV/atom), Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;25&lt;/sub&gt; (61)
                                &lt;br&gt; -0.004 eV/atom (+0.287 eV/atom),
                                Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (93) &lt;br&gt; -0.171 eV/atom
                                (+0.154 eV/atom), Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (94)
                                &lt;br&gt; -0.237 eV/atom (+0.087 eV/atom),
                                Fe&lt;sub&gt;23&lt;/sub&gt;C&lt;sub&gt;6&lt;/sub&gt; (144) &lt;br&gt; 0.024 eV/atom
                                (+0.024 eV/atom), Fe&lt;sub&gt;25&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (146)
                                &lt;br&gt; -0.067 eV/atom (+0.246 eV/atom),
                                Fe&lt;sub&gt;32&lt;/sub&gt;O&lt;sub&gt;35&lt;/sub&gt; (98) &lt;br&gt; -0.01 eV/atom
                                (+0.282 eV/atom), Fe&lt;sub&gt;38&lt;/sub&gt;O&lt;sub&gt;39&lt;/sub&gt; (65)
                                &lt;br&gt; 0.042 eV/atom (+0.325 eV/atom),
                                Fe&lt;sub&gt;41&lt;/sub&gt;O&lt;sub&gt;56&lt;/sub&gt; (83) &lt;br&gt; -0.243 eV/atom
                                (+0.079 eV/atom), Fe&lt;sub&gt;43&lt;/sub&gt;O&lt;sub&gt;64&lt;/sub&gt; (30)
                                &lt;br&gt; -0.289 eV/atom (+0.045 eV/atom)],
                  'marker': {'color': [0.007, 0.113, 0.019, 0.016, 0.108, 0.03, 0.008,
                                       0.063, 0.141, 0.001, 0.133, 0.003, 0.003, 0.003,
                                       0.005, 0.004, 0.007, 0.136, 0.135, 0.004, 0.137,
                                       0.136, 0.135, 0.008, 0.003, 0.003, 0.003, 0.003,
                                       0.237, 0.33, 0.104, 0.323, 0.007, 0.015, 0.01,
                                       0.301, 0.033, 0.016, 0.022, 0.01, 0.229, 0.103,
                                       0.173, 0.199, 0.247, 0.113, 0.022, 0.108, 0.344,
                                       0.387, 0.436, 0.436, 0.436, 0.335, 0.443, 0.577,
                                       0.335, 0.335, 0.323, 0.215, 0.324, 0.198, 0.303,
                                       0.331, 0.251, 0.228, 0.234, 0.323, 0.329, 0.321,
                                       0.228, 0.319, 0.327, 0.01, 0.092, 0.091, 0.043,
                                       0.014, 0.055, 0.037, 0.059, 0.079, 0.056, 0.004,
                                       0.047, 0.159, 0.04, 0.089, 0.102, 0.217, 0.199,
                                       0.216, 0.193, 0.041, 0.174, 0.177, 0.047, 0.058,
                                       0.065, 0.046, 0.038, 0.13, 0.087, 0.08, 0.187,
                                       0.086, 0.204, 0.079, 0.167, 0.135, 0.131, 0.085,
                                       0.086, 0.082, 0.16, 0.084, 0.156, 0.069, 0.074,
                                       0.075, 0.076, 0.079, 0.13, 0.077, 0.08, 0.079,
                                       0.046, 0.042, 0.114, 0.238, 0.144, 0.15, 0.185,
                                       0.275, 0.111, 0.144, 0.193, 0.107, 0.13, 0.122,
                                       0.044, 0.046, 0.252, 0.251, 0.277, 0.063, 0.046,
                                       0.254, 0.266, 0.33, 0.261, 0.206, 0.273, 0.283,
                                       0.281, 0.286, 0.293, 0.253, 0.3, 0.332, 0.334,
                                       0.303, 0.281, 0.109, 0.287, 0.154, 0.087, 0.024,
                                       0.246, 0.282, 0.325, 0.079, 0.045],
                             'colorbar': {'len': 0.5,
                                          'thickness': 0.02,
                                          'thicknessmode': 'fraction',
                                          'title': {'text': 'Energy Above Hull&lt;br&gt;(eV/atom)'},
                                          'x': 0,
                                          'xpad': 0,
                                          'y': 1,
                                          'yanchor': 'top',
                                          'ypad': 0},
                             'colorscale': [[0.0, '#fff5e3'], [0.5, '#f24324'], [1.0,
                                            '#c40000']],
                             'line': {'color': 'black', 'width': 1},
                             'opacity': 0.8,
                             'size': 4,
                             'symbol': 'diamond'},
                  'mode': 'markers',
                  'name': 'Above Hull',
                  'showlegend': True,
                  'type': 'scatterternary',
                  'uid': '7245e698-c61e-4c91-b384-b500f4c68d27'}],
        'layout': {'autosize': True,
                   'coloraxis': {'colorbar': {'x': 1, 'y': 0.05, 'yanchor': 'top'}},
                   'height': 400,
                   'margin': {'b': 40, 'l': 20, 'r': 20, 't': 40},
                   'paper_bgcolor': 'rgba(0.9.,0.9.,0.9.,0.9)',
                   'plot_bgcolor': 'rgba(0,0,0,0)',
                   'showlegend': True,
                   'template': '...',
                   'ternary': {'aaxis': {'layer': 'below traces',
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'color': 'rgba(0,0,0,0)', 'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'Fe'}},
                               'baxis': {'layer': 'below traces',
                                         'min': 0,
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'color': 'rgba(0,0,0,0)', 'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'C'}},
                               'caxis': {'layer': 'below traces',
                                         'min': 0,
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'color': 'rgba(0,0,0,0)', 'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'O'}},
                               'sum': 1},
                   'width': 800}
    })


    Original Materials Project



    FigureWidget({
        'data': [{'a': [0, 1.0, None, 0, 2.0, None, 4.0, 0, None, 0, 1.0, None, 0, 0,
                        None, 2.0, 4.0, None, 0, 4.0, None, 4.0, 1.0, None, 0, 0, None,
                        2.0, 1.0, None],
                  'b': [4.0, 0, None, 4.0, 2.0, None, 0, 0, None, 4.0, 0, None, 4.0,
                        4.0, None, 2.0, 0, None, 4.0, 0, None, 0, 0, None, 4.0, 0,
                        None, 2.0, 0, None],
                  'c': [8.0, 0, None, 8.0, 6.0, None, 6.0, 8.0, None, 0, 0, None, 8.0,
                        0, None, 6.0, 6.0, None, 8.0, 6.0, None, 6.0, 0, None, 8.0,
                        8.0, None, 6.0, 0, None],
                  'hoverinfo': 'none',
                  'line': {'color': 'black', 'width': 1.5},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '9b995702-3cda-4bd0-9478-d1d7eaf11c07'},
                 {'a': [0, 0, 1.0],
                  'b': [4.0, 4.0, 0],
                  'c': [0, 8.0, 0],
                  'fill': 'toself',
                  'fillcolor': '#2E91E5',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'C—CO&lt;sub&gt;2&lt;/sub&gt;—Fe',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': 'db1fcf87-f06a-44ee-b6dc-a087a16b44e3'},
                 {'a': [0, 4.0, 0],
                  'b': [4.0, 0, 0],
                  'c': [8.0, 6.0, 8.0],
                  'fill': 'toself',
                  'fillcolor': '#E15F99',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'CO&lt;sub&gt;2&lt;/sub&gt;—Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;—O&lt;sub&gt;2&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '7a66ab6b-7312-448b-a43d-80f0eb4b8eb1'},
                 {'a': [1.0, 4.0, 2.0],
                  'b': [0, 0, 2.0],
                  'c': [0, 6.0, 6.0],
                  'fill': 'toself',
                  'fillcolor': '#1CA71C',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'Fe—Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;—FeCO&lt;sub&gt;3&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': 'cd1670d5-4418-455e-b8c2-ca92c082c5d6'},
                 {'a': [0, 1.0, 2.0],
                  'b': [4.0, 0, 2.0],
                  'c': [8.0, 0, 6.0],
                  'fill': 'toself',
                  'fillcolor': '#FB0D0D',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'CO&lt;sub&gt;2&lt;/sub&gt;—Fe—FeCO&lt;sub&gt;3&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '032458ee-4777-4c6f-926a-caea56809b2d'},
                 {'a': [0, 4.0, 2.0],
                  'b': [4.0, 0, 2.0],
                  'c': [8.0, 6.0, 6.0],
                  'fill': 'toself',
                  'fillcolor': '#DA16FF',
                  'hovertemplate': '&lt;extra&gt;&lt;/extra&gt;',
                  'line': {'width': 0},
                  'marker': {'size': 4},
                  'mode': 'lines',
                  'name': 'CO&lt;sub&gt;2&lt;/sub&gt;—Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;—FeCO&lt;sub&gt;3&lt;/sub&gt;',
                  'opacity': 0.15,
                  'showlegend': False,
                  'type': 'scatterternary',
                  'uid': '53dd7f46-4dc2-4797-94fd-968649e98550'},
                 {'a': [0, 1.0, 2.0, 4.0, 0, 0],
                  'b': [4.0, 0, 2.0, 0, 0, 4.0],
                  'c': [8.0, 0, 6.0, 6.0, 8.0, 0],
                  'cliponaxis': False,
                  'hoverinfo': 'text',
                  'hoverlabel': {'font': {'size': 14}},
                  'hovertext': [CO&lt;sub&gt;2&lt;/sub&gt; (18) &lt;br&gt; -1.312 eV/atom (Stable), Fe
                                (0) &lt;br&gt; 0.0 eV/atom (Stable), FeCO&lt;sub&gt;3&lt;/sub&gt; (15)
                                &lt;br&gt; -0.927 eV/atom (Stable),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (17) &lt;br&gt; -0.392 eV/atom
                                (Stable), O&lt;sub&gt;2&lt;/sub&gt; (11) &lt;br&gt; 0.0 eV/atom (Stable),
                                C (37) &lt;br&gt; 0.0 eV/atom (Stable)],
                  'marker': {'color': 'green', 'line': {'color': 'black', 'width': 2.0}, 'size': 4, 'symbol': 'circle'},
                  'mode': 'markers',
                  'name': 'Stable',
                  'showlegend': True,
                  'type': 'scatterternary',
                  'uid': 'e56abde6-da1f-4463-bceb-3d1305fce635'},
                 {'a': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        2.0, 1.0, 28.0, 100.0, 2.0, 4.0, 4.0, 2.0, 2.0, 8.0, 8.0, 4.0,
                        4.0, 4.0, 64.0, 8.0, 16.0, 16.0, 64.0, 16.0, 8.0, 16.0, 16.0,
                        4.0, 16.0, 4.0, 8.0, 32.0, 32.0, 32.0, 8.0, 4.0, 8.0, 8.0, 4.0,
                        8.0, 8.0, 6.0, 4.0, 6.0, 6.0, 6.0, 12.0, 12.0, 6.0, 24.0, 24.0,
                        12.0, 6.0, 12.0, 12.0, 12.0, 3.0, 12.0, 12.0, 12.0, 24.0, 24.0,
                        24.0, 24.0, 12.0, 12.0, 12.0, 12.0, 6.0, 12.0, 12.0, 12.0,
                        12.0, 12.0, 8.0, 4.0, 8.0, 5.0, 10.0, 20.0, 5.0, 5.0, 10.0,
                        20.0, 14.0, 14.0, 28.0, 16.0, 8.0, 9.0, 9.0, 20.0, 10.0, 11.0,
                        13.0, 21.0, 21.0, 23.0, 23.0, 23.0, 23.0, 41.0, 43.0],
                  'b': [0, 0, 0, 0, 0, 0, 0, 0, 4.0, 4.0, 2.0, 2.0, 4.0, 4.0, 2.0,
                        8.0, 10.0, 14.0, 12.0, 8.0, 12.0, 16.0, 2.0, 2.0, 4.0, 2.0,
                        4.0, 2.0, 8.0, 8.0, 71.0, 16.0, 2.0, 4.0, 2.0, 8.0, 2.0, 4.0,
                        4.0, 8.0, 48.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.0, 8.0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 4.0,
                        4.0, 4.0, 8.0, 8.0, 6.0, 6.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 4.0, 4.0, 4.0,
                        4.0, 4.0, 0, 1.0, 2.0, 0, 0, 0, 0, 0, 4.0, 8.0, 0, 6.0, 12.0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.0, 0, 0],
                  'c': [2.0, 2.0, 2.0, 2.0, 4.0, 8.0, 8.0, 4.0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0, 8.0,
                        4.0, 16.0, 4.0, 8.0, 8.0, 16.0, 36.0, 0, 0, 0, 0, 0, 0, 0, 4.0,
                        4.0, 16.0, 16.0, 8.0, 24.0, 24.0, 96.0, 12.0, 24.0, 24.0, 96.0,
                        24.0, 12.0, 24.0, 24.0, 6.0, 24.0, 6.0, 12.0, 48.0, 48.0, 48.0,
                        12.0, 0, 0, 20.0, 14.0, 28.0, 28.0, 21.0, 18.0, 8.0, 8.0, 8.0,
                        16.0, 16.0, 8.0, 32.0, 32.0, 16.0, 8.0, 16.0, 16.0, 16.0, 4.0,
                        16.0, 16.0, 16.0, 32.0, 32.0, 32.0, 32.0, 16.0, 16.0, 16.0,
                        16.0, 0, 0, 0, 0, 0, 0, 10.0, 0, 0, 7.0, 14.0, 32.0, 8.0, 8.0,
                        0, 0, 16.0, 0, 0, 18.0, 9.0, 10.0, 13.0, 22.0, 11.0, 12.0,
                        15.0, 23.0, 32.0, 25.0, 32.0, 32.0, 0, 56.0, 64.0],
                  'cliponaxis': False,
                  'hoverinfo': 'text',
                  'hoverlabel': {'font': {'size': 14}},
                  'hovertext': [O&lt;sub&gt;2&lt;/sub&gt; (42) &lt;br&gt; 0.124 eV/atom (+0.124
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (43) &lt;br&gt; 0.012 eV/atom (+0.012
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (46) &lt;br&gt; 0.01 eV/atom (+0.01
                                eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (105) &lt;br&gt; 0.002 eV/atom
                                (+0.002 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (108) &lt;br&gt; 0.126
                                eV/atom (+0.126 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (115) &lt;br&gt;
                                0.018 eV/atom (+0.018 eV/atom), O&lt;sub&gt;2&lt;/sub&gt; (116)
                                &lt;br&gt; 0.001 eV/atom (+0.001 eV/atom), O&lt;sub&gt;2&lt;/sub&gt;
                                (126) &lt;br&gt; 0.03 eV/atom (+0.03 eV/atom), C (1) &lt;br&gt;
                                0.162 eV/atom (+0.162 eV/atom), C (2) &lt;br&gt; 0.006
                                eV/atom (+0.006 eV/atom), C (3) &lt;br&gt; 0.136 eV/atom
                                (+0.136 eV/atom), C (6) &lt;br&gt; 0.001 eV/atom (+0.001
                                eV/atom), C (34) &lt;br&gt; 0.006 eV/atom (+0.006 eV/atom), C
                                (35) &lt;br&gt; 0.006 eV/atom (+0.006 eV/atom), C (36) &lt;br&gt;
                                0.01 eV/atom (+0.01 eV/atom), C (38) &lt;br&gt; 0.028 eV/atom
                                (+0.028 eV/atom), C (39) &lt;br&gt; 0.145 eV/atom (+0.145
                                eV/atom), C (40) &lt;br&gt; 0.144 eV/atom (+0.144 eV/atom), C
                                (41) &lt;br&gt; 0.008 eV/atom (+0.008 eV/atom), C (44) &lt;br&gt;
                                0.146 eV/atom (+0.146 eV/atom), C (45) &lt;br&gt; 0.143
                                eV/atom (+0.143 eV/atom), C (49) &lt;br&gt; 0.141 eV/atom
                                (+0.141 eV/atom), C (50) &lt;br&gt; 0.012 eV/atom (+0.012
                                eV/atom), C (99) &lt;br&gt; 0.004 eV/atom (+0.004 eV/atom), C
                                (100) &lt;br&gt; 0.008 eV/atom (+0.008 eV/atom), C (101) &lt;br&gt;
                                0.007 eV/atom (+0.007 eV/atom), C (103) &lt;br&gt; 0.008
                                eV/atom (+0.008 eV/atom), C (106) &lt;br&gt; 0.008 eV/atom
                                (+0.008 eV/atom), C (113) &lt;br&gt; 0.266 eV/atom (+0.266
                                eV/atom), C (114) &lt;br&gt; 0.299 eV/atom (+0.299 eV/atom),
                                C (118) &lt;br&gt; 0.116 eV/atom (+0.116 eV/atom), C (141)
                                &lt;br&gt; 0.291 eV/atom (+0.291 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt;
                                (10) &lt;br&gt; -1.305 eV/atom (+0.007 eV/atom),
                                CO&lt;sub&gt;2&lt;/sub&gt; (32) &lt;br&gt; -1.297 eV/atom (+0.015
                                eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (51) &lt;br&gt; -1.301 eV/atom
                                (+0.011 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (102) &lt;br&gt; -1.015
                                eV/atom (+0.297 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (110) &lt;br&gt;
                                -1.281 eV/atom (+0.031 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt; (120)
                                &lt;br&gt; -1.295 eV/atom (+0.017 eV/atom), CO&lt;sub&gt;2&lt;/sub&gt;
                                (134) &lt;br&gt; -1.294 eV/atom (+0.018 eV/atom),
                                CO&lt;sub&gt;2&lt;/sub&gt; (142) &lt;br&gt; -1.303 eV/atom (+0.009
                                eV/atom), C&lt;sub&gt;4&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (104) &lt;br&gt; -0.607
                                eV/atom (+0.236 eV/atom), Fe (4) &lt;br&gt; 0.098 eV/atom
                                (+0.098 eV/atom), Fe (5) &lt;br&gt; 0.148 eV/atom (+0.148
                                eV/atom), Fe (145) &lt;br&gt; 0.175 eV/atom (+0.175 eV/atom),
                                Fe (160) &lt;br&gt; 0.279 eV/atom (+0.279 eV/atom), Fe (161)
                                &lt;br&gt; 0.099 eV/atom (+0.099 eV/atom), Fe (162) &lt;br&gt;
                                0.017 eV/atom (+0.017 eV/atom), Fe (163) &lt;br&gt; 0.055
                                eV/atom (+0.055 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (21) &lt;br&gt;
                                -0.045 eV/atom (+0.282 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (97)
                                &lt;br&gt; -0.067 eV/atom (+0.26 eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt;
                                (131) &lt;br&gt; -0.032 eV/atom (+0.295 eV/atom),
                                FeO&lt;sub&gt;2&lt;/sub&gt; (149) &lt;br&gt; -0.026 eV/atom (+0.301
                                eV/atom), FeO&lt;sub&gt;2&lt;/sub&gt; (175) &lt;br&gt; -0.033 eV/atom
                                (+0.294 eV/atom), Fe(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;2&lt;/sub&gt; (78)
                                &lt;br&gt; -0.876 eV/atom (+0.107 eV/atom),
                                Fe(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;2&lt;/sub&gt; (167) &lt;br&gt; -0.874
                                eV/atom (+0.109 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt;
                                (7) &lt;br&gt; -0.273 eV/atom (+0.119 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (25) &lt;br&gt; -0.262 eV/atom
                                (+0.131 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (31) &lt;br&gt;
                                -0.277 eV/atom (+0.116 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (33) &lt;br&gt; -0.321 eV/atom
                                (+0.071 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (54) &lt;br&gt;
                                -0.167 eV/atom (+0.225 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (59) &lt;br&gt; -0.233 eV/atom
                                (+0.16 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (70) &lt;br&gt;
                                -0.259 eV/atom (+0.134 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (74) &lt;br&gt; -0.315 eV/atom
                                (+0.077 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (95) &lt;br&gt;
                                -0.269 eV/atom (+0.123 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (111) &lt;br&gt; -0.156 eV/atom
                                (+0.236 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (125)
                                &lt;br&gt; -0.28 eV/atom (+0.113 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (133) &lt;br&gt; -0.224 eV/atom
                                (+0.168 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (148)
                                &lt;br&gt; -0.215 eV/atom (+0.177 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (156) &lt;br&gt; -0.103 eV/atom
                                (+0.29 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (157) &lt;br&gt;
                                -0.119 eV/atom (+0.274 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (159) &lt;br&gt; -0.104 eV/atom
                                (+0.289 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;3&lt;/sub&gt; (171)
                                &lt;br&gt; -0.117 eV/atom (+0.275 eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;C
                                (8) &lt;br&gt; 0.059 eV/atom (+0.059 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C (121) &lt;br&gt; 0.178 eV/atom (+0.178
                                eV/atom), Fe&lt;sub&gt;2&lt;/sub&gt;CO&lt;sub&gt;5&lt;/sub&gt; (147) &lt;br&gt; -0.51
                                eV/atom (+0.227 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (77) &lt;br&gt;
                                -0.806 eV/atom (+0.088 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (84) &lt;br&gt;
                                -0.798 eV/atom (+0.096 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (86) &lt;br&gt;
                                -0.818 eV/atom (+0.077 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (88) &lt;br&gt;
                                -0.827 eV/atom (+0.067 eV/atom),
                                Fe&lt;sub&gt;2&lt;/sub&gt;(CO&lt;sub&gt;3&lt;/sub&gt;)&lt;sub&gt;3&lt;/sub&gt; (85) &lt;br&gt;
                                -0.922 eV/atom (+0.062 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (14) &lt;br&gt; -0.177 eV/atom
                                (+0.196 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (16) &lt;br&gt;
                                -0.27 eV/atom (+0.104 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (22) &lt;br&gt; -0.269 eV/atom
                                (+0.105 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (24) &lt;br&gt;
                                -0.119 eV/atom (+0.254 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (26) &lt;br&gt; -0.211 eV/atom
                                (+0.163 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (47) &lt;br&gt;
                                -0.108 eV/atom (+0.265 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (53) &lt;br&gt; -0.226 eV/atom
                                (+0.147 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (56) &lt;br&gt;
                                -0.14 eV/atom (+0.234 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (69) &lt;br&gt; -0.172 eV/atom
                                (+0.202 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (72) &lt;br&gt;
                                -0.175 eV/atom (+0.198 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (73) &lt;br&gt; -0.207 eV/atom
                                (+0.167 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (75) &lt;br&gt;
                                -0.207 eV/atom (+0.167 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (76) &lt;br&gt; -0.21 eV/atom
                                (+0.163 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (127)
                                &lt;br&gt; -0.162 eV/atom (+0.212 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (128) &lt;br&gt; -0.21 eV/atom
                                (+0.164 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (129)
                                &lt;br&gt; -0.153 eV/atom (+0.221 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (130) &lt;br&gt; -0.245 eV/atom
                                (+0.128 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (132)
                                &lt;br&gt; -0.237 eV/atom (+0.136 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (135) &lt;br&gt; -0.161 eV/atom
                                (+0.213 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (136)
                                &lt;br&gt; -0.233 eV/atom (+0.14 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (139) &lt;br&gt; -0.227 eV/atom
                                (+0.147 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (143)
                                &lt;br&gt; -0.177 eV/atom (+0.197 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (164) &lt;br&gt; -0.262 eV/atom
                                (+0.112 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (168)
                                &lt;br&gt; -0.257 eV/atom (+0.117 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;O&lt;sub&gt;4&lt;/sub&gt; (169) &lt;br&gt; -0.224 eV/atom
                                (+0.15 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (12) &lt;br&gt; 0.055
                                eV/atom (+0.055 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (27) &lt;br&gt;
                                0.056 eV/atom (+0.056 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (48)
                                &lt;br&gt; 0.123 eV/atom (+0.123 eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C
                                (138) &lt;br&gt; 0.277 eV/atom (+0.277 eV/atom),
                                Fe&lt;sub&gt;3&lt;/sub&gt;C (150) &lt;br&gt; 0.16 eV/atom (+0.16
                                eV/atom), Fe&lt;sub&gt;3&lt;/sub&gt;C (151) &lt;br&gt; 0.16 eV/atom
                                (+0.16 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;O&lt;sub&gt;5&lt;/sub&gt; (140) &lt;br&gt;
                                -0.135 eV/atom (+0.229 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;C (109)
                                &lt;br&gt; 0.278 eV/atom (+0.278 eV/atom), Fe&lt;sub&gt;4&lt;/sub&gt;C
                                (112) &lt;br&gt; 0.132 eV/atom (+0.132 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (117) &lt;br&gt; -0.159 eV/atom
                                (+0.222 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt; (172)
                                &lt;br&gt; -0.117 eV/atom (+0.265 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (28) &lt;br&gt; -0.188 eV/atom
                                (+0.19 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (154) &lt;br&gt;
                                -0.166 eV/atom (+0.212 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (155) &lt;br&gt; -0.169 eV/atom
                                (+0.208 eV/atom), Fe&lt;sub&gt;5&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt; (9) &lt;br&gt;
                                0.059 eV/atom (+0.059 eV/atom),
                                Fe&lt;sub&gt;5&lt;/sub&gt;C&lt;sub&gt;2&lt;/sub&gt; (52) &lt;br&gt; 0.061 eV/atom
                                (+0.061 eV/atom), Fe&lt;sub&gt;7&lt;/sub&gt;O&lt;sub&gt;8&lt;/sub&gt; (71) &lt;br&gt;
                                -0.075 eV/atom (+0.274 eV/atom),
                                Fe&lt;sub&gt;7&lt;/sub&gt;C&lt;sub&gt;3&lt;/sub&gt; (13) &lt;br&gt; 0.08 eV/atom
                                (+0.08 eV/atom), Fe&lt;sub&gt;7&lt;/sub&gt;C&lt;sub&gt;3&lt;/sub&gt; (19) &lt;br&gt;
                                0.064 eV/atom (+0.064 eV/atom),
                                Fe&lt;sub&gt;8&lt;/sub&gt;O&lt;sub&gt;9&lt;/sub&gt; (64) &lt;br&gt; -0.074 eV/atom
                                (+0.272 eV/atom), Fe&lt;sub&gt;8&lt;/sub&gt;O&lt;sub&gt;9&lt;/sub&gt; (90) &lt;br&gt;
                                -0.061 eV/atom (+0.285 eV/atom),
                                Fe&lt;sub&gt;9&lt;/sub&gt;O&lt;sub&gt;10&lt;/sub&gt; (87) &lt;br&gt; -0.072 eV/atom
                                (+0.272 eV/atom), Fe&lt;sub&gt;9&lt;/sub&gt;O&lt;sub&gt;13&lt;/sub&gt; (174)
                                &lt;br&gt; -0.095 eV/atom (+0.292 eV/atom),
                                Fe&lt;sub&gt;10&lt;/sub&gt;O&lt;sub&gt;11&lt;/sub&gt; (62) &lt;br&gt; -0.057 eV/atom
                                (+0.286 eV/atom), Fe&lt;sub&gt;10&lt;/sub&gt;O&lt;sub&gt;11&lt;/sub&gt; (89)
                                &lt;br&gt; -0.045 eV/atom (+0.298 eV/atom),
                                Fe&lt;sub&gt;11&lt;/sub&gt;O&lt;sub&gt;12&lt;/sub&gt; (57) &lt;br&gt; -0.044 eV/atom
                                (+0.297 eV/atom), Fe&lt;sub&gt;13&lt;/sub&gt;O&lt;sub&gt;15&lt;/sub&gt; (82)
                                &lt;br&gt; -0.07 eV/atom (+0.28 eV/atom),
                                Fe&lt;sub&gt;21&lt;/sub&gt;O&lt;sub&gt;23&lt;/sub&gt; (66) &lt;br&gt; -0.045 eV/atom
                                (+0.297 eV/atom), Fe&lt;sub&gt;21&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (29)
                                &lt;br&gt; -0.188 eV/atom (+0.201 eV/atom),
                                Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;25&lt;/sub&gt; (61) &lt;br&gt; -0.04 eV/atom
                                (+0.301 eV/atom), Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (93)
                                &lt;br&gt; -0.136 eV/atom (+0.244 eV/atom),
                                Fe&lt;sub&gt;23&lt;/sub&gt;O&lt;sub&gt;32&lt;/sub&gt; (94) &lt;br&gt; -0.216 eV/atom
                                (+0.165 eV/atom), Fe&lt;sub&gt;23&lt;/sub&gt;C&lt;sub&gt;6&lt;/sub&gt; (144)
                                &lt;br&gt; 0.052 eV/atom (+0.052 eV/atom),
                                Fe&lt;sub&gt;41&lt;/sub&gt;O&lt;sub&gt;56&lt;/sub&gt; (83) &lt;br&gt; -0.221 eV/atom
                                (+0.157 eV/atom), Fe&lt;sub&gt;43&lt;/sub&gt;O&lt;sub&gt;64&lt;/sub&gt; (30)
                                &lt;br&gt; -0.269 eV/atom (+0.122 eV/atom)],
                  'marker': {'color': [0.124, 0.012, 0.01, 0.002, 0.126, 0.018, 0.001,
                                       0.03, 0.162, 0.006, 0.136, 0.001, 0.006, 0.006,
                                       0.01, 0.028, 0.145, 0.144, 0.008, 0.146, 0.143,
                                       0.141, 0.012, 0.004, 0.008, 0.007, 0.008, 0.008,
                                       0.266, 0.299, 0.116, 0.291, 0.007, 0.015, 0.011,
                                       0.297, 0.031, 0.017, 0.018, 0.009, 0.236, 0.098,
                                       0.148, 0.175, 0.279, 0.099, 0.017, 0.055, 0.282,
                                       0.26, 0.295, 0.301, 0.294, 0.107, 0.109, 0.119,
                                       0.131, 0.116, 0.071, 0.225, 0.16, 0.134, 0.077,
                                       0.123, 0.236, 0.113, 0.168, 0.177, 0.29, 0.274,
                                       0.289, 0.275, 0.059, 0.178, 0.227, 0.088, 0.096,
                                       0.077, 0.067, 0.062, 0.196, 0.104, 0.105, 0.254,
                                       0.163, 0.265, 0.147, 0.234, 0.202, 0.198, 0.167,
                                       0.167, 0.163, 0.212, 0.164, 0.221, 0.128, 0.136,
                                       0.213, 0.14, 0.147, 0.197, 0.112, 0.117, 0.15,
                                       0.055, 0.056, 0.123, 0.277, 0.16, 0.16, 0.229,
                                       0.278, 0.132, 0.222, 0.265, 0.19, 0.212, 0.208,
                                       0.059, 0.061, 0.274, 0.08, 0.064, 0.272, 0.285,
                                       0.272, 0.292, 0.286, 0.298, 0.297, 0.28, 0.297,
                                       0.201, 0.301, 0.244, 0.165, 0.052, 0.157, 0.122],
                             'colorbar': {'len': 0.5,
                                          'thickness': 0.02,
                                          'thicknessmode': 'fraction',
                                          'title': {'text': 'Energy Above Hull&lt;br&gt;(eV/atom)'},
                                          'x': 0,
                                          'xpad': 0,
                                          'y': 1,
                                          'yanchor': 'top',
                                          'ypad': 0},
                             'colorscale': [[0.0, '#fff5e3'], [0.5, '#f24324'], [1.0,
                                            '#c40000']],
                             'line': {'color': 'black', 'width': 1},
                             'opacity': 0.8,
                             'size': 4,
                             'symbol': 'diamond'},
                  'mode': 'markers',
                  'name': 'Above Hull',
                  'showlegend': True,
                  'type': 'scatterternary',
                  'uid': '2b1538a5-3ea7-4999-95d4-dd411c69079f'}],
        'layout': {'autosize': True,
                   'coloraxis': {'colorbar': {'x': 1, 'y': 0.05, 'yanchor': 'top'}},
                   'height': 400,
                   'margin': {'b': 40, 'l': 20, 'r': 20, 't': 40},
                   'paper_bgcolor': 'rgba(0.9.,0.9.,0.9.,0.9)',
                   'showlegend': True,
                   'template': '...',
                   'ternary': {'aaxis': {'layer': 'below traces',
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'Fe'}},
                               'baxis': {'layer': 'below traces',
                                         'min': 0,
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'C'}},
                               'caxis': {'layer': 'below traces',
                                         'min': 0,
                                         'showgrid': True,
                                         'showticklabels': True,
                                         'tickfont': {'size': 7},
                                         'title': {'font': {'size': 24}, 'text': 'O'}},
                               'sum': 1},
                   'width': 800}
    })


**Atomsオブジェクトアクセス方法**

dfに入っているAtomsオブジェクトはプリミティブセルです。


```python
cands.loc[5].atoms
```




    Atoms(symbols='Fe', pbc=True, cell=[[2.577899217605591, 0.0, 0.0], [1.2889496088027954, 2.2325263023376465, 0.0], [1.2889496088027954, 0.744175374507904, 2.1048457622528076]], calculator=ASECalculator(...))



conventional は以下のとおりアクセスします。


```python
cs.get_conventional(cands.loc[[5]])
```




    [Atoms(symbols='Fe4', pbc=True, cell=[3.6457000078888497, 3.6457000078888497, 3.6457000078888497])]




```python
cs.df
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
      &lt;th&gt;material_id&lt;/th&gt;
      &lt;th&gt;pretty_formula&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;e_above_hull&lt;/th&gt;
      &lt;th&gt;is_hubbard&lt;/th&gt;
      &lt;th&gt;elements&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;mp-1&lt;/td&gt;
      &lt;td&gt;Cs&lt;/td&gt;
      &lt;td&gt;-0.856633&lt;/td&gt;
      &lt;td&gt;0.038770&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Cs]&lt;/td&gt;
      &lt;td&gt;(Atom('Cs', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;mp-2&lt;/td&gt;
      &lt;td&gt;Pd&lt;/td&gt;
      &lt;td&gt;-5.179882&lt;/td&gt;
      &lt;td&gt;0.000000&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Pd]&lt;/td&gt;
      &lt;td&gt;(Atom('Pd', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;mp-3&lt;/td&gt;
      &lt;td&gt;Cs&lt;/td&gt;
      &lt;td&gt;-1.598010&lt;/td&gt;
      &lt;td&gt;0.096397&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Cs]&lt;/td&gt;
      &lt;td&gt;(Atom('Cs', [3.6706762313842773, 3.38793683052...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;mp-4&lt;/td&gt;
      &lt;td&gt;Nd&lt;/td&gt;
      &lt;td&gt;-4.628184&lt;/td&gt;
      &lt;td&gt;0.139959&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Nd]&lt;/td&gt;
      &lt;td&gt;(Atom('Nd', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;mp-7&lt;/td&gt;
      &lt;td&gt;S&lt;/td&gt;
      &lt;td&gt;-24.431696&lt;/td&gt;
      &lt;td&gt;0.064501&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[S]&lt;/td&gt;
      &lt;td&gt;(Atom('S', [-2.1066184043884277, -2.1441278457...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139362&lt;/th&gt;
      &lt;td&gt;mvc-16821&lt;/td&gt;
      &lt;td&gt;CaCr2O4&lt;/td&gt;
      &lt;td&gt;-105.019547&lt;/td&gt;
      &lt;td&gt;0.107777&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Cr, O, Ca]&lt;/td&gt;
      &lt;td&gt;(Atom('Ca', [-1.5691694021224976, 1.5241652727...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139363&lt;/th&gt;
      &lt;td&gt;mvc-16832&lt;/td&gt;
      &lt;td&gt;V2ZnO4&lt;/td&gt;
      &lt;td&gt;-96.660103&lt;/td&gt;
      &lt;td&gt;0.161231&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Zn, O, V]&lt;/td&gt;
      &lt;td&gt;(Atom('V', [3.009662389755249, 0.0, 0.0], inde...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139364&lt;/th&gt;
      &lt;td&gt;mvc-16833&lt;/td&gt;
      &lt;td&gt;AlV2O4&lt;/td&gt;
      &lt;td&gt;-106.922005&lt;/td&gt;
      &lt;td&gt;0.207305&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[O, V, Al]&lt;/td&gt;
      &lt;td&gt;(Atom('Al', [-1.4248672723770142, -1.519302725...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139365&lt;/th&gt;
      &lt;td&gt;mvc-16834&lt;/td&gt;
      &lt;td&gt;MgMn2O4&lt;/td&gt;
      &lt;td&gt;-98.566040&lt;/td&gt;
      &lt;td&gt;0.105475&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Mg, O, Mn]&lt;/td&gt;
      &lt;td&gt;(Atom('Mg', [0.0, 0.0, 0.0], index=0), Atom('M...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139366&lt;/th&gt;
      &lt;td&gt;mvc-16835&lt;/td&gt;
      &lt;td&gt;Mg(FeO2)2&lt;/td&gt;
      &lt;td&gt;-89.655792&lt;/td&gt;
      &lt;td&gt;0.142727&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, Mg, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Mg', [0.0, 0.0, 0.0], index=0), Atom('M...&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;139367 rows × 7 columns&lt;/p&gt;
&lt;/div&gt;



DB全部にアクセスするには以下で可能です。


```python
cs.df
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
      &lt;th&gt;material_id&lt;/th&gt;
      &lt;th&gt;pretty_formula&lt;/th&gt;
      &lt;th&gt;energy&lt;/th&gt;
      &lt;th&gt;e_above_hull&lt;/th&gt;
      &lt;th&gt;is_hubbard&lt;/th&gt;
      &lt;th&gt;elements&lt;/th&gt;
      &lt;th&gt;atoms&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;mp-1&lt;/td&gt;
      &lt;td&gt;Cs&lt;/td&gt;
      &lt;td&gt;-0.856633&lt;/td&gt;
      &lt;td&gt;0.038770&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Cs]&lt;/td&gt;
      &lt;td&gt;(Atom('Cs', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;mp-2&lt;/td&gt;
      &lt;td&gt;Pd&lt;/td&gt;
      &lt;td&gt;-5.179882&lt;/td&gt;
      &lt;td&gt;0.000000&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Pd]&lt;/td&gt;
      &lt;td&gt;(Atom('Pd', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;mp-3&lt;/td&gt;
      &lt;td&gt;Cs&lt;/td&gt;
      &lt;td&gt;-1.598010&lt;/td&gt;
      &lt;td&gt;0.096397&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Cs]&lt;/td&gt;
      &lt;td&gt;(Atom('Cs', [3.6706762313842773, 3.38793683052...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;mp-4&lt;/td&gt;
      &lt;td&gt;Nd&lt;/td&gt;
      &lt;td&gt;-4.628184&lt;/td&gt;
      &lt;td&gt;0.139959&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[Nd]&lt;/td&gt;
      &lt;td&gt;(Atom('Nd', [0.0, 0.0, 0.0], index=0))&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;mp-7&lt;/td&gt;
      &lt;td&gt;S&lt;/td&gt;
      &lt;td&gt;-24.431696&lt;/td&gt;
      &lt;td&gt;0.064501&lt;/td&gt;
      &lt;td&gt;False&lt;/td&gt;
      &lt;td&gt;[S]&lt;/td&gt;
      &lt;td&gt;(Atom('S', [-2.1066184043884277, -2.1441278457...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139362&lt;/th&gt;
      &lt;td&gt;mvc-16821&lt;/td&gt;
      &lt;td&gt;CaCr2O4&lt;/td&gt;
      &lt;td&gt;-105.019547&lt;/td&gt;
      &lt;td&gt;0.107777&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Cr, O, Ca]&lt;/td&gt;
      &lt;td&gt;(Atom('Ca', [-1.5691694021224976, 1.5241652727...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139363&lt;/th&gt;
      &lt;td&gt;mvc-16832&lt;/td&gt;
      &lt;td&gt;V2ZnO4&lt;/td&gt;
      &lt;td&gt;-96.660103&lt;/td&gt;
      &lt;td&gt;0.161231&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Zn, O, V]&lt;/td&gt;
      &lt;td&gt;(Atom('V', [3.009662389755249, 0.0, 0.0], inde...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139364&lt;/th&gt;
      &lt;td&gt;mvc-16833&lt;/td&gt;
      &lt;td&gt;AlV2O4&lt;/td&gt;
      &lt;td&gt;-106.922005&lt;/td&gt;
      &lt;td&gt;0.207305&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[O, V, Al]&lt;/td&gt;
      &lt;td&gt;(Atom('Al', [-1.4248672723770142, -1.519302725...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139365&lt;/th&gt;
      &lt;td&gt;mvc-16834&lt;/td&gt;
      &lt;td&gt;MgMn2O4&lt;/td&gt;
      &lt;td&gt;-98.566040&lt;/td&gt;
      &lt;td&gt;0.105475&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Mg, O, Mn]&lt;/td&gt;
      &lt;td&gt;(Atom('Mg', [0.0, 0.0, 0.0], index=0), Atom('M...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;139366&lt;/th&gt;
      &lt;td&gt;mvc-16835&lt;/td&gt;
      &lt;td&gt;Mg(FeO2)2&lt;/td&gt;
      &lt;td&gt;-89.655792&lt;/td&gt;
      &lt;td&gt;0.142727&lt;/td&gt;
      &lt;td&gt;True&lt;/td&gt;
      &lt;td&gt;[Fe, Mg, O]&lt;/td&gt;
      &lt;td&gt;(Atom('Mg', [0.0, 0.0, 0.0], index=0), Atom('M...&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;139367 rows × 7 columns&lt;/p&gt;
&lt;/div&gt;


