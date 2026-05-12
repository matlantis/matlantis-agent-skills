Copyright ENEOS, Corp., Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project

# moldraw

JSMEを用いて分子を描画しASE　atomsを返すUIツールです。  #eneos ibuka 2022/12/8  
#v0.1.0 2023/10/25 　Pandasバージョンアップ対応 パッケージ化対応,　引数からcalculator 設定可能に変更　3.8,3.9で動作確認


## セットアップ

### pfcc-extrasのインストール

Matlantis内のPackage Launcherから `pfcc-extras` をインストールしてください。

### JSMEインストール ・・・下セルで実行できます

JSMEはサーバ通信なしで分子構造描画ができJ-GLOBALでも採用されているBSD licenseソフトです。（https://jsme-editor.github.io/)&lt;br/&gt;
以下ではgitを用いてソースコードをダウンロードし使用しています。

### autodEのインストール　・・下セルで実行できます


SMILESからASE Atomsへの変換に、`autode` packageを使います。
 - https://github.com/duartegroup/autodE




```python
#JSMEはサーバ通信なしで分子構造描画ができJ-GLOBALでも採用されているBSD licenseソフトです。（https://jsme-editor.github.io/)
#以下ではgitを用いてソースコードをダウンロードしています。
#パッケージファイルにもコピーするため既に実施している場合も以下実行をお勧めします（2分以内です）
!git clone https://github.com/jsme-editor/jsme-editor.github.io.git
# Change from "master" branch to "2022-02-26" release
!cd jsme-editor.github.io &amp;&amp; git reset --hard 198ef9be65dd1e05846affff7a2dee17de3e267b
!mkdir -p ../moldraw/data/jsme/
!cp -r ./jsme-editor.github.io/dist/jsme/* ../moldraw/data/jsme/
#必要ない部分は削除します。
!rm -rf ./jsme-editor.github.io
```


```python
#autodEはsmiles から分子構造作成する際に利用します。
!pip install git+https://github.com/duartegroup/autodE.git
```


```python
!pip install ..  # To enable import from other folders
```


```python
#動作確認
import moldraw
moldraw.__version__ , moldraw.__file__ #パッケージとして保存されていることを確認（pythonバージョン別）
```




    ('0.1.0', '/home/jovyan/.py38/lib/python3.8/site-packages/moldraw/__init__.py')


