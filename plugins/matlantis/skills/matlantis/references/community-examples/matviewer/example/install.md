Copyright ENEOS, Corp. and  Preferred Computational Chemistry as contributors to Matlantis contrib project
# How to install matviewer

execute below commands. 　(Python3.8 or 3.9)

pfcc-extrasのインストール
　・Matlantis内のPackage Launcherから `pfcc-extras` をインストールしてください。
JSMEインストール ・・・下セルで実行できます
　・JSMEはサーバ通信なしで分子構造描画ができるBSD licenseソフトです。
autodEのインストール ・・・下セルで実行できます
　・SMILESからASE Atomsへの変換に、`autode` packageを使います。


```python
#JSMEはサーバ通信なしで分子構造描画ができJ-GLOBALでも採用されているBSD licenseソフトです。（https://jsme-editor.github.io/)
!git clone https://github.com/jsme-editor/jsme-editor.github.io.git
# Change from "master" branch to "2022-02-26" release
!cd jsme-editor.github.io &amp;&amp; git reset --hard 198ef9be65dd1e05846affff7a2dee17de3e267b
!mkdir -p ../matviewer/data/jsme/
!cp -r ./jsme-editor.github.io/dist/jsme/* ../matviewer/data/jsme/
!rm -rf ./jsme-editor.github.io
#autodEはsmiles から分子構造作成する際に利用します。
!pip install git+https://github.com/duartegroup/autodE.git
```





```python
!pip install  -e .. # To enable import from other folders
!pip install 'ase&gt;=3.23.0'
```



## Check install
**Please restart this kernel before executing following cells.**&lt;br&gt;
It is correctly installed when version is printed.


```python
import matviewer
print("matviewer version: ", matviewer.__version__)
```





    matviewer version:  0.1.2


##  For uninstall



```python
#!pip uninstall matviewer -y
```
