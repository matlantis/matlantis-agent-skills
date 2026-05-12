## Matlantis Group Drive Client のインストール
`matlantis-group-drive-client` を pip install します。


```python
!pip install matlantis-group-drive-client
```


## Group Drive へのアクセス

下記のように `GroupDriveFileSystem` クラスを使って Group Drive にアクセスします。
`drive_fs.drive_list()` を使用すると、現在設定されている Drive の一覧を取得することができます。
一般的なテナントでは Drive は一つだけ設定されていますが、 Group Drive では、複数の Drive の利用が可能になっています。
`drive['name']` はディレクトリ名を、 `drive['drive_name']` はディレクトリ名とは別にドライブの表示名を指します。


```python
from matlantis_group_drive_client import GroupDriveFileSystem

drive_fs = GroupDriveFileSystem()
drive_list = drive_fs.drive_list()
print("Count of drives:", len(drive_list))

drive = drive_list[0]
print("Directory Name:", drive['name'])
print("Drive Name:", drive['drive_name'])
```


ここでは、一つ目の Drive(`drive`) に対して操作を行います。
Group Drive では、各 Drive はルート直下に置かれたディレクトリとして表現されます。そして、`drive["name"]` というのは、その Drive の名前であり、ディレクトリ名でもあります。
たとえば、 `n1q0el9v2d58m2jg` という名前の Drive がある場合、Group Drive のルート直下には `n1q0el9v2d58m2jg` というディレクトリが存在します。
（`n1q0el9v2d58m2jg` というのは一例で、実際には環境によって異なります。）
この Drive の直下にあるファイル/ディレクトリの一覧は `drive_fs.ls("n1q0el9v2d58m2jg")` を実行することにより取得できます。
次のコードでは、先ほど取得した `drive` 変数を使って、ファイルの一覧を取得して表示しています。


```python
files = drive_fs.ls(drive["name"])
for f in files:
    print(f["name"])
```


このように、 この Drive の中の `hoge.txt` というファイルにアクセスするためには、`"f{drive['name']}/hoge.txt"` というパスを参照する必要があります。

## ディレクトリの作成
作業用のディレクトリとして `example_group_drive_client` を `drive` の直下に作成します。
都合が悪い場合は実行前に適宜修正してください。



```python
dir_name = f"{drive['name']}/example_group_drive_client"
drive_fs.mkdir(dir_name)
```


## JupyterLab 上のファイルを Group Drive へアップロードする
JupyterLab 上のファイルシステムにあるデータを Group Drive にアップロードします。すでに存在するファイルパスを指定すると上書きされます。


```python
open("./README.md", "a").close()  # create example file

uploaded_readme_path = dir_name + "/README.md"
drive_fs.upload(local_path="./README.md", drive_path=uploaded_readme_path)
```


## Group Drive 上のファイルを JupyterLab へダウンロードする
Group Drive 上のファイルを JupyterLab のファイルシステムにダウンロードします。すでに存在するファイルパスを指定すると上書きされます。


```python
downloaded_readme_path = "./README_downloaded.md"
drive_fs.download(drive_path=uploaded_readme_path, local_path=downloaded_readme_path)
```


## ファイルの作成
次に、Group Drive 上にファイルを新規作成してみます。
次の例では、`f"{drive['name']}/example_group_drive_client/new_file.txt"` というパスにファイルを作成します。
もしこのパスで都合が悪い場合は実行前に適宜修正してください。


```python
new_file_path = dir_name + "/new_file.txt"

with drive_fs.open(new_file_path, "w") as f:
    f.write("new file")
```


以上で `new file` という文字列が書かれた `new_file.txt` というファイルが作成されました。
一度、Group Drive の UI から、該当のファイルがあることを確認してみてください。

## ファイルの閲覧
次に、先ほどの `new_file.txt` を python 上でメモリ上にダウンロードしてみます。


```python
# use cat() to get file contents as bytes
cat_response = drive_fs.cat(new_file_path)
print(f"cat: {cat_response}")

# use open() to get a file handle
with drive_fs.open(new_file_path, "r") as f:
    print(f"open: {f.read()}")
```


## ファイルの削除
基本的な操作の最後として、今作成したファイルを削除します。


```python
try:
    drive_fs.rm(new_file_path)
except Exception as e:
    # the file did not exist
    print(e)
drive_fs.ls(dir_name)
```


このように、 python 上から Group Drive のファイルを扱うことが可能になります。

次は、少し応用的な用途として、ase との連携について説明します。

## `ase.io` との連携
下記の例では、`Atoms()` によって新しい atom を作成し、それを `cif` ファイルとして Group Drive 上に保存しています。
そして、そのファイルをダウンロードし、出力しています。
このように、`ase.io.write()` や `ase.io.read()` に Group Drive の IO object をそのまま渡すことができます。

注意点としては、

- `mode` は `rb` や `wb` のように `b` をつける
- `format` は適切に指定する
    - 現在動作確認している format は `cif` のみ
    - `traj` は file seek を必要とするため対応が難しいことがわかっています

があります。


```python
import ase.io
from ase import Atoms

ase_file_path = dir_name + "/atom.cif"
atoms = Atoms("H2", [[0, 0, 0], [1.0, 0, 0]])

with drive_fs.open(ase_file_path, mode="wb") as f:
    ase.io.write(f, atoms, format="cif")

with drive_fs.open(ase_file_path, mode="rb") as f:
    atoms = ase.io.read(f, format="cif")
    print(atoms.symbols.formula)
```


## GroupDriveのディレクトリごとのディスク使用量の確認
実用的な応用例として、ディレクトリごとのディスク使用量とファイル数を表示する、duのようなプログラムを挙げます。
このcellがcopy and pasteで実行できるように、今まで定義されていた変数を再定義しています。
ここでは、Group driveのroot直下にあるファイル・ディレクトリについて表示するようにしています。pathを適宜書き換えることで他のディレクトリについて使用量を確認することができます。


```python
!pip install humanize tqdm matlantis-group-drive-client
```


```python
from matlantis_group_drive_client import GroupDriveFileSystem
from tqdm.notebook import tqdm
import humanize

drive_fs = GroupDriveFileSystem()
drive = drive_fs.drive_list()[0]
path = drive['name']

with tqdm() as pbar:
    for e in drive_fs.ls(path):
        pbar.update(1)
        count = 0
        size = 0
        for _, _, files in drive_fs.walk(e['name'], detail=True):
            for _, info in files.items():
                pbar.update(1)
                count += 1
                size += info['size']
        name = e['name'].split('/')[1]
        print(f"{humanize.naturalsize(size, gnu=True):5} {count:5} {name}")
```


## Group Drive の特徴と制限について
Group Drive は高い可用性やスケーラビリティを実現するために、一般的なファイルシステムとは異なる性質を持っています。  
そのため、 NFS のようなネットワークを介して接続するファイルシステムとは使い勝手が異なることがあります。  
例えば、ランダムアクセスを必要とする file seek などは対応していません。
また、Group Drive 上のファイルパスをそのまま python コードやサードパーティーライブラリなどで使用することはできません。Group Drive 上のファイルを操作したい場合は、このチュートリアルで行っているようにファイルハンドラーを使用してください。
