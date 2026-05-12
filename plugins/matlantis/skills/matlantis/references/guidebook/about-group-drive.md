# Group Drive について ​

## Group Drive とは ​

あなたが所属するグループのユーザー間で、ファイルを自由に共有することができるストレージです。

Group Drive は Notebook, Dashboard どちらからでも利用することができます。

また、[matlantis-group-drive-client](./about-matlantis-other-packages.html#matlantis-group-drive-client) というライブラリを使うことで、Python コードから Group Drive を操作することができます。

## 主な使い方 ​

### 計算結果や成果をメンバー間で共有して活用 ​

  * 作業ファイルや計算結果を Group Drive で管理することで、メンバー間での分担作業や共同作業が行いやすくなります。
  * よく使う関数などを python ファイルに書いて共有することもできます。



### チームの成果を整理して保管 ​

  * チーム毎やプロジェクト毎にフォルダを作ることで、成果を一元管理することができます。



### 他の人のスクリプトから学ぶ、フィードバックを返す ​

  * 共有された便利なスクリプトや計算方法などを学ぶことができます。
  * Notebook などをダウンロードすることなくプレビューができるので、効率的にフィードバックをすることが可能です。



## 各機能詳細 ​

### ローカルファイルをアップロード ​

ユーザーのローカルマシンからファイルをアップロードすることができます。アップロードには2通りの方法があります。

  * 該当ファイルを Group Drive にドラッグ＆ドロップ 
    *   * アップロードボタンを押して該当ファイルを選択 
    * 


### Notebook の個人のワークスペースと Group Drive 間でファイルをコピー ​

ファイルをドラッグ＆ドロップすることで、Notebook の個人のワークスペースと Group Drive 間でファイルをコピーすることができます。

### ファイルプレビュー ​

Group Drive 上のファイルをプレビューすることができます。プレビューはファイルをダブルクリックするか、ファイル上で「右クリック

> View」で実行できます。

### リンクを使って共有 ​

Group Drive 上のファイルを他のユーザーに見てもらいたい際に「Copy sharable link > Copy link for Notebook」を選択することで共有用 URL をコピーすることができます。その URL を他のユーザーに共有することで、対象のファイルを直接を開いてもらうことができます。

### Python コードからのアクセス ​

[matlantis-group-drive-client](./about-matlantis-other-packages.html#matlantis-group-drive-client) を使うことで、Python コードから Group Drive 上のファイルに直接アクセスすることができます。 このライブラリは [fsspec](https://filesystem-spec.readthedocs.io/en/latest/index.html) というファイルシステムの仕様に従っており、 下記のように Group Drive 上のファイル操作を可能にします。

python
    
    
    from matlantis_group_drive_client import GroupDriveFileSystem
    drive_fs = GroupDriveFileSystem()
    drive = drive_fs.drive_list()[0]
    with drive_fs.open(f"{drive['name']}/new_file.txt", "w") as f:
        f.write("file contents")

詳しい使い方については、Example Launcher にある “How to access Group Drive from Python” というチュートリアルを参照ください。

## ヘルプドキュメントについて ​

Group Drive のヘルプドキュメントは以下のボタンをクリックすることで確認できます。

## 注意点 ​

  * Group Drive ではファイルを直接編集することはできません。編集したい場合は一度 Notebook の個人のワークスペースにコピーしてきて編集するか、ローカルマシンにダウンロードした上で編集してください。
  * ローカルマシンから Group Drive に直接フォルダをアップロードすることはできません。
  * Group Drive にあるフォルダを直接ローカルマシンにダウンロードすることはできません。ダウンロードしたい場合は該当のフォルダを一度 Notebook の個人のワークスペースにコピーしてきて、そこからダウンロードしてください。
  * 操作可能なファイルサイズの上限は**125 GiB** です。
  * Group Drive では Drive 毎に容量が制限されています。ご利用中の容量と最大使用可能量は Group Drive 右上の UI から確認することができます。容量追加をご希望のお客様はお問い合わせください。 
  * 表示されている利用量は基本的に即時反映されますが、まれに反映が遅れることがあり、実際の利用量と異なる場合があります。


