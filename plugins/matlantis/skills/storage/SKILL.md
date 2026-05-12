---
name: mt-storage
description: >
  Group Drive永続ストレージを扱うスキルです。
  GroupDriveFileSystem, matlantis-group-drive-client, Group Drive,
  fs.upload, fs.download, ファイル転送, ファイル共有,
  永続化, persistent storage, チーム共有,
  オブジェクトストレージ
  に関するコード生成時に使用してください。
---

# Group Drive ストレージ

## 概要

Matlantis の計算環境（コンテナ）は一時的なものであり、セッション終了後にデータが消える可能性があります（`/home/jovyan` 下はある程度保持されますが、確実な永続化やチーム共有には向きません）。**Group Drive** は永続化されたオブジェクトストレージであり、専用ライブラリ `matlantis-group-drive-client` を介してアクセスします。

本ガイドでは、Group Drive への接続、ファイルの転送、ディレクトリ管理、およびチーム運用におけるベストプラクティスを説明します。Group Drive はオブジェクトストレージベースのため、通常のファイルシステムとは異なる制約があります。これらの制約を理解した上で運用設計を行うことが重要です。

### Group Drive の主な特徴

- チーム間でのデータ共有が可能
- セッション終了後もデータが永続化される
- オブジェクトストレージベース（ランダムアクセス不可）
- Python 標準の `open()` では直接アクセスできない
- GUI ファイルブラウザ（Basic Launcher > Share > Group Drive）でも確認可能

## ワークフロー

Group Drive 運用の標準的なワークフローは以下の通りです。

1. **ライブラリのインストール**: `matlantis-group-drive-client` のインストール（初回のみ）
2. **Group Drive への接続**: クライアントを初期化し、利用可能なドライブ情報を取得する
3. **ドライブの選択**: 割り当てられたドライブリストから操作対象を選択する
4. **ファイル転送**: 計算結果のアップロード、または過去データのダウンロードを実行する
5. **ディレクトリ管理**: 必要なディレクトリの作成、ファイル一覧の確認、不要ファイルの削除を行う

## 実装パターン

### パターン A: 接続とドライブ選択

Group Drive に接続し、操作対象のドライブ情報を取得します。Matlantis では複数のドライブが割り当てられる場合があるため、リストから適切なドライブを選択する手順が必要です。

```python
# 必要に応じてインストール（初回のみ）
# !pip install matlantis-group-drive-client

from matlantis_group_drive_client import GroupDriveFileSystem
from typing import Tuple

def connect_to_drive() -> Tuple[GroupDriveFileSystem, str]:
    """
    Group Driveに接続し、最初のドライブのルートパスを取得する。

    Returns:
        fs: ファイルシステムクライアント
        drive_root: ドライブの識別名（パスの一部として使用）
    """
    fs = GroupDriveFileSystem()
    drives = fs.drive_list()

    if not drives:
        raise RuntimeError("No Group Drive assigned to this user.")

    # 利用可能なドライブの一覧を表示
    print("Available drives:")
    for i, drive in enumerate(drives):
        print(f"  [{i}] {drive['drive_name']} (ID: {drive['name']})")

    # 通常はリストの最初のドライブを使用
    # drive['name']: システム上の識別子（パスの一部、UUID形式の場合あり）
    # drive['drive_name']: GUI上の表示名
    target_drive = drives[0]
    drive_root = target_drive['name']

    print(f"\nConnected to: {target_drive['drive_name']} (ID: {drive_root})")
    return fs, drive_root
```

#### パスの構造について

`drive['name']` で取得される文字列は、一見すると UUID のような長い文字列の場合があります。この文字列をルートディレクトリとして、すべてのリモートパスを結合してください。

```python
# パス構築の例
drive_root = "a1b2c3d4-e5f6-7890-abcd-ef1234567890"
remote_path = f"{drive_root}/results/optimization/output.traj"
```

### パターン B: ファイルのアップロード

ローカル（JupyterLab 環境）のファイルを Group Drive にアップロードします。

```python
from pathlib import Path

def upload_results(
    fs: GroupDriveFileSystem,
    drive_root: str,
    local_path: str,
    remote_folder: str
):
    """
    ローカルファイルをGroup Driveにアップロードする。

    Args:
        fs: GroupDriveFileSystemインスタンス
        drive_root: ドライブのルートパス
        local_path: アップロードするローカルファイルのパス
        remote_folder: アップロード先のリモートフォルダ名
    """
    # リモートディレクトリパスの構築
    remote_dir = f"{drive_root}/{remote_folder}"

    # ディレクトリが存在しない場合は作成
    if not fs.exists(remote_dir):
        fs.mkdir(remote_dir)
        print(f"Created directory: {remote_dir}")

    filename = Path(local_path).name
    remote_path = f"{remote_dir}/{filename}"

    # アップロード実行
    fs.upload(local_path=local_path, drive_path=remote_path)
    print(f"Uploaded: {local_path} -> {remote_path}")
```

#### 複数ファイルのアップロード

```python
def upload_directory(
    fs: GroupDriveFileSystem,
    drive_root: str,
    local_dir: str,
    remote_folder: str
):
    """
    ローカルディレクトリ内の全ファイルをアップロードする。
    """
    local_path = Path(local_dir)
    for file_path in local_path.glob("*"):
        if file_path.is_file():
            upload_results(fs, drive_root, str(file_path), remote_folder)
```

### パターン C: ファイルのダウンロード

Group Drive からローカル環境にファイルをダウンロードします。

**重要**: Group Drive 上のパスを `pandas.read_csv("/path/to/drive/file.csv")` のように直接指定することはできません。必ず一度 `download` でローカルにコピーしてから読み込んでください。

```python
def download_data(
    fs: GroupDriveFileSystem,
    drive_path: str,
    local_path: str
):
    """
    Group Driveからファイルをダウンロードする。

    注意: Group Driveのパスはオブジェクトストレージ上の仮想パスであり、
    Python標準の open() では直接開けない。
    必ず download -> ローカル読み込みの2段階で処理する。
    """
    if not fs.exists(drive_path):
        raise FileNotFoundError(f"Remote file not found: {drive_path}")

    fs.download(drive_path=drive_path, local_path=local_path)
    print(f"Downloaded: {drive_path} -> {local_path}")
```

#### ダウンロード後のデータ読み込みパターン

```python
import pandas as pd
from ase.io import read

# CSV データの読み込み
download_data(fs, f"{drive_root}/data/results.csv", "/tmp/results.csv")
df = pd.read_csv("/tmp/results.csv")

# トラジェクトリの読み込み
download_data(fs, f"{drive_root}/structures/opt.traj", "/tmp/opt.traj")
atoms = read("/tmp/opt.traj")
```

### パターン D: ディレクトリ管理

ファイル一覧の取得、ディレクトリの作成・削除などの管理操作です。

```python
def list_files(fs: GroupDriveFileSystem, remote_dir: str):
    """
    指定ディレクトリ内のファイル一覧を表示する。
    """
    if not fs.exists(remote_dir):
        print(f"Directory not found: {remote_dir}")
        return []

    files = fs.ls(remote_dir)
    print(f"--- Files in {remote_dir} ---")
    for f in files:
        # fは辞書形式 {'name': 'path/to/file', ...}
        print(f"  {f['name']}")
    return files


def create_directory(fs: GroupDriveFileSystem, remote_dir: str):
    """
    リモートディレクトリを作成する。
    """
    if not fs.exists(remote_dir):
        fs.mkdir(remote_dir)
        print(f"Created: {remote_dir}")
    else:
        print(f"Already exists: {remote_dir}")
```

### パターン E: 上書き防止とタイムスタンプ退避

成果物の上書き事故を防ぐため、既存の出力先がある場合は日時サフィックスで退避してから新規作成する標準手順です。

```python
from datetime import datetime

def safe_create_output_dir(
    fs: GroupDriveFileSystem,
    drive_root: str,
    folder_name: str
) -> str:
    """
    出力ディレクトリが存在する場合はタイムスタンプ付きで退避し、
    新しいディレクトリを作成する。
    """
    remote_dir = f"{drive_root}/{folder_name}"

    if fs.exists(remote_dir):
        # 既存ディレクトリをタイムスタンプ付きで退避
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = f"{remote_dir}_backup_{timestamp}"
        # 注: Group Drive にはrename/moveがないため、
        # 必要に応じてコピー+削除で対応
        print(f"Warning: {remote_dir} already exists.")
        print(f"Creating new directory with timestamp suffix.")
        remote_dir = f"{drive_root}/{folder_name}_{timestamp}"

    fs.mkdir(remote_dir)
    print(f"Output directory: {remote_dir}")
    return remote_dir
```

### パターン F: 構造化された保存

ログファイルと構造ファイルを分離して保存する運用パターンです。

```python
def save_calculation_results(
    fs: GroupDriveFileSystem,
    drive_root: str,
    job_id: str,
    structure_files: list,
    log_files: list
):
    """
    計算結果を構造化して保存する。
    ログと構造ファイルを分離保存する。
    """
    base_dir = f"{drive_root}/calculations/{job_id}"

    # ディレクトリ構造の作成
    for subdir in ["structures", "log_files"]:
        remote_dir = f"{base_dir}/{subdir}"
        if not fs.exists(remote_dir):
            fs.mkdir(remote_dir)

    # 構造ファイルのアップロード
    for f in structure_files:
        filename = Path(f).name
        fs.upload(local_path=f, drive_path=f"{base_dir}/structures/{filename}")

    # ログファイルのアップロード
    for f in log_files:
        filename = Path(f).name
        fs.upload(local_path=f, drive_path=f"{base_dir}/log_files/{filename}")

    print(f"Results saved to: {base_dir}")
```

### パターン G: 冪等な再実行運用

Notebook バッチ実行時に「出力が既にあるジョブは SKIP」する運用を取り入れることで、失敗ジョブのみ再投入できるようにします。

```python
def run_with_skip(
    fs: GroupDriveFileSystem,
    drive_root: str,
    job_id: str,
    computation_fn,
    *args,
    **kwargs
):
    """
    出力が既に存在する場合はSKIPし、
    存在しない場合のみ計算を実行する冪等な実行パターン。
    """
    output_path = f"{drive_root}/results/{job_id}"

    if fs.exists(output_path):
        print(f"SKIPPED: Output already exists for {job_id}")
        return "SKIPPED"

    print(f"RUNNING: {job_id}")
    result = computation_fn(*args, **kwargs)

    # 結果の保存
    # ...

    print(f"COMPLETED: {job_id}")
    return "COMPLETED"
```

### パターン H: 命名規約による追跡

ジョブ ID やタイムスタンプを出力ファイル名に含め、Group Drive 上で「どの条件・どの実行」かをファイル名だけで追跡可能にします。

```python
from datetime import datetime

def generate_output_filename(
    prefix: str,
    job_id: str,
    extension: str = ".traj",
    include_timestamp: bool = True
) -> str:
    """
    追跡可能なファイル名を生成する。

    例: "opt_job001_20260316_143000.traj"
    """
    parts = [prefix, job_id]
    if include_timestamp:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        parts.append(timestamp)

    return "_".join(parts) + extension
```

## ベストプラクティス

### 直接アクセスの禁止を徹底する

Group Drive のパスは仮想的なものであり、Python 標準の `open()` や `pandas.read_csv()` では直接指定できません。必ず `download` でローカルにコピーしてから読み込んでください。

```python
# 誤ったパターン（動作しない）
# df = pd.read_csv(f"{drive_root}/data.csv")

# 正しいパターン
fs.download(drive_path=f"{drive_root}/data.csv", local_path="/tmp/data.csv")
df = pd.read_csv("/tmp/data.csv")
```

### ランダムアクセス・seek の不可を理解する

オブジェクトストレージベースのため、ファイルの途中だけを読み書きする `seek` 操作はサポートされていません。ファイル全体を一括でアップロード/ダウンロードする設計にしてください。

### 排他制御への対応

複数人が同時に同じファイル名でアップロードした場合、後勝ちで上書きされます。チーム運用時はファイル名にタイムスタンプやユーザー名を含めるなどの工夫を推奨します。

```python
import os
username = os.environ.get("USER", "unknown")
filename = f"result_{username}_{timestamp}.traj"
```

### 逐次・並列の使い分け

共有環境では、まず逐次実行で出力形式を固め、次に並列投入へ切り替えてください。これにより保存競合の原因切り分けが容易になります。

### 運用ログの可観測性

`READY/WAITING` 状態、CPU/メモリ/Token 率を定期ログ化し、保存先の負荷と投入レートを後追いで説明できる状態を保つことを推奨します。

### GUI での確認

Matlantis 画面左上の「Basic Launcher」>「Share」>「Group Drive」から、GUI ファイルブラウザでドライブの中身を確認できます。アップロード結果の確認や、チームメンバーの成果物の確認に便利です。

### 命名規約の固定化

ジョブ ID やタイムスタンプを出力ファイル名に含め、Group Drive 上で「どの条件・どの実行」かをファイル名だけで追跡可能にしてください。

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| `No Group Drive assigned to this user` | ドライブが割り当てられていない | Matlantis 管理者にドライブの割り当てを依頼する |
| `open()` でファイルが開けない | Group Drive のパスを直接指定している | `fs.download()` でローカルにコピーしてから読み込む |
| ファイルが見つからない | パスの構築ミス（drive_root の欠如など） | `fs.ls()` でディレクトリ内容を確認する。drive_root を先頭に含めているか確認する |
| アップロードしたファイルが上書きされた | 同名ファイルの同時アップロード（後勝ち） | ファイル名にタイムスタンプやユーザー名を含める |
| ディレクトリが作成できない | 親ディレクトリが存在しない | 階層的にディレクトリを作成する |
| ダウンロードが遅い | ファイルサイズが大きい | 必要なファイルのみを個別にダウンロードする。大規模データは分割保存を検討する |
| `seek` 操作でエラー | オブジェクトストレージはランダムアクセス非対応 | ファイル全体をダウンロードしてからローカルで処理する |
| 冪等再実行で SKIP されない | 出力パスの生成ロジックが一致していない | `run_with_skip` で使用するパス生成を統一する |
| パスが長すぎて読みにくい | UUID 形式の drive_root | drive_root を変数に格納し、一貫して使用する |

## 関連ガイド

- **物性計算** (property-analysis/SKILL.md): 大規模計算結果の保存先として Group Drive を使用
- **反応経路探索** (reaction/SKILL.md): 探索結果のバックアップ保存
- **外部データ連携** (external-data/SKILL.md): 取得した構造データの共有保存
- **可視化** (visualization/SKILL.md): 可視化画像の保存・共有
