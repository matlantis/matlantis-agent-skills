---
name: mt-ssh
description: >
  Matlantis への SSH 接続設定を扱うスキルです。
  SSH, ssh-keygen, websocat, ssh_matlantis.sh, ssh_matlantis.bat,
  ProxyCommand, VSCode Remote-SSH, WinSCP,
  SSH接続, リモートアクセス, authorized_keys
  に関するコード生成時に使用してください。
---
# SSH 接続

## 概要

Matlantis への SSH 接続をセットアップすることで、ローカルのターミナルや VSCode から直接 Matlantis 環境にアクセスできます。SSH 接続は WebSocket 経由（`websocat` を使用）で確立されます。

### 前提条件

- admin ユーザーが Dashboard > Users で SSH 接続を許可済みであること
- Matlantis が起動していること
- Matlantis ターミナルで `/tmp/sshd.pid` が存在することを確認（なければ Matlantis を再起動）

### 接続の仕組み

SSH 接続は `ssh_matlantis.sh`（Mac/Linux）または `ssh_matlantis.bat`（Windows）スクリプトを ProxyCommand として使用し、`websocat` 経由で Matlantis の WebSocket エンドポイント (`wss://<tenant>.matlantis.com/nb/<user_id>/default/api/ssh-over-ws`) に接続します。

## ワークフロー

```
1. 前提ツールをインストールする（websocat 等）
   ↓
2. SSH 鍵を生成する
   ↓
3. 公開鍵を Matlantis にアップロードする
   ↓
4. ローカルの SSH を設定する（~/.ssh/config）
   ↓
5. 接続をテストする
```

## 実装パターン

### パターン A: Mac/Linux でのセットアップ

#### ステップ 1: websocat のインストール

```bash
brew install websocat
```

#### ステップ 2: SSH 鍵の生成

```bash
ssh-keygen -t rsa -b 4096
# VSCode 利用の場合はパスフレーズなし推奨
```

`~/.ssh/id_rsa`（秘密鍵）と `~/.ssh/id_rsa.pub`（公開鍵）が生成されます。

#### ステップ 3: 公開鍵のアップロード

`~/.ssh/id_rsa.pub` の内容を Matlantis の `/home/jovyan/.ssh/authorized_keys` に保存します。

Matlantis ターミナルで実行します。

```bash
cat /home/jovyan/id_rsa.pub >> /home/jovyan/.ssh/authorized_keys
```

**重要**: アップロード後、**Matlantis を再起動**してください。再起動しないと `Permission denied (publickey)` エラーになります。

#### ステップ 4: ssh_matlantis.sh の取得と設定

以下の URL から `ssh_matlantis.sh` をダウンロードします（ブラウザでクリック）。

```
https://redirect.matlantis.com/api/me/ssh
```

ダウンロードしたファイルに実行権限を付与します。

```bash
chmod +x ~/ssh_matlantis.sh
```

このスクリプトにはユーザー固有の情報（`MATLANTIS_USER_ID`, `NOTEBOOK_PRE_SHARED_KEY`, `MATLANTIS_DOMAIN`）が自動で埋め込まれています。

スクリプトの内容:

```bash
#! /bin/bash
set -e

MATLANTIS_USER_ID=<matlantis_user_id>
NOTEBOOK_PRE_SHARED_KEY=<notebook_pre_shared_key>
MATLANTIS_DOMAIN=<tenant_name>.matlantis.com
cd $(dirname $0)

if ! command -v websocat &> /dev/null; then
    if [ ! -f ./websocat ]; then
        echo "websocat not found, please download and save it to $PATH"
        exit 1
    fi
    WEBSOCAT=./websocat
else
    WEBSOCAT=websocat
fi

$WEBSOCAT --binary \
  -H="cookie: matlantis-notebook-pre-shared-key=${NOTEBOOK_PRE_SHARED_KEY}" \
  wss://${MATLANTIS_DOMAIN}/nb/${MATLANTIS_USER_ID}/default/api/ssh-over-ws
```

#### ステップ 5: ~/.ssh/config の設定

`~/.ssh/config` に以下を追記します（`[USERNAME]` は Mac のユーザー名に置き換え）。

```
Host matlantis
    ProxyCommand /Users/[USERNAME]/ssh_matlantis.sh
    User jovyan
    IdentityFile ~/.ssh/id_rsa
```

パーミッションを設定します。

```bash
chmod 600 ~/.ssh/config
```

#### ステップ 6: 接続確認

```bash
ssh -o StrictHostKeyChecking=accept-new matlantis
# または
ssh matlantis echo "SSH connection successful"
```

#### ファイル転送（scp）

```bash
scp local_file.txt matlantis:/home/jovyan/
```

### パターン B: Windows でのセットアップ

#### ステップ 1: SSH の動作確認

PowerShell で `ssh.exe` が使えるか確認します。

```powershell
cd $HOME
ssh.exe -V
```

バージョンが表示されれば OK です。

#### ステップ 2: SSH 鍵の生成

```powershell
ssh-keygen.exe -t rsa -b 4096
# VSCode 利用の場合はパスフレーズなし推奨
```

`C:\Users\[USERNAME]\.ssh\id_rsa` と `id_rsa.pub` が生成されます。

#### ステップ 3: 公開鍵のアップロード

`C:\Users\[USERNAME]\.ssh\id_rsa.pub` の内容を Matlantis の `/home/jovyan/.ssh/authorized_keys` に保存します。

Matlantis ターミナルで実行します。

```bash
cat /home/jovyan/id_rsa.pub >> /home/jovyan/.ssh/authorized_keys
rm /home/jovyan/id_rsa.pub
```

**重要**: アップロード後、**Matlantis を再起動**してください。

#### ステップ 4: ssh_matlantis.bat の取得

以下の URL から `ssh_matlantis.bat` をダウンロードします。

```
https://redirect.matlantis.com/api/me/ssh?script_type=bat
```

例えば `C:\Users\[USERNAME]\MyScripts\` フォルダを作成して保存します（場所は任意）。

ダウンロード時に「危険なファイル」の警告が出る場合はセキュリティチェックを外してダウンロードしてください。

bat ファイルの内容（ユーザー固有情報が自動入力されている）:

```bat
@echo off
setlocal enabledelayedexpansion

set MATLANTIS_USER_ID=<matlantis_user_id>
set NOTEBOOK_PRE_SHARED_KEY=<notebook_pre_shared_key>
set MATLANTIS_DOMAIN=<tenant_name>.matlantis.com
cd /d %~dp0

set BIN_NAME=websocat.x86_64-pc-windows-gnu.exe
set BIN_URL=https://github.com/vi/websocat/releases/latest/download/%BIN_NAME%
set WEBSOCAT_PATH=C:\Program Files\Websocat\%BIN_NAME%

where %BIN_NAME% >nul 2>&1
if %errorlevel% equ 0 (
    set WEBSOCAT=%BIN_NAME%
) else (
    if exist "%WEBSOCAT_PATH%" (
        set WEBSOCAT="%WEBSOCAT_PATH%"
    ) else (
        curl -s -L -o %BIN_NAME% %BIN_URL%
        set WEBSOCAT=%BIN_NAME%
    )
)

%WEBSOCAT% --binary -H="cookie: matlantis-notebook-pre-shared-key=%NOTEBOOK_PRE_SHARED_KEY%" wss://%MATLANTIS_DOMAIN%/nb/%MATLANTIS_USER_ID%/default/api/ssh-over-ws
```

websocat は bat ファイルと同じフォルダに自動ダウンロードされます。失敗した場合は `C:\Program Files\Websocat\` に手動で配置してください。

#### ステップ 5: ~/.ssh/config の設定

`C:\Users\[USERNAME]\.ssh\config` に以下を追記します（**全て絶対パスで指定**）。

```
Host matlantis
    ProxyCommand C:\Users\[USERNAME]\MyScripts\ssh_matlantis.bat
    User jovyan
    IdentityFile C:\Users\[USERNAME]\.ssh\id_rsa
```

**注意**: Windows では相対パスが動作しない場合があるため、ProxyCommand と IdentityFile は必ず絶対パスで指定してください。

#### ステップ 6: 接続確認

```powershell
ssh matlantis
```

### パターン C: VSCode Remote-SSH 統合

SSH 接続が確立できたら、VSCode の Remote-SSH 拡張機能を使って Matlantis に接続できます。

1. VSCode に **Remote-SSH** 拡張機能をインストールします
2. リモートエクスプローラーから **SSH > matlantis** に接続します
3. 必要に応じて以下の拡張機能もインストールします:
   - Jupyter Kernels
   - Python Environments

VSCode Remote-SSH を使うことで、ローカルの VSCode から直接 Matlantis 上のファイルを編集し、Notebook を実行できます。

### パターン D: WinSCP 統合（Windows）

Windows ユーザーは WinSCP を使ったファイル転送も可能です。

1. 新しいサイトを作成: ホスト名 `matlantis`、ユーザー名 `jovyan`
2. 設定 > プロキシ: ローカルプロキシを選択し、プロキシコマンドに `ssh_matlantis.bat` のパスを指定
3. 設定 > 認証: `id_rsa` 秘密鍵ファイルを指定（PuTTY 形式 `.ppk` に変換されます）

## ベストプラクティス

### SSH 鍵の管理

- VSCode Remote-SSH を使用する場合は、パスフレーズなしの鍵を推奨します（接続のたびにパスフレーズ入力が求められるため）
- 複数の Matlantis 環境に接続する場合は、環境ごとに別の Host エントリーを `~/.ssh/config` に追加してください

### 接続の安定性

- Matlantis が停止すると SSH 接続も切断されます。Matlantis の自動停止設定に注意してください
- 長時間接続を維持する場合は、`~/.ssh/config` に `ServerAliveInterval 60` を追加するとタイムアウトを防げます

```
Host matlantis
    ProxyCommand /Users/[USERNAME]/ssh_matlantis.sh
    User jovyan
    IdentityFile ~/.ssh/id_rsa
    ServerAliveInterval 60
    StrictHostKeyChecking accept-new
```

### セキュリティ

- `ssh_matlantis.sh` / `ssh_matlantis.bat` にはプリシェアードキーが含まれています。他者と共有しないでください
- `~/.ssh/config` のパーミッションは `600` に設定してください（Mac/Linux）

## よくあるエラーと対処

| エラー | 原因 | 対処 |
|---|---|---|
| `Permission denied (publickey)` | 公開鍵のアップロード後に Matlantis を再起動していない | Matlantis を再起動してください。再起動しないと `authorized_keys` が反映されません |
| `302 ERROR` | Matlantis が起動していない | ブラウザで Matlantis を起動してから SSH 接続を再試行してください |
| `Host key verification failed` | ホスト鍵が変更された（Matlantis 再起動等） | `ssh -o StrictHostKeyChecking=accept-new matlantis` で再接続するか、`~/.ssh/known_hosts` から該当行を削除してください。恒久対策として config に `StrictHostKeyChecking=accept-new` を追加してください |
| `websocat not found`（Mac） | websocat がインストールされていない | `brew install websocat` を実行してください |
| `websocat not found`（Windows） | bat ファイルによる自動ダウンロードが失敗 | websocat を手動でダウンロードし、`C:\Program Files\Websocat\` に配置してください |
| 接続がタイムアウトする | Matlantis の自動停止 / ネットワーク問題 | `ServerAliveInterval 60` を config に追加してください。Matlantis が起動中か確認してください |
| Windows でパスが認識されない | 相対パスを使用している | ProxyCommand と IdentityFile は**絶対パス**で指定してください |
| `/tmp/sshd.pid` が存在しない | SSH デーモンが起動していない | Matlantis を再起動してください |

## 関連ガイド

- **バックグラウンドジョブ** (`background-job/SKILL.md`): SSH 経由でのジョブ監視
- **統合ワークフロー** (`integrated-workflows/SKILL.md`): リモートからのパイプライン実行
