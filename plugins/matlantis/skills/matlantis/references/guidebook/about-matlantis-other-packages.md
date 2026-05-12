# Matlantis のその他のパッケージについて ​

Matlantis では、[pfp-api-client](./about-pfp-api-client.html) や [matlantis-features](./about-matlantis-features.html) 以外にも計算をサポートするライブラリを提供しています。

## matlantis-group-drive-client とは ​

`matlantis-group-drive-client` は Python コードから Group Drive を操作するためのライブラリです。

### サンプルコード ​

python
    
    
    from matlantis_group_drive_client import GroupDriveFileSystem
    
    drive_fs = GroupDriveFileSystem()
    drive = drive_fs.drive_list()[0]
    
    # Create directory on Group Drive
    dir_name = f"{drive['name']}/example_group_drive_client"
    drive_fs.mkdir(dir_name)
    
    # Upload a file to Group Drive
    drive_fs.upload(local_path="./test.txt", drive_path=f"{drive['name']}/example_group_drive_client/test.txt")
    
    # Download a file from Group Drive
    drive_fs.download(drive_path=f"{drive['name']}/example_group_drive_client/test.txt", local_path="./test.txt")
    
    # Use open()
    with drive_fs.open(f"{drive['name']}/example_group_drive_client/test.txt", "r") as f:
        print(f"open: {f.read()}")
    
    # Remove a file
    drive_fs.rm(f"{drive['name']}/example_group_drive_client/test.txt")

### チュートリアル ​

`matlantis-group-drive-client` のチュートリアルを提供しています。Example Launcher から利用可能です。

### API Reference ​

`matlantis-group-drive-client` の [API Reference](/api/docs/ja/matlantis-group-drive-client/0.1.0/index.html) を提供しています。このリファレンスにはライブラリで使用される各クラスや関数の引数の種類が記載されています。

## Experimental Packages とは ​

Matlantis では実験的なステータスのライブラリを Experimental Packages として提供しています。これらはバージョンの互換性を保証するものではなく、将来的に予告なく変更・削除される可能性があります。

Experimental Packages は Package Launcher から利用可能です。
