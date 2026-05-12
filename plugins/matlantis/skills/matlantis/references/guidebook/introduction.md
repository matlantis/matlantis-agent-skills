# はじめに ​

Matlantisへようこそ! Matlantisは、汎用的で高速なニューラルネットワークポテンシャル（NNP）であるPFPに基づいた原子シミュレータを提供します。このシミュレータはブラウザのJupyterLab Notebook環境でコードを書くだけで利用することができ、ローカルマシンにソフトウェアをインストールする必要がありません。

## Getting started ​

Matlantisのnotebook環境を起動すると、Launcherページが自動的にロードされます。ここから、新しいnotebook、ターミナルウィンドウを作成することができます。また、Matlantisが提供するパッケージの使用方法を示す様々なチュートリアルや使用例も提供しています。

## Matlantisパッケージの全体構成 ​

notebook内でシミュレーションを実行する場合、計算のリクエストがPFP API serverに送られ、NNP計算が実行されます。Matlantisでは `pfp-api-client` というパッケージを提供しており、ユーザーはこのパッケージを通してサーバーからPFPによって計算されたエネルギーや力を得ることができるようになっています。

この他に、`matlantis-features` というパッケージがあり、`Features` という高レベルの関数が提供されています。これらのFeaturesは、内部で `pfp-api-client` を利用したシミュレーションを自動的に実行し、様々な物性値を簡単に取得できるようにするものです。

以下は、Matlantisが提供するパッケージ間の関係を示す図です。

notebook環境の詳しい仕様については [Notebook の仕様](/api/docs/ja/matlantis-guidebook/notebooks.html) をご覧ください。

## Documents and references ​

ユーザーのMatlantisへの理解を深める役に立つリファレンスを提供しています。これらのリファレンスはヘルプメニューからアクセスすることができます。

### About PFP ​

[About PFP](/api/docs/ja/pfp-description/2022.02/index.html) ではニューラルネットワークポテンシャル PFP のアーキテクチャと、 PFP の訓練に使用したデータセットについて詳しく説明しています。また、PFP の計算モード、適用範囲、および制限について説明し、さらに NNP や計算化学全般の広い分野での位置付けを説明します。

### About DFT-D3 ​

[About DFT-D3](/api/docs/ja/d3-description/2021.11/index.html) ではPFPにおけるDFT-D3分散力補正の実装と使用方法に関する詳細について説明しています。 `pfp-api-client` では分散力補正を含む計算モードも提供されています。

### Matlantis Guidebook ​

Matlantis Guidebook（本文書）では `pfp-api-client` と `matlantis-features` パッケージで提供される機能の使用方法とその背景の詳細や、Matlantis を利用する上で知っておくと便利な機能やライブラリついて説明しています。

### API References ​

[pfp-api-client](/api/docs/ja/pfp-api-client/latest/pfp_api_client.html) と [matlantis-features](/api/docs/ja/matlantis-features/latest/matlantis_features.html) のAPIリファレンスを提供しています。これらのリファレンスには、ライブラリで使用される各クラスや関数の引数の種類が記載されています。

### Tutorials and Examples ​

`pfp-api-client` と `matlantis-features` の使い方や、構造データの生成や外部ツールとのやり取りの方法を紹介するチュートリアルと使用例をnotebookとして提供しています。これらのnotebookは、Launcher画面から利用することができます。
