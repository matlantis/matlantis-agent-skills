# PFP Versions ​

このページでは、 Matlantis内で利用可能なPFPの各バージョンの仕様(サポート元素、 計算モード)などの情報を確認できます。

## Matlantisで提供されているPFPのバージョン ​

  * [v8.0.0](./pfp_versions/v8.0.0.html)
  * [v7.0.0](./pfp_versions/v7.0.0.html)
  * [v6.0.0](./pfp_versions/v6.0.0.html)
  * [v5.0.0](./pfp_versions/v5.0.0.html)
  * [v4.0.0](./pfp_versions/v4.0.0.html)
  * [v3.0.0](./pfp_versions/v3.0.0.html)
  * [v2.0.0](./pfp_versions/v2.0.0.html)
  * [v1.1.0](./pfp_versions/v1.1.0.html)
  * [v1.0.0](./pfp_versions/v1.0.0.html)
  * [v0.0.0](./pfp_versions/v0.0.0.html)



## 計算モードと各PFPバージョンの対応状況 ​

PFPは4種類のDFT計算条件による学習データによって訓練されています。これに対応して、計算モードと呼ばれる機能を提供しています 計算モードを切り替えることによってどちらの学習データに合った挙動をするかを選べます。 計算モードのデフォルトはPFPのバージョンごとに異なります。

#### 計算モードの一覧 ​

計算モード| 説明  
---|---  
PBE| PBE汎関数  
PBE_U| PBE汎関数(Hubbard補正あり)  
PBE_PLUS_D3| PBE汎関数+DFT-D3分散補正  
PBE_U_PLUS_D3| PBE汎関数(Hubbard補正あり)+DFT-D3分散補正  
R2SCAN| r2SCAN汎関数  
R2SCAN_PLUS_D3| r2SCAN汎関数+DFT-D3分散補正  
WB97XD| ωB97xd汎関数  
  
PBE、PBE_U、PBE_PLUS_D3、PBE_U_PLUS_D3、WB97XDが現在の計算モード名です。旧名称ではそれぞれCRYSTAL_U0、CRYSTAL、CRYSTAL_U0_PLUS_D3、CRYSTAL_PLUS_D3、MOLECULEに対応します。 PFPのバージョンによって利用可能な計算モードとデフォルト計算モードは異なります。

#### 各PFPバージョンで利用可能な計算モードの一覧とデフォルト計算モード ​

Version| PBE_U| WB97XD| PBE_U_PLUS_D3| PBE| PBE_PLUS_D3| R2SCAN| R2SCAN_PLUS_D3  
---|---|---|---|---|---|---|---  
v8.0.0| ✓| ✓| ✓| ✓ (default)| ✓| ✓| ✓  
v7.0.0| ✓| ✓| ✓| ✓ (default)| ✓| |   
v6.0.0| ✓| ✓| ✓| ✓ (default)| ✓| |   
v5.0.0| ✓| ✓| ✓| ✓ (default)| ✓| |   
v4.0.0| ✓| ✓| ✓| ✓ (default)| ✓| |   
v3.0.0| ✓ (default)| ✓| ✓| ✓| ✓| |   
v2.0.0| ✓ (default)| ✓| ✓| ✓| ✓| |   
v1.1.0| ✓ (default)| ✓| ✓| | | |   
v1.0.0| ✓ (default)| ✓| | | | |   
v0.0.0| ✓ (default)| ✓| | | | | 
