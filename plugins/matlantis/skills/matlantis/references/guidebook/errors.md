# エラーの一覧 ​

このページでは pfp-api と matlantis-features のエラーの内容に関する説明を確認できます。

## pfp-api ​

### エラーの一覧 ​

例外名| 発生原因| 対処方法  
---|---|---  
PFPAPIError| pfp-api のサーバー側でエラーが発生しました。| 以下の PFPAPIError の詳細 を参照してください。  
RetriesExceeded| 再実行数が Estimator で指定された max_retries を超えています。| サーバーの利用状況を確認するか、max_retries の値を高く設定してください。  
ConcurrentUseDetected| 一つのASECalculatorに対して並列アクセスしたことが検出されました。| スレッドごとに個別のASECalculatorインスタンスを生成し、使用してください。  
MultiCalculatorUseDetected| 一つのEstimatorを複数のASECalculator で共有して使用したことが検出されました。| ASECalculatorごとに個別のEstimatorインスタンスを生成し、使用してください。  
  
### PFPAPIError の詳細 ​

メッセージ内容| 発生原因| 対処方法  
---|---|---  
Atoms are far away from the primitive cell| 原子とセルとの距離が遠すぎます。| 入力のセルと原子座標を確認してください。  
Cell is too small| 入力構造のセルが pbc には小さすぎます。| 入力のセルを確認してください。  
Failed to decode JWT string| 入力メタデータの JWT 文字列のデコードに失敗しました。| 入力の EstimateRequest を確認してください。  
Illegal atomic number was detected| サポートされていない元素が指定されました。| 入力の元素情報を確認してください。  
Internal server error| サーバー側で予期せぬエラーが発生しました。| このエラーが続く場合は、お問い合わせフォームにてご連絡ください。  
Invalid input value is detected| 不正な入力が指定されました。| 入力を確認してください。  
No element information is provided| 指定された入力には元素情報が含まれていません。| 入力の元素情報を確認してください。  
Not implemented| gRPC で実装されていないリクエストが受け付けられました。| 入力の EstimateRequest を確認してください。  
Operation incomplete due to timeout| タイムアウトにより操作が失敗しました。| pfp-api-server の利用状況を確認してしばらく経ってから再度操作してください。  
Too big for GPU memory| GPU メモリが不足しています。| しばらく経ってから再度操作してください。このエラーが続く場合、入力構造を削減する必要がある可能性があります。  
Too many atoms| 入力構造の原子数が上限を超えています。| 入力の原子数を確認してください。  
Too many atoms considering periodic boundary conditions| 周期境界条件で推論を実行する際に必要な隣接セル内の仮想原子（ゴースト原子）を含めた入力構造の原子数が上限を超えています。| 入力の原子数を確認してください。  
Too many neighbors| 入力構造の近傍数が上限を超えています。これは密度の大きい(原子同士が接近しすぎている)系で原子数がある程度大きい場合に発生します。| 入力の近傍数を確認してください。近傍数を減らすためには原子数を減らすか、密度を小さくする必要があります。  
  
## matlantis-features ​

### エラーの一覧 ​

例外名| 発生原因| 対処方法  
---|---|---  
DeprecationWarning| 実行された操作方法が Deprecated となりました。将来のサポートが保証されません。| メッセージ内容にて推奨される操作方法に従って操作してください。  
MatlantisError| matlantis-features のサーバー側でエラーが発生しました。| 以下の MatlantisError の詳細 を参照してください。  
PermutationIndexError| フォノンの symmetry で原子位置の順列の検出が失敗しました。| 入力構造を確認してください。  
RetriesExceeded| 再実行数が Estimator で指定された max_retries を超えました。| サーバーの利用状況を確認するか、max_retries の値を高く設定してください。  
RNEMDIntegratorError| rNEMD extension を用いた Feature では、integrator 引数が LangevinIntegrator に設定されした。これにより、不正確な結果が得られると考えられます。| integrator の種類を変更してください。  
RNEMDTrajectoryError| rNEMD を用いた PostFeature では、rNEMD extension を使用せずに作成された Trajectory が渡されました。| rNEMD extension を使って Trajectory ファイルを作り直してください。  
UserWarning| 指定された入力では不正確な結果が得られる可能性があります。| メッセージ内容を参照して、入力の設定を確認してください。  
  
### MatlantisError の詳細 ​

メッセージ内容| 発生原因| 対処方法  
---|---|---  
Internal server error| サーバー側で予期せぬエラーが発生しました。| このエラーが続く場合は、お問い合わせフォームにてご連絡ください。  
Not implemented| gRPC で実装されていないリクエストが受け付けられました。| 入力の EstimateRequest を確認してください。  
Operation incomplete due to timeout| タイムアウトにより操作が失敗しました。| サーバーの利用状況を確認してしばらく経ってから再度操作してください。
