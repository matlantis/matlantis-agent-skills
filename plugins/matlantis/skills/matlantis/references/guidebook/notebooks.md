# Notebook の仕様 ​

このページでは notebook や各ユーザーに割り振られた Workspace の容量などを確認できます。また、PFPの推論のための計算リソースは各テナントに別途割り振られています。

## Notebookのスペック ​

Notebook のスペックは以下となります。

#### Notebook のスペック一覧表 ​

CPU| 2 コア  
---|---  
RAM| 8 GiB  
Workspace の容量| 100 GiB  
  
Notebook 環境は以下の条件のいずれかを満たす際、自動的にシャットダウンされます。

  * Notebook 環境を60日間連続で起動している
  * 30日以上ブラウザからのアクセスがない
  * 以下の3つの条件をすべて満たす 
    * 90分以上ブラウザからのアクセスがない
    * 全カーネルが 90 分間以上 Idle 状態
    * Notebook と Kernel 以外に python, papermill, pt_main_thread, lmp_serial, lmp, lmp_mpi プロセスが存在しない



## 接続可能ドメイン ​

接続可能ドメインは以下となります。

  * api.business.githubcopilot.com:443
  * api.enterprise.githubcopilot.com:443
  * api.github.com:443
  * api.githubcopilot.com:443
  * api.materialsproject.org:443
  * api.wandb.ai:443
  * az764295.vo.msecnd.net:443
  * codeload.github.com:443
  * conda.anaconda.org:443
  * contribs-api.materialsproject.org:443
  * copilot-proxy.githubusercontent.com:443
  * copilot-telemetry.githubusercontent.com:443
  * download-cdn.jetbrains.com:443
  * download.jetbrains.com:443
  * download.pytorch.org:443
  * downloads.marketplace.jetbrains.com:443
  * figshare.com:443
  * files.pythonhosted.org:443
  * github-cloud.githubusercontent.com:443
  * github.com:443
  * gitlab.com:443
  * icsd.fiz-karlsruhe.de:443
  * json.schemastore.org:443
  * legacy.materialsproject.org:443
  * marketplace.visualstudio.com:443
  * materialsproject-parsed.s3.amazonaws.com:443
  * materialsproject.org:443
  * ndownloader.figshare.com:443
  * objects.githubusercontent.com:443
  * open-catalyst-api.metademolab.com:443
  * oqmd.org:443
  * origin-tracker.githubusercontent.com:443
  * pkg-containers.githubusercontent.com:443
  * plugins.jetbrains.com:443
  * providers.optimade.org:443
  * proxy.business.githubcopilot.com:443
  * proxy.enterprise.githubcopilot.com:443
  * pubchem.ncbi.nlm.nih.gov:443
  * pypi.nvidia.com:443
  * pypi.org:443
  * pypi.python.org:443
  * release-assets.githubusercontent.com:443
  * repo.anaconda.com:443
  * resources.jetbrains.com:443
  * telemetry.business.githubcopilot.com:443
  * update.code.visualstudio.com:443
  * vscode.blob.core.windows.net:443
  * www.crystallography.net:443
  * www.google-analytics.com:443
  * www.jetbrains.com:443
  * www.schemastore.org:443



## マルチプロセス処理における Estimator の使用について ​

マルチプロセス処理で Estimator を使用する場合は、 joblib.Parallel で backend に loky を指定する、 もしくは、何も指定しないでデフォルトのままでお使いください（デフォルト backend は loky ）。

Pythonの標準ライブラリの multiprocessing を使用する場合は、メインプロセスで Estimator を初期化した状態で、 マルチプロセス処理で生成した子プロセスで Estimator を使用することはできません。 ただし、以下の例に示すように Estimator を一切初期化していない状態で子プロセスを生成した場合は、子プロセスで Estimator を使用することができます。 なお、 joblib.Parallel を使用した場合でも backend に multiprocessing を指定した場合は同様の制限を受けます。 マルチスレッド処理や新しいプロセスを生成しない処理の場合は上記の制限を受けません。

使用例:

python
    
    
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    
    import ase
    from ase.build import bulk
    from pfp_api_client import ASECalculator, Estimator
    
    
    def calculate_energy():
        with ASECalculator(Estimator()) as calculator:
            atoms = bulk("Si")
            atoms.calc = calculator
            return atoms.get_potential_energy()
    
    NUM_ACCESS=2
    
    # OK: Running joblib.Parallel with default loky backend
    from joblib import Parallel, delayed
    energy_list = Parallel(NUM_ACCESS)(
        delayed(calculate_energy)() for _ in range(NUM_ACCESS))
    
    # OK: No Estimator is instantiated in this process
    with ProcessPoolExecutor(max_workers=NUM_ACCESS) as executor:
        futures = [executor.submit(calculate_energy) for _ in range(NUM_ACCESS)]
        energy_list = [future.result() for future in futures]
    
    # OK: Running multi-threaded computations
    # even after a Estimator is instantiated by `calculate_energy()` call in this process.
    _ = calculate_energy()
    energy_list = Parallel(NUM_ACCESS, backend="threading")(
        delayed(calculate_energy)() for _ in range(NUM_ACCESS))
    
    with ThreadPoolExecutor(max_workers=NUM_ACCESS) as executor:
        futures = [executor.submit(calculate_energy) for _ in range(NUM_ACCESS)]
        energy_list = [future.result() for future in futures]

エラーとなる例:

python
    
    
    # NOT OK: A Estimator is instantiated by `calculate_energy()` call in this process
    # before multi-process execution
    _ = calculate_energy()
    with ProcessPoolExecutor(max_workers=NUM_ACCESS) as executor:
        futures = [executor.submit(calculate_energy) for _ in range(NUM_ACCESS)]
        energy_list = [future.result() for future in futures]
    
    # NOT OK: a Estimator is instantiated in this process
    # before running joblib.Parallel with multiprocessing backend
    energy_list = Parallel(NUM_ACCESS, backend="multiprocessing")(
        delayed(calculate_energy)() for _ in range(NUM_ACCESS))

各プロセスで繰り返し PFP の計算をおこなう場合、同一プロセス内で `Estimator` オブジェクトを共有できます。 これにより、 `Estimator` オブジェクトの作成回数が削減され計算パフォーマンスの悪化を防げます。

python
    
    
    import ase
    from ase.build import bulk
    from pfp_api_client import ASECalculator, Estimator
    
    
    # Not sharing an estimator.
    def calculate_energy():
        energies_by_repeats = {}
        for n_repeats in range(1, 20):
            atoms = bulk("Si") * (n_repeats, n_repeats, n_repeats)
            with ASECalculator(Estimator()) as calculator:
                atoms.calc = calculator
                energy = atoms.get_potential_energy()
                energies_by_repeats[n_repeats] = energy
        return energies_by_repeats
    
    
    # Sharing an estimator during for-loop.
    def calculate_energy_share():
        energies_by_repeats = {}
        with ASECalculator(Estimator()) as calculator:
            for n_repeats in range(1, 20):
                atoms = bulk("Si") * (n_repeats, n_repeats, n_repeats)
                atoms.calc = calculator
                energy = atoms.get_potential_energy()
                energies_by_repeats[n_repeats] = energy
        return energies_by_repeats
    
    NUM_ACCESS=2
    
    from joblib import Parallel, delayed
    
    energies_list = Parallel(NUM_ACCESS)(
        delayed(calculate_energy_share)() for _ in range(NUM_ACCESS))
