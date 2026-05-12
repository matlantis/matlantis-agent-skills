# Simulation of interfacial structures of Pt (111)/benzene

## Abstract
The main purpose of this example is to demonstrate how to build and validate the LightPFP model to simulate the solid-liquid interface. The PreFerred Potentential (PFP) is a powerful tool, however, its limited simulation size can pose challenges for complex systems like solid-liquid interfaces. To address this limitation, we developed LightPFP, which retains the accuracy of PFP while enabling much larger simulation boxes with faster computational speed for specific materials system. Using the Pt/benzene interface as a benchmark system, we illustrate the process of creating a LightPFP model and validate its performance by comparing the interfacial structure predicted by LightPFP and PFP models. Our results showcase the excellent agreement between the two models. At last, we demonstrate LightPFP's superior computational efficiency by runing molecular dynamics (MD) simulation with a large scale Pt/benzene interface.


## Introduction
Machine learning (ML) potentials have revolutionized atomic-level simulations by accurately capturing complex interatomic interactions with much faster computational speed than density functional theory (DFT). The Preferred Potential (PFP) is a versatile and universal ML potential. The version v5.0 is capable of handling any combination of 72 elements without additional training. However, its large model size makes it difficult to handle large simulation boxes due to GPU memory limitations. To overcome this, we developed LightPFP, a toolkit for training smaller ML potentials from PFP, offering faster computational speed and a focus on specific materials systems.

LightPFP utilizes the self-developed moment tensor potential (MTP) as the ML component[1]. MTP's simpler architecture and significantly smaller number of parameters compared to PFP enable knowledge distillation, where the larger PFP model acts as the teacher for training the smaller LightPFP (MTP) student model. Briefly, a set of training structures for the target material system is collected, and their energies, forces, and stresses are labeled using PFP as reference data. The LightPFP model is then trained by minimizing the difference in these properties between PFP and LightPFP. It's important to note that LightPFP is no longer a universal model; it is only valid for the specific material it is trained for.

In this study, we demonstrate LightPFP's capacity and reliability in simulating solid-liquid interfaces using the Pt/benzene interface as a model system. By leveraging LightPFP's enhanced speed and accuracy, researchers can gain insights into solid-liquid interface behavior and expand its applications in materials science, chemistry, and nanotechnology.


## Results and discussion
### Dataset 
To train the LightPFP model, we compiled a diverse set of structures with their energies, forces, and stresses labeled using the PFP. The training dataset encompassed various components: bulk Pt, Pt (111) slab, bulk benzene, and Pt (111)/benzene interface structures.

The dataset generation process was streamlined by a powerful tool within the LightPFP package, the dataset-generator. This tool facilitated the creation of multiple new structures from one or more initial structures, employing methods such as compressing, stretching, and deforming cell shapes, creating vacancies and surface structures, randomly displacing and substituting atoms, and crucially, conducting molecular dynamics (MD) simulations. The bulk Pt dataset was generated using all these methods, while the bulk benzene database was generated through cell compression, stretching, and MD simulations. The Pt (111) slab and Pt (111)/benzene interface datasets were generated solely using MD simulations. All calculations were performed using the PFP v5.0.0, with the PBE_PLUS_D3 calculation mode employed to capture Van der Waals interactions.

Typically, dataset collection is the most time-consuming step in building ML potentials, often taking several days to weeks. However, thanks to PFP's faster computational speed compared to DFT, we could swiftly create a robust training dataset within minutes to hours. It's worth noting that the dataset's accuracy is slightly worse than its DFT counterparts due to the inherent errors of the PFP model.

The dataset was randomly split into training and testing datasets, with 90% comprising the training dataset and the remaining 10% constituting the testing dataset. The script used to generate the training database can be found in the Appendix.

**Table I: Composition of the training dataset for LightPFP models**

| | Methods | Num of structures | Num of Atoms |
| - | - | -: | -: |
| Pt | Cell compression/Streching/Deforming; Atom Displacement; Vacancy; Surface and MD (NPT 500~1500K) | 2249 | 111,290 |
| Pt (111) | MD (NPT 500\~1000K) | 640 | 43,680 |
| Benzene | Cell compression/Streching; MD (NPT 300\~400K; NVT 500\~1500K) | 984 | 177,120 |
| Interface | MD (NPT 300\~800K; NVT 500\~1750K) | 5400 | 2,632,800 |


### Model training

We trained various LightPFP models by varying the hyperparameters, such as the cutoff radius, maximum radial basis function level (levmax), and μ cost. The detailed explaination of these parameters can be found in the original paper of MTP[1] and our guidebook. We trained six models with the same datasets and model settings. The energy, forces, and stress mean absolute errors (MAE) in the test dataset are listed in Table II. All models achieved high accuracy in the prediction of energy, forces, and stress.


**Table II: The accuracy of LightPFP models in the prediction of energy, force and stress in the test dataset.**

| model settings      | cutoff   | Energy MAE | Force MAE | Stress MAE |
| ------------------- | -------: | ---------: | --------: | ---------: |
|  levmax=8, μ cost=1 |      5.0 |       1.29 |      86.4 |      0.582 |
|  levmax=8, μ cost=1 |      6.0 |       1.37 |      90.8 |      0.489 |
|  levmax=8, μ cost=1 |      7.0 |       1.65 |      97.6 |      0.596 |
| levmax=16, μ cost=4 |      5.0 |       1.57 |      87.4 |      0.565 |
| levmax=16, μ cost=4 |      6.0 |       1.34 |      88.0 |      0.499 |
| levmax=16, μ cost=4 |      7.0 |       1.25 |      92.6 |      0.529 |



&lt;img src="assets/pt_benzene_interface/fig1.png"&gt;

**Fig 1: The prediction accruacy of the LightPFP model (levmax=8, μ cost=1, cutoff distance 6.0 Å)**


### Validation of LightPFP

To validate the performance of the LightPFP models, we performed MD simulations on a small Pt(111)/benzene interface system and compared the trajectories with PFP. The LightPFP model with levmax=8, μ cost=1, cutoff=6.0 Å was used for validation.

The initial simulation box has dimensions of 38.58 x 33.10 x 101.32 Å, with cell angles of 90&lt;sup&gt;°&lt;/sup&gt;, 90&lt;sup&gt;°&lt;/sup&gt;, and 120&lt;sup&gt;°&lt;/sup&gt;. It consists of a total of 9,072 atoms, including 1,512 Pt atoms and 630 benzene molecules. We have chosen a relatively smaller structure in order to compare the results with PFP.

The initial structure undergoes a 20 ps equilibrium at a temperature of 300.0 K at first, using the NVT ensemble. Subsequently, MD simulations are conducted for 100 ps at 300.0 K and 1 bar, utilizing the NPT ensemble. A snapshot is saved for future analysis.


#### Density

In Fig 2, the plotted density of the system can be observed. The density gradually decreases to approximately 5.08 g/cm&lt;sup&gt;3&lt;/sup&gt; after undergoing a 100 ps molecular dynamics (MD) simulation. The density of LightPFP aligns well with the density of PFP, which is 5.06 g/cm&lt;sup&gt;3&lt;/sup&gt;.

However, it is important to note that the accurate density should ideally be obtained from an even longer MD simulation as the system may not have fully equilibrated within the 100 ps timeframe. Since the primary purpose of this test is to validate the performance of LightPFP and considering the time required for PFP calculation, we have made the decision not to conduct a longer MD simulation.

&lt;img src="assets/pt_benzene_interface/fig2.png"&gt;

**Fig 2: Comparison of density between MTP and PFP MD trajectories.**


#### Radial distribution function

To calculate the radial distribution function, we used the snapshots of the MD trajectory taken after 50 ps. The resulting radial distribution function is shown in Fig 3. We discovered that the results obtained from both the MTP and PFP trajectories exhibit a high degree of agreement.

&lt;img src="assets/pt_benzene_interface/rdf.png"&gt;

**Fig 3. Comparison of the radial distribution function between the MTP and PFP trajectories.**


#### Distribution of atomic positions

To analyze the spatial distribution of atoms along the Z direction (perpendicular to the interface), we captured snapshots of the MD trajectory after 50 ps, and statistic the Z positions of each atom by applying a Gaussian function with a width of 0.25 Å. The resulting spatial distribution of atoms is plotted in Fig4.

&lt;table&gt;
  &lt;tr&gt;
    &lt;td&gt;&lt;img src="assets/pt_benzene_interface/z_distribution.png"&gt;&lt;/td&gt;
    &lt;td&gt;&lt;img src="assets/pt_benzene_interface/z_distribution_zoom.png"&gt;&lt;/td&gt;
  &lt;/tr&gt;
&lt;/table&gt;

**Fig 4: The spatial distribution of Pt, C and O atoms along the Z direction. (a) Whole simulation box. (b) Zoom in the interface.**


The curve in Fig4 shows several peaks in the H and C elements near the Pt surface (around 20 Å), indicating that the benzene structure is significantly different from the uniform liquid phase due to adsorption with the Pt surface. As we move further away from the Pt surface, the interfacial liquid structure gradually transitions towards a uniform liquid phase. According to the figure, the thickness of the interfacial layer is approximately 15 Å.

The position and intensity of the peaks in the LightPFP model were compared with the PFP result, and they matched each other across most regions. This alignment demonstrates that the LightPFP model effectively captures and represents the structural characteristics of the solid/liquid interface.

Symmetric peaks can be observed around z=100 Å, indicating the presence of the same liquid-solid interface due to the periodic boundary condition. However, these peaks exhibit slight fluctuations in shape due to noise originating from cell shape changes in the NPT MD simulation.


#### Benzene adsorption on Pt surface

&lt;table&gt;
  &lt;tr&gt;
    &lt;td&gt;&lt;img src="assets/pt_benzene_interface/benzene_adsorption_light_pfp.png" alt="Image 1" width="300" height="200"&gt;&lt;/td&gt;
    &lt;td&gt;&lt;img src="assets/pt_benzene_interface/benzene_adsorption_pfp.png" alt="Image 2" width="300" height="200"&gt;&lt;/td&gt;
    &lt;td&gt;&lt;img src="assets/pt_benzene_interface/benzene_adsorption_single.png" alt="Image 3" width="300" height="200"&gt;&lt;/td&gt;
  &lt;/tr&gt;
&lt;/table&gt;

**Fig 5: The benzene adsorption on Pt surface (a) LightPFP (b) PFP. (c) The bri30 adsorption site of the benzene molecule.**


In Fig 5, we can observe the benzene molecules that have adsorbed onto the Pt surface. Specifically, in the LightPFP MD simulation, 18 benzene molecules were found to be adsorbed onto the Pt surface within the specified area. In the PFP MD simulation, on the other hand, there were a total of 21 benzene molecules observed to be adsorbed onto the Pt surface within the same area. In conclusion, the LightPFP MD simulation reproduced the surface coverage rate of benzene on the Pt surface well.

Figure 5(c) illustrates the adsorption structure of a single benzene molecule. The bri30 conformation, in which the center of the benzene molecule is located on the bridge site of the Pt surface, was found to be the most energetically stable in first-principles calculations [2]. Interestingly, we observed that almost all the molecules were adsorbed in the bri30 site in both the MTP and PFP simulations. This result is consistent with the findings of the first-principle calculations.


### Large scale simulation

We utilized a large simulation box containing 84,600 atoms and employed the LightPFP model. The system initially underwent equilibrium for 20 ps at 300.0 K under the NVT ensemble. Subsequently, MD simulations were performed for 100 ps at 300.0 K and 1 bar using the NPT ensemble. The final structure of the system can be observed in Fig 6.

&lt;img src="assets/pt_benzene_interface/Pt_benzene.png"&gt;

**Fig 6 The snapshot of the Pt/benzene interface (84,600 atoms) after 100 ps molecular dynamics simulation.**

The entire MD simulation, lasting a total of 120 ps, took approximately 9.2 hours to complete on the Matlantis platform. The MD speed was calculated to be about 3.19 μs/atom/step. In comparison, we estimated the speed of PFP in the MD simulation, based on the linear extrapolation of the small-scale MD simulation speed, to be around 147 μs/atom/step. Thus, the LightPFP demonstrated a speed advantage of approximately 46 times compared to PFP for inference in MD simulations.

Please note that when using LightPFP, additional time is required for dataset collection and model training. Based on rough estimates, the time cost for dataset collection is approximately 6.0 hours, while model training requires an additional 1 hour. The overall time cost for PFP and LightPFP simulations, for a system size of about 10000 atoms, can be found in the table III. Generally, the LightPFP does not exhibit a significant advantage for short MD simulations (e.g., less than 10 ps). However, its speed advantage becomes apparent when running MD simulations for over 100 ps, even when considering the time cost of dataset collection and model training.

Additionally, the upper limitation of the simulation box is significantly larger when using LightPFP. Taking the model parameters (levmax=8, μ cost=1, cutoff distance of 6.0 Å) as an example, the upper limitation of the simulation box is estimated to be around 0.5 million atoms. This is much larger than the limitation of PFP, which is around 10000 atoms. Therefore, users are encouraged to utilize LightPFP when simulating large systems.


**Table III: The comparison of overall timecost (in hours) between PFP and LightPFP. (Number of atoms: 10,000)**

| MD steps   | PFP      | LightPFP |
| ---------: | -------: | -------: |
|      1,000 |     0.41 |     7.01 |
|     10,000 |     4.08 |     7.09 |
|    100,000 |    40.83 |     7.89 |
|  1,000,000 |   408.33 |    15.86 |


The present example focuses on demonstrating the LightPFP model for simulating the clean Pt/benzene interface, a system that can also be modeled using the original PFP method. However, many research topics in solid-liquid interfaces involve more complex scenarios that necessitate much larger simulation scales and longer time scales, which can be challenging for the conventional PFP approach.

One such area of interest is the effect of surface defects, such as adatoms, stacking faults, and surface steps, on molecular adsorption and interfacial properties. Capturing the influence of these defects requires simulating a sufficiently large surface area, which can quickly become computationally intractable for PFP as the system size increases.

Additionally, investigating dynamical processes like atom/molecule diffusion, adsorption, and desorption events demands simulations over extended time scales to adequately sample these rare events. While PFP simulations can access relatively short time scales, probing longer time scales becomes increasingly difficult due to the computational cost.

In these research areas, the LightPFP model can be a powerful tool, enabling large-scale simulations of complex solid-liquid interfaces with surface defects and the study of rare dynamical events over long time scales. By retaining the accuracy of PFP while offering superior computational efficiency, LightPFP can open up new avenues for exploring intricate phenomena at solid-liquid interfaces, ultimately advancing our understanding of these ubiquitous systems.


## Conclusion

In conclusion, the LightPFP model has been successfully trained and validated for simulating the Pt/benzene interface. It shows excellent accuracy and reliability, capturing the complex interatomic interactions at the solid-liquid interface. The density, radial distribution function, and atomic position distribution agree well with the PFP results. Moreover, the LightPFP model offers faster computational speed and can handle large-scale simulations, making it a valuable tool for studying solid-liquid interfaces in materials science


## Reference

[1] Novikov, I. S., Gubaev, K., Podryabinkin, E. V, &amp; Shapeev, A. V. (2021). The MLIP package: moment tensor potentials with MPI and active learning. Machine Learning: Science and Technology, 2, 025002. https://doi.org/10.1088/2632-2153/abc9fe

[2] Morin, C., Simon, D., &amp; Sautet, P. (2004). Chemisorption of benzene on Pt(111), Pd(111), and Rh(111) metal surfaces: A structural and vibrational comparison from first principles. Journal of Physical Chemistry B, 108(18), 5653–5665. https://doi.org/10.1021/jp0373503




# Script and Code
# (Please run the following command on CPU notebook)

## Dataset Generation

Collect the training dataset with different types of structures. 
It will take about 6.0 hours to finish all these tasks.
Note that `assets/pt_benzene_interface_datasets/Pt.cif` (mp-126) used in `assets/pt_benzene_interface_datasets/bulk.json` is downloaded from the Materials Project under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)


```
dataset-generation -c assets/pt_benzene_interface_datasets/bulk.json
dataset-generation -c assets/pt_benzene_interface_datasets/surface.json
dataset-generation -c assets/pt_benzene_interface_datasets/liquid_small.json
dataset-generation -c assets/pt_benzene_interface_datasets/liquid_large.json
dataset-generation -c assets/pt_benzene_interface_datasets/interface_small.json
dataset-generation -c assets/pt_benzene_interface_datasets/interface_large.json
dataset-generation -c assets/pt_benzene_interface_datasets/interface_short_small.json
dataset-generation -c assets/pt_benzene_interface_datasets/interface_short_middle.json
dataset-generation -c assets/pt_benzene_interface_datasets/interface_short_large.json
```

Validate the datasets and check whether there is some abnormal structures which forces larger than 30 eV/Å.
(You can also use the dataset saved in "./assets".)

```
validate-dataset -d *.h5 --max_forces 30.0
```


## Submit training job

Run a training job with the dataset collected from the last step. 
The dataset

```
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-small-6.json --model-name pt-benzene-interface-small-6
```

The user can also try different model settings. However, it will cost additional computational resourch.

```
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-large-5.json --model-name pt-benzene-interface-large-5
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-large-6.json --model-name pt-benzene-interface-large-6
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-large-7.json --model-name pt-benzene-interface-large-7
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-small-5.json --model-name pt-benzene-interface-small-5
submit-light-pfp-training -d *.h5 -c assets/pt_benzene_interface/job-small-7.json --model-name pt-benzene-interface-small-7
```



## Preparing packages
* Install the necessary packages.


```python
# If you have already installed these packages, you may skip this step.
%pip install -U light-pfp-client light-pfp-data
```



## MD simulation

### Run MD simulation with PFP

The MD simulation: 
1. Run structural optimization.
2. Run MD simulation with NVT ensemble to equilibrium the system.
3. Run MD simulation with NPT ensemble and save the trajectory.


```python
import os
import json
import numpy as np
from time import perf_counter

from ase import units
from ase.io import Trajectory
from ase.optimize import FIRE
from ase.filters import FrechetCellFilter
from ase.md.npt import NPT
from ase.md.langevin import Langevin

class PrintDyn:
    def __init__(self, dyn, logfile=None):
        self.dyn = dyn
        self.st = perf_counter()
        self.logfile = logfile
        if self.logfile is not None:
            with open(self.logfile, "w") as fd:
                fd.write("# step E_tot E_pot density T elapsed_time\n")
    def __call__(self):
        dyn = self.dyn
        msg = (
            f"{dyn.get_number_of_steps(): 6d} {atoms.get_total_energy():.3f} {atoms.get_potential_energy():.3f} "
            f"{atoms.get_masses().sum() / units.kg / atoms.get_volume() * 1e27:.5f} "
            f"{atoms.get_temperature():.1f} {perf_counter() - self.st:.2f}"
        )
        print(msg)
        if self.logfile is not None:
            with open(self.logfile, "a") as fd:
                fd.write(msg+"\n")


def run_md(atoms, save_dir, temperature, init_steps, steps):
    print(f"Number of atoms: {len(atoms)}")

    os.makedirs(f"{save_dir}", exist_ok=True)

    fire = FIRE(atoms)
    fire.run(fmax=0.1, steps=100) # set max_step = 100
    fcf = FrechetCellFilter(atoms, hydrostatic_strain=True)
    fire = FIRE(fcf)
    fire.run(fmax=0.1, steps=1000) # set max_step = 1000
    atoms.write(f"{save_dir}/opt.cif")

    md = Langevin(
        atoms, 
        timestep=units.fs, 
        temperature_K=300.0,
        friction=0.1
    )
    print_dyn=PrintDyn(md)
    md.attach(print_dyn, interval=100)
    md.run(steps=init_steps)
    atoms.write(f"{save_dir}/start.cif")

    traj = Trajectory(f"{save_dir}/md.traj", "w", atoms = atoms)
    md = NPT(
        atoms, 
        timestep=units.fs, 
        temperature_K=temperature,
        externalstress=101325 * units.Pascal,
        mask = np.eye(3),
        ttime = 20.0 * units.fs,
        pfactor = 2e6 * units.GPa*(units.fs**2),
        logfile =None,
    )

    print_dyn=PrintDyn(md, logfile=f"{save_dir}/md.log")
    md.attach(traj.write, interval=10000)
    md.attach(print_dyn, interval=100)
    md.run(steps=steps)
    atoms.write(f"{save_dir}/final.cif")

```



Run with the PFP calculator.
It will cost about 45 hours to finish. 

The run command is shown for reference purposes only.
To run it succesfully, a professional plan is required due to the high number of atoms/neighbors.

The user can also use the result saved in `./assets/pt_benzene_interface_pfp_md`.


```python
from ase.io import read

from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator(model_version="v5.0.0", calc_mode=EstimatorCalcMode.PBE_PLUS_D3)
calc = ASECalculator(estimator)
```


```python
atoms = read("assets/pt_benzene_interface/interface_init.xyz") * (7, 6, 1)
atoms.calc = calc

#run_md(atoms, "pfp_md", 300.0, 20000, 100000)
```



### Run MD simulation with LightPFP
# (Please run the following script on LightPFP notebook)

Run MD simulation again with the LightPFP model. It will takes about 2 hours to finish.


```python
from ase.io import read

from light_pfp_client import Estimator, ASECalculator

# Please replace it with the model_id of the model you trained
model_id = "examplesv1ptbnzn"

estimator = Estimator(model_id=model_id, book_keeping=True)
calc = ASECalculator(estimator)
```


```python
atoms = read("assets/pt_benzene_interface/interface_init.xyz") * (7, 6, 1)
atoms.calc = calc

run_md(atoms, "light_pfp_md", 300.0, 20000, 100000)
```



### LightPFP validation
#### Density

Compare the density between the PFP trajectory and LightPFP trajectory.


```python
import numpy as np
import matplotlib.pyplot as plt

pfp_log = np.loadtxt("assets/pt_benzene_interface_pfp_md/md.log")
light_pfp_log = np.loadtxt("light_pfp_md/md.log")
 
fig, ax = plt.subplots(1,1)
ax.plot(pfp_log[:,0]/1000, pfp_log[:,3], label="PFP")
ax.plot(light_pfp_log[:,0]/1000, light_pfp_log[:,3], label="LightPFP")
plt.legend()
plt.xlabel("time (ps)")
plt.ylabel("density (g/cm$^3$)")
```



#### Radial distribution function


```python
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm
from ase.io import Trajectory


def get_distance_matrix(atoms):
    structure = AseAtomsAdaptor.get_structure(atoms)
    return structure.distance_matrix
    

def get_rdf(traj, cutoff, bin_width):
    steps = len(traj)
    elements = traj[0].get_atomic_numbers()
    elements_unique = np.unique(elements)
    n_elements = len(elements_unique)
    rdf_bins = np.append(np.arange(0, cutoff, bin_width), cutoff)
    bin_center = (rdf_bins[1:] + rdf_bins[:-1]) / 2.0
    rdf_values = np.zeros([n_elements, n_elements, len(rdf_bins)-1])
    for atoms in tqdm(traj):
        distance = get_distance_matrix(atoms)
        volume = atoms.get_volume()
        for j in range(n_elements):
            for k in range(n_elements):
                indices_1 = elements == elements_unique[j]
                indices_2 = elements == elements_unique[k]
                rho = np.sum(indices_2) / volume
                if j == k:
                    triu_indices = np.triu_indices(np.sum(indices_1), k=1)
                    dist = distance[indices_1][:, indices_2][triu_indices]
                    bin_counts= np.histogram(dist, rdf_bins)[0] * 2 / np.sum(indices_1)
                else:
                    bin_counts= np.histogram(distance[indices_1][:, indices_2], rdf_bins)[0]  / np.sum(indices_1)
                volume_shell = 4.0 * np.pi * bin_center**2 * bin_width
                rdf_values[j, k] += bin_counts / volume_shell / rho 
    rdf_values /= steps
    return bin_center, rdf_values


pfp_traj = Trajectory("assets/pt_benzene_interface_pfp_md/md.traj")
light_pfp_traj = Trajectory("light_pfp_md/md.traj")

bin_center, pfp_rdf = get_rdf(pfp_traj[5:], 8.0, 0.01)
bin_center, light_pfp_rdf = get_rdf(light_pfp_traj[5:], 8.0, 0.01)
```


```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 3, figsize=[15, 10])

def ax_plot(ax, pfp_rdf, light_pfp_rdf, i, j):
    ax.plot(bin_center, pfp_rdf[i, j], label="PFP")
    ax.plot(bin_center, light_pfp_rdf[i, j], label="LightPFP")
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("g(r)")
    ax.legend()


ax_plot(ax[0, 0], pfp_rdf, light_pfp_rdf, 0, 0)
ax[0,0].set_title("H-H")
ax_plot(ax[0, 1], pfp_rdf, light_pfp_rdf, 0, 1)
ax[0,1].set_title("H-C")
ax_plot(ax[0, 2], pfp_rdf, light_pfp_rdf, 0, 2)
ax[0,2].set_title("H-Pt")
ax_plot(ax[1, 0], pfp_rdf, light_pfp_rdf, 1, 1)
ax[1,0].set_title("C-C")
ax_plot(ax[1, 1], pfp_rdf, light_pfp_rdf, 1, 2)
ax[1,1].set_title("C-Pt")
ax_plot(ax[1, 2], pfp_rdf, light_pfp_rdf, 2, 2)
ax[1,2].set_title("Pt-Pt")
```



#### Distribution of atomic positions


```python
import scipy
import numpy as np
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from ase.io import Trajectory


x = np.arange(-5.0, 110.0, 0.1)
z_cliab = 10.0

def gaussian_smoothing(data, bandwidth, x):
    y = np.zeros_like(x)
    for d in tqdm(data, leave=False):
        kernel = scipy.stats.norm(d, bandwidth).pdf(x)
        y += kernel
    y /= len(data)
    return y

def analyze_atoms(atoms, numbers=[1, 6, 78]):
    z = atoms.positions[:,2]
    ind  = (atoms.numbers == 78 ) &amp; (z &gt; 80.0)
    z[ind] -= atoms.cell.array[2,2]
    z_shift = z[atoms.numbers==78].mean() - z_cliab
    z -= z_shift
    y = np.empty([len(x), len(numbers)])
    for i, n in tqdm(enumerate(numbers)):
        y[:, i] = gaussian_smoothing(z[atoms.numbers==n], 0.25, x)
    return y

pfp_traj = Trajectory("assets/pt_benzene_interface_pfp_md/md.traj")
light_pfp_traj = Trajectory("light_pfp_md/md.traj")

n_atoms_pfp = len(pfp_traj[0])
n_atoms_light_pfp = len(light_pfp_traj[0])

distribution_pfp = np.mean([analyze_atoms(atoms) for atoms in pfp_traj[5:]], axis=0)
distribution_light_pfp = np.mean([analyze_atoms(atoms) for atoms in light_pfp_traj[5:]], axis=0)

plt.figure(figsize=[15, 8])
plt.plot(x, distribution_pfp[:,0], c="b", label="H PFP")
plt.plot(x, distribution_pfp[:,1], c="r", label="C PFP")
plt.plot(x, distribution_pfp[:,2], c="g", label="Pt PFP")

plt.plot(x, distribution_light_pfp[:,0], c="b", ls="--", label="H MTP")
plt.plot(x, distribution_light_pfp[:,1], c="r", ls="--", label="C MTP")
plt.plot(x, distribution_light_pfp[:,2], c="g", ls="--", label="Pt MTP")

plt.xlabel("z position (A)")
plt.legend()
```



#### Benzene adsorption on Pt surface


```python
from ase.io import Trajectory


def get_interface_1_layer(atoms):
    atoms_benzene = atoms[atoms.numbers!= 78]
    atoms_pt = atoms[(atoms.numbers == 78) &amp; (atoms.positions[:,2] &gt; 17.0) &amp; (atoms.positions[:,2] &lt; 50.0)]
    n_mol = len(atoms_benzene) // 12
    for i in range(n_mol):
        mol = atoms_benzene[i*12:(i+1)*12]
        com = mol.get_center_of_mass()
        if com[2] &lt; 21:
            atoms_pt += mol
    return atoms_pt


pfp_traj = Trajectory("assets/pt_benzene_interface_pfp_md/md.traj")
light_pfp_traj = Trajectory("light_pfp_md/md.traj")

atoms_ad_1_pfp = get_interface_1_layer(pfp_traj[-1])
atoms_ad_1_pfp.write("adsorption_1_pfp.cif")

atoms_ad_1_mtp = get_interface_1_layer(light_pfp_traj[-1])
atoms_ad_1_mtp.write("adsorption_1_light_pfp.cif")
```



### Large scale simulation


```python
from ase.io import read

from light_pfp_client import Estimator, ASECalculator

# Please replace it with the model_id of the model you trained if you want for reproduction purposes.
# By default, a pretrained model is used that can be found in Model Library -&gt; Example tab.
model_id = "examplesv1ptbnzn"

estimator = Estimator(model_id=model_id, book_keeping=True)
calc = ASECalculator(estimator)
```


```python
atoms = read("assets/pt_benzene_interface/interface_init.xyz") * (20, 20, 1)
atoms.calc = calc
print(f"Number of atoms: {len(atoms)}")

run_md(atoms, "light_pfp_md_large", 300.0, 20000, 100000)
```
