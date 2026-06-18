# Matlantis Agent Skills

[Japanese version / 日本語版](README.ja.md)

A plugin marketplace providing Matlantis / PFP atomistic simulation skills for Claude Code and Codex.

## Overview

Matlantis Agent Skills is a collection of skills that enables Claude Code and Codex to assist with atomistic simulations on the [Matlantis](https://matlantis.com/) platform. 16 skill modules guide Python code generation for structure optimization, molecular dynamics, property calculations, and more — all powered by PFP (Preferred Potential) neural network potentials.

Each skill systematically covers implementation patterns, best practices, and common errors with solutions, all built on ASE (Atomic Simulation Environment).

## Skills Catalog

| Category | Skill | Summary |
|----------|-------|---------|
| **Foundation** | [matlantis](#matlantis) | Platform overview, PFP calculation mode selection, library architecture |
| **Foundation** | [setup](#setup) | Calculator initialization, calc_mode configuration, structure file I/O |
| **Foundation** | [ase-basics](#ase-basics) | ASE fundamentals, Atoms object, molecule & crystal construction |
| **Structure Building** | [modeling](#modeling) | Surface slab generation, SMILES→3D conversion, substitution/doping, liquids & amorphous |
| **Structure Building** | [external-data](#external-data) | Structure retrieval from PubChem, Materials Project, etc.; file format conversion |
| **Computation** | [optimization](#optimization) | Fixed-cell & variable-cell relaxation, symmetry preservation, multi-stage optimization |
| **Computation** | [dynamics](#dynamics) | NVE/NVT/NPT molecular dynamics, metadynamics, transport properties |
| **Computation** | [reaction](#reaction) | NEB reaction path search, transition state search, activation energy |
| **Computation** | [adsorption](#adsorption) | Adsorption structure construction, adsorption energy, adsorption site search |
| **Computation** | [property-analysis](#property-analysis) | Elastic constants, phonons, QHA, diffusion coefficients, viscosity, thermal conductivity |
| **Workflow** | [integrated-workflows](#integrated-workflows) | Structure screening, adsorption campaigns, defect analysis pipelines |
| **Workflow** | [background-job](#background-job) | Background execution of long-running calculations, checkpoint management |
| **Acceleration** | [light-pfp](#light-pfp) | Light-PFP custom models (15–40x faster) training & inference |
| **Utility** | [pfcc-extras](#pfcc-extras) | Visualization, structure processing, and job execution control utilities |
| **Utility** | [visualization](#visualization) | Atomic structure & trajectory visualization, nglview, POV-Ray |
| **Infrastructure** | [storage](#storage) | Group Drive persistent storage, team data sharing |
| **Infrastructure** | [ssh](#ssh) | SSH connection setup, VSCode Remote-SSH, WinSCP integration |

## Getting Started

### Prerequisites

- **Claude Code** CLI installed when using Claude Code
  - A version supporting the plugin feature is required
  - Installation: https://claude.ai/code
- **Codex** CLI installed when using Codex plugins
  - Installation: https://developers.openai.com/codex/cli

### Installation

#### Claude Code

Run the following commands in the Claude Code REPL:

```
# 1. Register the marketplace
/plugin marketplace add https://github.com/matlantis/matlantis-agent-skills
```

```
# 2. Open the plugin manager
/plugin

# 3. Enable the matlantis plugin from the "Marketplaces" tab
```

#### Codex

Run the following commands from your shell:

```bash
# 1. Register the marketplace
codex plugin marketplace add matlantis/matlantis-agent-skills

# 2. Install the plugin from the marketplace
codex plugin add matlantis@matlantis-agent-skills
```

Once installed or enabled, the agent will automatically reference the appropriate skills when you work with Matlantis-related code or questions.

## Skills Details and Usage Examples

### matlantis

Platform overview skill. Covers PFP calculation mode selection guidelines (PBE, R2SCAN, WB97XD, etc.), library hierarchy, and critical constraints.

```
"I want to calculate the surface energy of copper. Which calc_mode should I use?"
```

### setup

Calculator initialization and environment setup. Covers Estimator/ASECalculator configuration patterns, structure file I/O, and batch execution settings.

```
"Write code to initialize the PFP Calculator in R2SCAN mode"
```

### ase-basics

ASE fundamentals. Covers Atoms object creation/editing, molecule & crystal construction, periodic boundary conditions, and file I/O.

```
"Load a crystal structure from a CIF file and create a 2x2x2 supercell"
```

### modeling

Structure modeling. Covers surface slab generation (ASE builders / pymatgen), 3D molecule generation from SMILES (RDKit), atomic substitution/doping, and construction of liquids, amorphous materials, and high-entropy alloys.

```
"Create a 4-layer FCC copper (111) surface slab with 15 Å of vacuum"
```

### external-data

Structure retrieval from external databases. Covers integration with PubChem (organic molecules), Materials Project (inorganic crystals), AFLOW, COD, and file format conversion.

```
"Retrieve the BaTiO3 crystal structure from Materials Project and optimize it"
```

### optimization

Structure optimization (relaxation). Covers fixed-cell optimization (atomic positions only) and variable-cell optimization (including lattice vectors), LBFGS/BFGS/FIRE optimizers, and symmetry-preserving optimization.

```
"Relax the crystal structure with variable-cell optimization, preserving symmetry"
```

### dynamics

Molecular dynamics (MD) simulation. Covers NVE/NVT/NPT ensembles, Langevin thermostat, PLUMED metadynamics, trajectory analysis, and transport property calculations (diffusion, viscosity, thermal conductivity).

```
"Run a 100 ps NVT MD simulation at 300 K"
```

### reaction

Reaction path and transition state search. Covers NEB (Nudged Elastic Band), CI-NEB, and reaction barrier calculations using ReactionStringFeature.

```
"Calculate the activation energy for CO dissociation on a surface using NEB"
```

### adsorption

Molecular adsorption on surfaces. Covers adsorption structure construction with add_adsorbate, adsorption energy calculation with D3 dispersion correction, and adsorption site search (Optuna black-box optimization recommended).

```
"Find the most stable adsorption site for CO on Pt(111) and compute the adsorption energy"
```

### property-analysis

Property calculations and MD post-processing. Covers elastic constants/moduli, phonon band structures & density of states, quasi-harmonic approximation (QHA), diffusion coefficients, viscosity, thermal conductivity, and radial distribution functions (RDF).

```
"Calculate the phonon band structure and density of states for silicon"
```

### integrated-workflows

Integrated workflow (pipeline) templates. Provides structure screening (multi-structure batch optimization + energy ranking), adsorption campaigns (site search + ranking), and defect analysis pipelines.

```
"Batch-screen 10 alloy compositions and rank them by stability"
```

### background-job

Background execution of long-running calculations. Covers detached computation from notebook sessions, checkpoint/restart, and job monitoring.

```
"Run a large-scale MD in the background and monitor progress"
```

### light-pfp

Light-PFP custom models. Covers the workflow for training MTP (Moment Tensor Potential) models specialized for specific materials, achieving 15–40x faster inference than PFP.

```
"Create a Light-PFP model specialized for copper-based materials"
```

### pfcc-extras

pfcc_extras utility toolkit. Covers the visualization layer (show_gui, view_ngl), structure processing layer (partial occupancy, molecule wrapping), MD event scheduling, and job execution control.

```
"Run multiple structure optimization jobs in parallel with resource limits"
```

### visualization

Atomic structure & trajectory visualization. Covers simple display with pfcc_extras.show_gui, interactive visualization with nglview, trajectory animation, and high-quality POV-Ray rendering.

```
"Display an MD trajectory as an animation and render specific frames in high quality"
```

### storage

Group Drive persistent storage. Covers file upload/download to cloud object storage and team data sharing.

```
"Save calculation results to Group Drive and share with team members"
```

### ssh

SSH connection setup. Covers WebSocket-based SSH tunneling, VSCode Remote-SSH integration, and file transfer via WinSCP.

```
"Connect to the Matlantis environment from VSCode via Remote-SSH"
```

## License

Copyright 2026 Matlantis Corporation.

This project is licensed under the [Apache License 2.0](LICENSE).

## Third-Party Licenses

This repository contains documentation from Matlantis, which is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0).
Copyright 2021-present, Preferred Networks, Inc. and ENEOS Corporation.
