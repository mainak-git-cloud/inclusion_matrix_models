# Read Me

## Overview
This repository contains the necessary scripts for developing simulation files and analyzing experimental images. These scripts pertain to inclusion (cell)-matrix systems with multiple inclusions contracting and recovering, enabling the assessment of the mechanical properties of the matrix after inclusion recovery.

**Author:** Mainak Sarkar, University of Illinois at Urbana-Champaign (2024-2025)

---

## Simulations

### Step 1: Generating a Quasi-2D Fibrous Matrix
- Find the required codes in the repository: [fiber_network_model](https://github.com/jknotbohm/fiber_network_model).
- These scripts were developed by Mainak Sarkar and others in Prof. Jacob Notbohm’s research group at the University of Wisconsin-Madison.
- Navigate to the `1-NetworkGeneration` folder and locate the subfolder `simulated_annealing_2D`.
- Follow the associated README files.
- Run the scripts to generate and save the matrix model data in a `.MAT` file.

### Step 2: Introducing Inclusions
- Locate the `.MAT` file generated in Step 1.
- Load the file into the script `make_punctures_in_net.m`.
- Follow in-script comments to specify the file path.
- Customize the locations of circular holes by modifying lines #41-51.
- Run the script to generate a new `.MAT` file containing the updated model data with inclusions.

### Step 3: Generating the Abaqus Input File
- Run `boundary_condition_script.m`.
- Load the updated model definition.
- Execute the script to generate the Abaqus input (`.inp`) file.
- Modify solver settings using the function file `inp_gen_incl_ctr_plus_recv.m`.
- Sample `.inp` and `.odb` files are provided in the `example_inp_file` folder and at [Box folder](https://uofi.box.com/s/pex2d2nfivc2pl0ckpn02sfkhags1966), respectively.

---

## Experiments

### Spectral Energy Analysis of Confocal Microscopy Images
- Use the script `image_spectral_energy_finder.m` to compute the total spectral energy of a confocal microscopic image of the matrix.
- Add the path of the image to be analyzed in line #10 of the script.
- A representative image stack of the cell-free matrix is provided as an example.
