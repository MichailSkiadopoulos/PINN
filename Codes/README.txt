# Physics-Informed Neural Networks (PINN) for Ultrasonic Nondestructive Evaluation

This repository contains codes related to the development of Physics-Informed Neural Networks (PINNs) and data-driven Neural Networks (NNs) for ultrasonic nondestructive evaluation (NDE).

The workflow combines finite element simulations, synthetic ultrasonic waveform generation, and machine learning approaches for training, validation, testing, and transfer learning.

# Workflow Overview

The repository follows the general pipeline below:

1. Generate synthetic ultrasonic data using 3D finite element simulations in Abaqus
2. Submit and run simulations on a computing cluster
3. Extract simulated ultrasonic waveforms
4. Use the generated data to train, validate, and test:
   - Physics-Informed Neural Networks (PINNs)
   - Data-driven Neural Networks (NNs)
5. Investigate:
   - synthetic-only training
   - experimental-only training
   - transfer learning from synthetic to experimental data

# Finite Element Simulation and Data Generation

## `Generation_of_Abaqus_input_files.py`

This code generates 3D finite element models in **Abaqus** to simulate **pulse-echo ultrasonic testing**.

The generated simulations are used to create synthetic ultrasonic waveform datasets for the pretraining and development of the neural network models.

## `PBS_file_for_submission_in_ROAR_cluster.py`

This code generates/submits PBS job files for running the Abaqus simulations on the **ROAR computing cluster**.

It automates the large-scale execution of the generated finite element models.

## `Waveform_extraction_and_saving_totxt.py`

This code extracts the predicted ultrasonic waveforms from the Abaqus simulation outputs.

The extracted waveforms are saved and organized for later use in:

- training
- validation
- testing

of the neural network models.

# Neural Network Models

## Data-Driven Neural Networks

### `data_driven_NN_simulated_only.ipynb`

Training of a purely data-driven neural network using only synthetic/simulated ultrasonic data.

### `data_driven_NN_experimental_only.ipynb`

Training of a purely data-driven neural network using only experimental ultrasonic data.

### `data_driven_transfer_learning.ipynb`

Transfer learning framework where the data-driven neural network is:

1. pretrained using synthetic/simulated data
2. subsequently adapted to experimental data using partial retraining/fine-tuning

# Physics-Informed Neural Networks (PINNs)

## `physics_informed_NN_simulated_only.ipynb`

Training of a Physics-Informed Neural Network (PINN) using only synthetic/simulated data.

The training process incorporates physics-based constraints alongside data fitting.

## `physics_informed_NN_experimental_only.ipynb`

Training of a Physics-Informed Neural Network (PINN) using only experimental ultrasonic data.

## `physics_informed_transfer_learning.ipynb`

Transfer learning framework for the PINN approach where the model is:

1. pretrained using synthetic/simulated data
2. transferred and partially retrained using experimental data

This approach leverages physics-informed learning together with synthetic-to-experimental adaptation.