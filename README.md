# PINNs and DeepONet for Tumor Growth Modeling

This repository contains the source code accompanying the work **"Data-Driven Parameter Identification for Tumor Growth Models"** (Liu et al.).

The project applies **Physics-Informed Neural Networks (PINNs)** and **Deep Operator Networks (DeepONets)** to solve inverse problems arising in tumor growth modeling. The primary objective is to infer biologically meaningful parameters directly from observed tumor growth data while enforcing the underlying physical laws described by partial differential equations (PDEs).

Our framework combines data-driven learning with mathematical modeling to estimate tumor proliferation dynamics from limited and noisy observations. In addition to parameter identification, the repository includes implementations for learning solution operators of tumor growth PDEs using DeepONets.

The corresponding paper is available on arXiv:

**Data-Driven Parameter Identification for Tumor Growth Models**
http://arxiv.org/abs/2511.15940

## Main Features

* Physics-Informed Neural Networks (PINNs) for inverse PDE problems
* Deep Operator Networks (DeepONets) for operator learning
* Parameter identification in tumor growth models
* Robust learning from sparse and noisy data
* Applications to both synthetic and experimental tumor datasets
* PyTorch implementation

## Citation

If you find this repository useful in your research, please consider citing our paper.

