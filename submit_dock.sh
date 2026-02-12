#!/bin/bash
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --job-name=anchored_dock
#SBATCH --output=dock_%j.log

module load pixi/0.56.0
module load AutoDock-GPU/1.5.3-CUDA
pixi run dock

