#!/bin/bash

rm -r Logs Data
python setup.py
sbatch runmeld.slurm
