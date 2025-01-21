# MELD Tutorials
This repo is a collection of minimal MELD tutorials for new and advanced users. I have implicit/explicit versions and also show how to run conventional MD, T-REMD and T,H-REMD all using MELD.

These tutorials were written to work with MELD v0.6.1 (the most recent version on conda). Note that you need access to as many GPUs as replicas you plan to run. For instance, if you want to run 20 replicas in one simulation, you will need access to 20 GPUs.

To install locally (and on some clusters) I typically have a local Anaconda installation and: `conda create -n meld -c conda-forge meld mpi4py`. That should do the trick. If you are installing on a cluster, you may need to `module load cuda` first or something like that.

If you run into any issues, please let me know and I will be happy to provide any assistance I can. I have experience also installing both OpenMM and MELD from source on a variety of HPC platforms.
