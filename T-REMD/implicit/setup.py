import numpy as np
import meld
from meld import unit as u
from meld.system.scalers import LinearRamp
from meld.system.indexing import AtomIndex

n_replicas = 20
tau = 50 # in ps
n_exchanges = 1000
block_size=10

protein = meld.AmberSubSystemFromPdbFile("4AKE.amber.pdb")

build_options = meld.AmberOptions(
    forcefield="ff14sbside",
    implicit_solvent_model="gbNeck2",
    cutoff=1.8*u.nanometer
)

builder = meld.AmberSystemBuilder(build_options)
system = builder.build_system([protein]).finalize()

system.temperature_scaler = meld.GeometricTemperatureScaler(0, 0.8, 300. * u.kelvin, 400. * u.kelvin)

options = meld.RunOptions(
    timesteps = int(tau/0.002),
    minimize_steps = 20000,
 )

remd = meld.setup_replica_exchange(system, n_replicas=n_replicas, n_steps=n_exchanges)

meld.setup_data_store(system, options, remd, block_size)
