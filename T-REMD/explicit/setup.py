import numpy as np
import meld
from meld import unit as u
from meld.system.scalers import LinearRamp
from meld.system.indexing import AtomIndex

n_replicas = 20
tau = 50 # in ps
n_exchanges = 100
block_size=10

protein = meld.AmberSubSystemFromPdbFile("4AKE.amber.pdb")

build_options = meld.AmberOptions(
    forcefield="ff15ipq",
    solvation="explicit",
    solvent_forcefield="spceb",
    solvent_distance=1.2,
    explicit_ions=True,
    p_ion="Na+",
    n_ion="Cl-",
    p_ioncount=0,
    n_ioncount=0,
    enable_pme=True,
    enable_pressure_coupling=True,
    use_big_timestep=False,
    cutoff=1.0*u.nanometer,
)

builder = meld.AmberSystemBuilder(build_options)
system = builder.build_system([protein]).finalize()

system.temperature_scaler = meld.ConstantTemperatureScaler(300*u.kelvin)
rest2_scaler = meld.ConstantTemperatureScaler(300*u.kelvin)
rest2_scaler = meld.GeometricTemperatureScaler(0, 0.8, 300. * u.kelvin, 400. * u.kelvin)

options = meld.RunOptions(
    timesteps = int(tau/0.002),
    minimize_steps = 20000,
    use_rest2=True,
    rest2_scaler=meld.REST2Scaler(300*u.kelvin, rest2_scaler),
 )

remd = meld.setup_replica_exchange(system, n_replicas=n_replicas, n_steps=n_exchanges)

meld.setup_data_store(system, options, remd, block_size)
