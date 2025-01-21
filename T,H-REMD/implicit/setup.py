import numpy as np
import meld
from meld import unit as u
from meld.system.scalers import LinearRamp
from meld.system.indexing import AtomIndex

kj = u.kilojoule/(u.mole*u.nanometer**2)
n_replicas = 20
tau = 50 # in ps
n_exchanges = 1000
block_size=10

protein = meld.AmberSubSystemFromPdbFile("4AKE.amber.pdb")

build_options = meld.AmberOptions(
    forcefield="ff14sbside",
    implicit_solvent_model="gbNeck2",
    cutoff=1.8 * u.nanometer
)
builder = meld.AmberSystemBuilder(build_options)

system = builder.build_system([protein]).finalize()
system.temperature_scaler = meld.GeometricTemperatureScaler(0, 0.8, 300. * u.kelvin, 375. * u.kelvin)

common_percent = 0.9
common_scaler = system.restraints.create_scaler('constant')
dists_common = []

# Collect MELD restraints into a list
with open("common.dat") as cfile:
    rest_group = []
    lines = cfile.read().splitlines()
    for line in lines:
        i1 = int(line.split()[0])
        i2 = int(line.split()[1])
        rest_group.append(system.restraints.create_restraint('distance', common_scaler, LinearRamp(0,100,0,1),
                                                             atom1=AtomIndex(i1), 
                                                             atom2=AtomIndex(i2),
                                                             r1=0.0*u.nanometer, r2=0.0*u.nanometer, 
                                                             r3=0.7*u.nanometer, r4=1.0*u.nanometer, 
                                                             k=500*kj))
    dists_common.append(system.restraints.create_restraint_group(rest_group, int(common_percent * len(rest_group))))

# Add MELD restraints to the system
system.restraints.add_selectively_active_collection(dists_common, 1)

options = meld.RunOptions(
    timesteps = int(tau/0.002),
    minimize_steps = 20000,
 )
remd = meld.setup_replica_exchange(system, n_replicas=n_replicas, n_steps=n_exchanges)

meld.setup_data_store(system, options, remd, block_size)
