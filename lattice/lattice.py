from __future__ import print_function
from openmm import app
import openmm as mm
from openmm import unit
from sys import stdout
from openmmtools import testsystems
import mdtraj as md
import numpy as np

N=108
m=3
N=4*m**3
test_sys = testsystems.LennardJonesFluid(nparticles=108, lattice=True,reduced_density=0.5)
system, positions = test_sys.system, test_sys.positions


r=3.72

f=open('lattice.xyz','w')
f.write(str(N-1)+'\n ')
f.write(' '+'\n')

for i in range(N):
    if i!=int((N-1)/2):
        f.write('Ar ')
        for a in range(3):
            f.write(str(r*((positions[i][a]._value-positions[int((N-1)/2)][a]._value))/positions[1][1]._value/np.sqrt(2))+' ')
        f.write('\n')

#integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,1.0*unit.femtoseconds)

#platform = mm.Platform.getPlatformByName('Reference')

#simulation = app.Simulation(test_sys.topology, system, integrator, platform)
#simulation.context.setPositions(test_sys.positions)

#app.PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
#traj = md.load('output.pdb')
#traj.save_xyz('output.xyz')