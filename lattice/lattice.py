from __future__ import print_function
from openmm import app
import openmm as mm
from openmm import unit
from sys import stdout
from openmmtools import testsystems
import mdtraj as md
import numpy as np

m=4
N=4*m**3
test_sys = testsystems.LennardJonesFluid(nparticles=N, lattice=True,reduced_density=0.5)
system, positions = test_sys.system, test_sys.positions


r=3.755

f=open('lattice.xyz','w')
f.write(str(N-1)+'\n ')
f.write(' '+'\n')


COM=np.zeros(3,float)

R=np.zeros((N,3),float)
for i in range(N):
    for a in range(3):
        COM[a]+=positions[i][a]._value
COM=(1.0/N)*COM

for i in range(N):
    for a in range(3):
        R[i][a]=positions[i][a]._value-COM[a]

index_vac=0

d_origin=10.0
for i in range(N):
    d=0.0
    for a in range(3):
        d+=R[i][a]**2
    d=np.sqrt(d)
    if d<d_origin:
        index_vac=i
        d_origin=d

dist=0.0
for a in range(3):
    dist+=(R[1][a]-R[0][a])**2
dist=np.sqrt(dist)

for i in range(N):
    if i!=int(index_vac):
        f.write('Ar ')
        for a in range(3):
            f.write(str(r*((R[i][a]-R[int(index_vac)][a]))/dist)+' ')
        f.write('\n')

#integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,1.0*unit.femtoseconds)

#platform = mm.Platform.getPlatformByName('Reference')

#simulation = app.Simulation(test_sys.topology, system, integrator, platform)
#simulation.context.setPositions(test_sys.positions)

#app.PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
#traj = md.load('output.pdb')
#traj.save_xyz('output.xyz')