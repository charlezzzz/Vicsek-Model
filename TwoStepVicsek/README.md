## Installation

In order to use the TwoStepVicsek integration method
1) install the **HOOMD-blue** base code from https://github.com/glotzerlab/hoomd-blue
2) cd into the hoomd/md directory
3) copy/overwrite the files in this folder into the md directory 

## Job scripts

A simple example of using TwoStepVicsek with Hoomd's neighbor list:

```python
import hoomd
from hoomd import md
from hoomd.md import nlist as nl 

hoomd.context.initialize()
unit_cell = hoomd.lattice.unitcell(
	N=1,
	a1=[2, 0, 0],
	a2=[0, 2, 0],
	a3=[0, 0, 1],  # Set to [0,0,1] for 2D Simulations
	dimensions=2,
	position=[[0,0,0]],
	type_name=['A'],
	orientation=[[0.678,0,0,1]]
	)
system = hoomd.init.create_lattice(unit_cell, n=5)
nl_cell = nl.cell(r_buff=0)
all = hoomd.group.all();
hoomd.md.integrate.mode_standard(dt=1);
hoomd.md.integrate.vicsek(all, nlist=nl_cell, r_cut=0.5, v0=1, eta=0.2, seed=12)
hoomd.run(1e2)
```
