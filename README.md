# Vicsek-Model
Fall 2019 Research Rotation

## Equations of motion

TwoStepVicsek is an integration method for **HOOMD-blue** that implements the following equations of motion,

<img src="http://latex.codecogs.com/gif.latex?v_i^{t+\Delta{t}}=(v_0\Delta{t})\,e^{i\,\theta_i^{\,t}}" border="0"/>
<img src="http://latex.codecogs.com/gif.latex?\theta_i^{t+\Delta{t}}=\text{Arg}\Big[\sum_{j\in\mathcal{N}_{\text{i}}}e^{\theta_j^{\,t}}+\big|\mathcal{N}_{\text{i}}\big|\eta\,e^{\phi_i^{\,t}}\Big]+\delta\,\varphi_i^t" border="0"/>

where  <img src="http://latex.codecogs.com/gif.latex?\delta,\eta\in[0,1]" border="0"/>  are scalar and vectorial noise intensities, and <img src="http://latex.codecogs.com/gif.latex?\phi,\varphi\in[-\pi,\pi]" border="0"/>  are random variables uniformly distributed on the unit circle. The summation is carried over the set of neighbors of particle <img src="http://latex.codecogs.com/gif.latex?i" border="0"/>, <img src="http://latex.codecogs.com/gif.latex?\mathcal{N}_{\text{i}}" border="0"/>, which is the set of all particles <img src="http://latex.codecogs.com/gif.latex?j" border="0"/> within some cut-off radius <img src="http://latex.codecogs.com/gif.latex?r_{\text{cut}}" border="0"/> of particle <img src="http://latex.codecogs.com/gif.latex?i" border="0"/>, such that
 <img src="http://latex.codecogs.com/gif.latex?|\bold{r}_i-\bold{r}_j|<r_{\text{cut}}" border="0"/>. 

## Initializing systems

Systems obeying the above equations of motion typically have random initial positions and orientations for each particle. The InitializeUnitCell class can be used to create these initial conditions for a system of linear length 
<img src="http://latex.codecogs.com/gif.latex?a" border="0"/> 
containing
<img src="http://latex.codecogs.com/gif.latex?N" border="0"/> 
particles. An example of using this class to create a system with particle density 
<img src="http://latex.codecogs.com/gif.latex?\rho=1" border="0"/>
with periodic boundary conditions is:
```python
from InitializeUnitCell import InitializeUnitCell

...

N, a = 25, 5
initialize = InitializeUnitCell(N, a)
unit_cell = hoomd.lattice.unitcell(
	N=N, 
	a1=[a, 0, 0],
	a2=[0, a, 0],
	a3=[0, 0, 1],  # Set to [0,0,1] for 2D Simulations
	dimensions=2,
	position=initialize.randomized_positions(),
	type_name=['A']*N,
	mass=initialize.masses(),
	charge=initialize.charges(),
	diameter=initialize.diameters(),
	moment_inertia=initialize.inertias(),
	orientation=initialize.randomized_orientations())
system = hoomd.init.create_lattice(unit_cell, n=1)

...
```
