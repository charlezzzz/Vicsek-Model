from RunHoomd import RunHoomd
from CellList import Cell, CellList
from CalculateGNF import CalculateGNF
import math
from math import ceil, sqrt, pow
import os
from os import system, path
import shutil
import pickle
import csv
import hoomd, hoomd.md
import gsd, gsd.hoomd
from hoomd.md import nlist as nl 
from InitializeUnitCell import InitializeUnitCell
import random


################################################################################################################


def erase_old_data(save_dir, density_plots_dir):
	if path.exists(save_dir):
		shutil.rmtree(save_dir)
	os.mkdir(save_dir)
	os.mkdir(save_dir +'/'+ density_plots_dir)

################################################################################################################

def run_simulation(L, rho, v0, delta, run_time, period, hoomd_output_path):

	hoomd.context.initialize()

	# Define the initial conditions of the new system
	N = int(rho*L*L)
	initialize = InitializeUnitCell(N, L)
	unit_cell = hoomd.lattice.unitcell(
	    N=N, 
	    a1=[L, 0, 0],
	    a2=[0, L, 0],
	    a3=[0, 0, 1],  # Set to [0,0,1] for 2D Simulations
	    dimensions=2,
	    position=initialize.randomized_positions(),
	    type_name=['A']*N,
	    orientation=initialize.randomized_orientations())
	system = hoomd.init.create_lattice(unit_cell, n=1)

	# Integrate the system forward in time
	seed = random.randint(1,10e4)
	nl_cell = nl.cell(r_buff=0)
	all = hoomd.group.all()
	hoomd.md.integrate.mode_standard(dt=1)
	hoomd.md.integrate.vicsek(all, nlist=nl_cell, r_cut=1, v0=v0, delta=delta, eta=0, seed=seed)
	d = hoomd.dump.gsd(hoomd_output_path, period=period, group=all, overwrite=True)
	hoomd.run(run_time)

################################################################################################################

def plot_current_positions(suffix, plot_path, system):
		import matplotlib
		import matplotlib.pyplot as plt
		import numpy as np
		from scipy.stats import gaussian_kde

		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)

		num_particles = len(system.particles.position)

		for i in range(0, num_particles):
			print("plotting particle: ", i, "/", num_particles, end="\r")
			position = system.particles.position[i]
			ax.scatter(position[0], position[1], color='k', marker='.', s=0.1)

		plt.savefig(plot_path + str(suffix) + '.png')
		plt.clf()


################################################################################################################

def main():

	L = 256
	rho = 0.5
	v0 = 0.5
	delta = 0.3
	
	save_dir = 'reproduce_chate_ginelli_2008_travelling_bands'
	plots_dir = 'plots'

	hoomd_output_name = 'trajectory.gsd'
	plot_name = 'system_t'
	

	#######################################################################################################

	hoomd_output_path = save_dir +'/'+ hoomd_output_name
	plot_path = save_dir +'/'+ plots_dir +'/'+ plot_name

	#######################################################################################################

	erase_old_data(save_dir, plots_dir)
	
	#######################################################################################################

	period = 5e3
	num_plots = 4
	run_time = num_plots*period

	run_simulation(L, rho, v0, delta, run_time, period, hoomd_output_path)
	snapshots = gsd.hoomd.open(name=hoomd_output_path, mode="rb")
	for t in range(0, len(snapshots)):
		plot_current_positions(period*t, plot_path, snapshots[t])

	#######################################################################################################


if __name__ == "__main__":
	system('clear')
	main()
