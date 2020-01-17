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

class TrajectoryDetails:
	
	def __init__(self, simulation_details):
		self.time_integrated = 1
		self.L = None
		self.rho = None
		self.v0 = None
		self.delta = None
		self.path = simulation_details

################################################################################################################


def erase_old_data(save_dir, states_dir, posplots_dir, gnf_data_dir, cell_list_dir):
	if path.exists(save_dir):
		shutil.rmtree(save_dir)
	os.mkdir(save_dir)
	os.mkdir(save_dir +'/'+ states_dir)
	os.mkdir(save_dir +'/'+ posplots_dir)
	os.mkdir(save_dir +'/'+ gnf_data_dir)
	os.mkdir(save_dir +'/'+ gnf_data_dir +'/'+ cell_list_dir)

def make_new_details(L, rho, v0, delta, details_path):
	# Create a new TrajectoryDetails object
	td = TrajectoryDetails(details_path)
	td.L = L
	td.rho = rho
	td.v0 = v0
	td.delta = delta
	file = open(details_path, 'wb')
	pickle.dump(td, file)
	file.close()

def make_new_system(L, rho, v0, delta, hoomd_output_path):
	hoomd.context.initialize()

	# Define the initial conditions of the new system
	N = rho*L*L
	seed = random.randint(1,10e4)
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

	# Take one step forward in time to save an initial .gsd file
	nl_cell = nl.cell(r_buff=0)
	all = hoomd.group.all()
	hoomd.md.integrate.mode_standard(dt=1)
	hoomd.md.integrate.vicsek(all, nlist=nl_cell, r_cut=1, v0=v0, delta=delta, eta=0, seed=seed)
	d = hoomd.dump.gsd(hoomd_output_path, period=1, group=all, overwrite=True)
	hoomd.run(1)	

################################################################################################################

def get_trajectory_details(details_path):
	file = open(details_path, 'rb')
	trajectory_details = pickle.load(file)
	return trajectory_details

def update_trajectory_details(details_path, run_time):
	td = get_trajectory_details(details_path)
	td.time_integrated += run_time
	os.remove(details_path)
	file = open(details_path, 'wb')
	pickle.dump(td, file)
	file.close()
	print('\n\nSystem has now been evolved, %d time-steps' %(int(td.time_integrated)))

def equilibrate_system(v0, delta, run_time, hoomd_output_path, hoomd_output_backup_path, details_path):
	hoomd.context.initialize()

	# Make a backup copy of the simulation file in case of crash during integration
	shutil.copyfile(hoomd_output_path, hoomd_output_backup_path) 
	os.remove(hoomd_output_path)

	# Restore the previous state of the system
	hoomd.init.read_gsd(hoomd_output_backup_path, frame=-1)

	# Integrate the system forward in time
	seed = random.randint(1,10e4)
	nl_cell = nl.cell(r_buff=0)
	all = hoomd.group.all()
	hoomd.md.integrate.mode_standard(dt=1)
	hoomd.md.integrate.vicsek(all, nlist=nl_cell, r_cut=1, v0=v0, delta=delta, eta=0, seed=seed)
	d = hoomd.dump.gsd(hoomd_output_path, period=run_time-1, group=all, overwrite=True)
	hoomd.run(run_time)

	# Once the system has been integrated, update the integration time
	update_trajectory_details(details_path, run_time)

################################################################################################################

def collect_gnf_data(v0, delta, hoomd_output_path, gnf_output_path, data_collection_time):
	hoomd.context.initialize()

	# Restore the previous state of the system
	hoomd.init.read_gsd(hoomd_output_path, frame=-1)

	# Integrate the system forward in time
	seed = random.randint(1,10e4)
	nl_cell = nl.cell(r_buff=0)
	all = hoomd.group.all()
	hoomd.md.integrate.mode_standard(dt=1)
	hoomd.md.integrate.vicsek(all, nlist=nl_cell, r_cut=1, v0=v0, delta=delta, eta=0, seed=seed)
	d = hoomd.dump.gsd(gnf_output_path, period=1, group=all, overwrite=True)
	hoomd.run(data_collection_time)

def perform_gnf_calculations(L, rho, gnf_output_path, cell_sizes, cell_list_path):
	
	snapshots = gsd.hoomd.open(name=gnf_output_path, mode="rb")
	data_points = []

	for cell_size in cell_sizes:
		
		print("Working on cell-size: [%d,%d]\n" %(cell_size[0], cell_size[1]))
		cell_area = pow(L,2)/(cell_size[0]*cell_size[1])
		langle_n_rangle = rho*cell_area
		print("<n>=%f\n" %(langle_n_rangle))
		n_sq_sum = 0
		num_cells = 0

		num_time_steps = int(len(snapshots)/12)
		for t in range(0, num_time_steps):

			print("dt: ", t, "/", num_time_steps, end="\r")
			system = snapshots[t]
			cell_list = CellList(system=system, L=L, n_x=cell_size[0], n_y=cell_size[1])
			cell_list.make_cell_list()
			cell_list.populate_cell_list()

			for row in cell_list.cell_list:
				for cell in row:
					n_sq_sum+=pow(cell.num_particles, 2)	
					num_cells += 1
		
		langle_n_sq_rangle = n_sq_sum / num_cells
		std_dev = sqrt(langle_n_sq_rangle - pow(langle_n_rangle, 2))
		data_points.append([langle_n_rangle, std_dev, num_cells, str(cell_size[0])+'x'+str(cell_size[1])])

	return data_points

def save_gnf_values(gnf_values_path, data_points):
	csv_save_location = gnf_values_path
	with open(csv_save_location, mode='a+') as num_fluctuations_file:
		for i in range(0, len(data_points)):
			data_point = data_points[i]
			num_fluc_writer = csv.writer(num_fluctuations_file, delimiter=',')
			num_fluc_writer.writerow(data_point)

################################################################################################################

def retrieve_csv_data(gnf_values_path):
	x, y = [], []
	with open(gnf_values_path, 'r') as csvfile:
		data = csv.reader(csvfile, delimiter=',')
		for data_point in data:
			x.append(float(data_point[0]))
			y.append(float(data_point[1]))
	return x, y

def plot_gnf_values(gnf_values_path):
	from scipy import stats
	import numpy as np
	import math
	from math import log, pow
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	x_min, y_min, x_max, y_max = 1, 1, 1e4, 1e4

	x, y = retrieve_csv_data(gnf_values_path)
	for i in range(0, len(x)):
		ax.scatter(x[i], y[i], color='k', marker='o')

	ax.plot([x_min, x_max], [pow(x_min, 0.8), pow(x_max, 0.8)], ls='--', color='r')

	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(x_min, x_max)
	ax.set_xlabel('$\\langle{n}\\rangle$', fontsize=14)
	ax.xaxis.set_label_coords(1.05, -0.025)
	ax.set_ylabel('$\\Delta{n}$', fontsize=14, rotation=360)
	ax.yaxis.set_label_coords(-0.035, 0.97)
	plt.show()



################################################################################################################

def main():

	L = 256
	rho = 2
	v0 = 0.5
	delta = 0.25
	
	save_dir = 'reproduce_chate_ginelli_2008_gnf'
	states_dir = 'saved_states'
	posplots_dir = 'pos_plots'
	gnf_data_dir = 'gnf_data'
	cell_list_dir = 'cell_lists'

	hoomd_output_name = 'trajectory.gsd'
	hoomd_output_backup_name = 'trajectory_backup.gsd'
	details_name = 'details.obj'
	gnf_name = 'gnf.gsd'
	gnf_calc_name = 'gnf_calc.obj'
	gnf_csv_name = 'gnf_plot_values.csv'
	state_name = 'state_t'
	plot_name = 'pos_t'
	cell_list_name = 'CellSize_'
	

	#######################################################################################################

	hoomd_output_path = save_dir +'/'+ hoomd_output_name
	hoomd_output_backup_path = save_dir +'/'+ hoomd_output_backup_name

	details_path = save_dir +'/'+ details_name
	saved_states_path = save_dir +'/'+ states_dir +'/'+ state_name
	pos_plots_path = save_dir +'/'+ states_dir +'/'+ plot_name

	gnf_path = save_dir +'/'+ gnf_data_dir
	gnf_output_path = gnf_path +'/'+ gnf_name
	gnf_calc_path = gnf_path +'/'+ gnf_calc_name
	gnf_values_path = gnf_path +'/'+ gnf_csv_name
	
	cell_lists_path = save_dir +'/'+ gnf_data_dir +'/'+ cell_list_dir 
	cell_list_path = cell_lists_path +'/'+ cell_list_name

	#######################################################################################################

	#erase_old_data(save_dir, states_dir, posplots_dir, gnf_data_dir, cell_lists_path)
	#make_new_details(L, rho, v0, delta, details_path)
	#make_new_system(L, rho, v0, delta, hoomd_output_path)

	#######################################################################################################

	run_time = 2e4
	equilibrate_system(v0, delta, run_time, hoomd_output_path, hoomd_output_backup_path, details_path)

	#######################################################################################################

	data_collection_time = 2e3
	collect_gnf_data(v0, delta, hoomd_output_path, gnf_output_path, data_collection_time)

	#######################################################################################################

	#cell_sizes = [ [15,15] ]
	#data_points = perform_gnf_calculations(L, rho, gnf_output_path, cell_sizes, cell_list_path)
	#save_gnf_values(gnf_values_path, data_points)

	
	###################################################################################################
	
	plot_gnf_values(gnf_values_path)
	
	###################################################################################################
	

if __name__ == "__main__":
	system('clear')
	main()
