from RunHoomd import RunHoomd
from CalculateOrderParameter import CalculateOrderParameter
import matplotlib.pyplot as plt
import os
from os import system, path
import shutil
import csv


def erase_old_data(save_directory, sn_dir, vn_dir):
	if path.exists(save_directory):
		shutil.rmtree(save_directory)
	os.mkdir(save_directory)
	os.mkdir(sn_dir)
	os.mkdir(vn_dir)

def retrieve_csv_data(open_location):
	data0 = []
	data1 = []
	with open(open_location, 'r') as csvfile:
		data = csv.reader(csvfile, delimiter=',')
		for data_point in data:
			data0.append(float(data_point[0]))
			data1.append(float(data_point[1]))
	return data0, data1

def plot_comparison(lit_dir_loc, sn_dir, vn_dir, taop_filename):
	# set up the plot variables
	import matplotlib
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	# define the locations of the files whose data is to be plotted
	literature_location_1 = lit_dir_loc + 'onset_of_collective_and_cohesive_motion/Lit_Scalar_TAOP.csv'
	literature_location_2 = lit_dir_loc + 'onset_of_collective_and_cohesive_motion/Lit_Vectorial_TAOP.csv'
	experiment_location_1 = sn_dir + taop_filename
	experiment_location_2 = vn_dir + taop_filename
	
	# retrieve the data to be plotted
	lit_scalar_noises, lit_scalar_taops = retrieve_csv_data(literature_location_1)
	lit_vectorial_noises, lit_vectorial_taops = retrieve_csv_data(literature_location_2)
	exp_scalar_noises, exp_scalar_taops = retrieve_csv_data(experiment_location_1)
	exp_vectorial_noises, exp_vectorial_taops = retrieve_csv_data(experiment_location_2)
	
	# add the data to be plotted
	ax.scatter(lit_scalar_noises, lit_scalar_taops, color='r', marker='x', label='Literature: Scalar')
	ax.scatter(lit_vectorial_noises, lit_vectorial_taops, color='r', marker='s', label='Literature: Vectorial')
	ax.plot(exp_scalar_noises, exp_scalar_taops, color='k', ls='solid', label='TwoStepVicsek: Scalar')
	ax.plot(exp_vectorial_noises, exp_vectorial_taops, color='k', ls='--', label='TwoStepVicsek: Vectorial')
	
	# configure and plot
	ax.set_xlim(0.20, 0.80)
	ax.set_ylim(0.00, 1.00)
	ax.set_xlabel('$\\delta$', fontsize=24)
	ax.set_ylabel('$\\langle{\\phi}\\rangle_{t}$', rotation=0, fontsize=20)
	ax.xaxis.set_label_coords(0.90, -0.05)
	ax.yaxis.set_label_coords(-0.1, 0.5)
	#ax.legend(prop={'size': 14})
	plt.show()


####################################################################################################################
	

def main():

	# System parameters used by Ginelli and Chate
	L = 32
	rho = 2
	v0 = 0.5

	#################################################################################################################

	# setting this variable to true will erase all previous data
	new_experiment = False

	# noise intensities can be added to this list, and corresponding taops will be appending to the saved data set
	noises =[]
	vectorial = True
	
	# how many time steps to integrate for each noise intensity
	run_time = int(5e4)

	# each time new experiment is run, these directories will be created 
	primary_dir = 'reproduce_gregoire_chate_2004_taop'
	sn_dir = primary_dir + '/scalar_noise'
	vn_dir = primary_dir + '/vectorial_noise'

	# define the save file names
	taop_filename = '/polar_order_parameter.csv'
	gsd_output_path = primary_dir + '/trajectory.gsd'

	# this variable holds the path to the directory 'Literature'
	lit_dir_loc = 'Literature/'

	if new_experiment:
		erase_old_data(primary_dir, sn_dir, vn_dir)

	#################################################################################################################
	
	for noise in noises:

		# run commands are dependent on type of noise being used
		simulation = RunHoomd(L=L, rho=rho, v0=v0, num_time_steps=run_time, filename=gsd_output_path)
		if vectorial:
			simulation.eta = noise
			save_loc = vn_dir 
		else:
			simulation.delta = noise
			save_loc = sn_dir 

		# get gsd snapshots
		simulation.collect_data()
		data = simulation.retrieve_data()

		# calculate the order parameter at each timestep, and the time-averaged order parameter
		calc_op = CalculateOrderParameter()
		order_parameter_list = calc_op.get_order_parameter_time_evolution(data)	
		taop = calc_op.calculate_time_averaged_order_parameter(order_parameter_list, run_time)

		# save the calculated taop value for the given noise
		with open(save_loc + taop_filename, 'a+', newline="") as file:
			writer = csv.writer(file, delimiter=',')
			writer.writerows([[noise, taop]])

		# for neatness purposes, delete the .gsd simulation file once all calculations are complete
		if path.exists(simulation.filename):
			os.remove(simulation.filename)
	
	#################################################################################################################

	plot_comparison(lit_dir_loc, sn_dir, vn_dir, taop_filename)

	#################################################################################################################

if __name__ == "__main__":
	system('clear')
	main()

