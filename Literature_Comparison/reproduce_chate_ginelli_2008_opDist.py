from RunHoomd import RunHoomd
from CalculateOrderParameter import CalculateOrderParameter
from CalculateBinderCumulant import CalculateBinderCumulant
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

def save_order_parameter_list(simulation, order_parameter_list_path):
	# get gsd snapshots
	simulation.collect_data()
	data = simulation.retrieve_data()

	# calculate the order parameter at each timestep, and the time-averaged order parameter
	calc_op = CalculateOrderParameter()
	order_parameter_list = calc_op.get_order_parameter_time_evolution(data)	

	data_points = []
	for dt in range(0, len(order_parameter_list)):
		data_points.append([dt, order_parameter_list[dt]])

	# save the calculated taop value for the given noise
	with open(order_parameter_list_path, 'w+', newline="") as file:
		writer = csv.writer(file, delimiter=',')
		writer.writerows(data_points)

	# to save hard-drive memory, delete the .gsd simulation file once all calculations are complete
	if path.exists(simulation.filename):
		os.remove(simulation.filename)

def save_order_parameter_distribution(order_parameter_list_path, order_parameter_distribution_path):
	#
	time_steps, order_parameter_list = retrieve_csv_data(order_parameter_list_path)
	calc_bc = CalculateBinderCumulant()
	calc_bc.data = order_parameter_list
	distribution, bin_edges = calc_bc.make_distribution(100)

	data_points = []
	for phi in range(0, len(bin_edges)-1):
		data_points.append([bin_edges[phi], distribution[phi]])

	# save the calculated taop value for the given noise
	with open(order_parameter_distribution_path, 'w+', newline="") as file:
		writer = csv.writer(file, delimiter=',')
		writer.writerows(data_points)

def plot_comparison(lit_dir_loc, vn_order_parameter_distribution_path, sn_order_parameter_distribution_path):
	# set up the plot variables
	import matplotlib
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	# define the locations of the files whose data is to be plotted
	literature_location_1 = lit_dir_loc + 'collective_motion_without_cohesion/Lit_Scalar_Dist.csv'
	literature_location_2 = lit_dir_loc + 'collective_motion_without_cohesion/Lit_Vectorial_Dist.csv'
	experiment_location_1 = sn_order_parameter_distribution_path
	experiment_location_2 = vn_order_parameter_distribution_path
	
	# retrieve the data to be plotted
	lit_scalar_phis, lit_scalar_probs = retrieve_csv_data(literature_location_1)
	lit_vectorial_phis, lit_vectorial_probs = retrieve_csv_data(literature_location_2)
	exp_scalar_phis, exp_scalar_probs = retrieve_csv_data(experiment_location_1)
	exp_vectorial_phis, exp_vectorial_probs = retrieve_csv_data(experiment_location_2)
	
	for i in range(0, len(exp_scalar_probs)):
		exp_scalar_probs[i] = exp_scalar_probs[i]*100
	for i in range(0, len(exp_vectorial_probs)):
		exp_vectorial_probs[i] = exp_vectorial_probs[i]*100

	# add the data to be plotted
	ax.scatter(lit_scalar_phis, lit_scalar_probs, color='r', marker='x', label='Literature: Scalar')
	ax.scatter(lit_vectorial_phis, lit_vectorial_probs, color='r', marker='s', label='Literature: Vectorial')
	ax.plot(exp_scalar_phis, exp_scalar_probs, color='k', ls='solid', label='TwoStepVicsek: Scalar')
	ax.plot(exp_vectorial_phis, exp_vectorial_probs, color='k', ls='--', label='TwoStepVicsek: Vectorial')
	
	# configure and plot
	ax.set_xlim(0.00, 0.90)
	ax.set_xlabel('$\\phi$', fontsize=18)
	ax.set_ylabel('$P(\\phi)$', rotation=0, fontsize=20)
	#ax.legend()
	plt.show()


def main():

	#################################################################################################################

	'''
	Note: In Chate, Ginelli 2008, run time was 3e5. 
	'''
	L = 64
	rho = 2
	v0 = 0.5
	run_time = int(1e5)
	eta_t = 0.61468 #0.4146 is complete ordered
	delta_t = 0.459

	#################################################################################################################

	primary_dir = 'reproduce_chate_ginelli_2008_opDist'
	sn_dir = primary_dir + '/scalar_noise'
	vn_dir = primary_dir + '/vectorial_noise'
	opList_filename = '/order_parameter_list.csv'
	opDist_filename = '/order_parameter_distribution.csv'
	lit_dir_loc = 'Literature/'
	gsd_output_path = primary_dir + '/trajectory.gsd'

	#################################################################################################################
	
	#erase_old_data(primary_dir, sn_dir, vn_dir)
	
	#################################################################################################################
	
	simulation_vectorial = RunHoomd(L=L, rho=rho, eta=eta_t, v0=v0, num_time_steps=run_time, filename=gsd_output_path)
	order_parameter_list_path = vn_dir + opList_filename
	vn_order_parameter_distribution_path = vn_dir + opDist_filename
	save_order_parameter_list(simulation_vectorial, order_parameter_list_path)
	save_order_parameter_distribution(order_parameter_list_path, vn_order_parameter_distribution_path)
	
	#################################################################################################################
	
	simulation_scalar = RunHoomd(L=L, rho=rho, delta=delta_t, v0=v0, num_time_steps=run_time, filename=gsd_output_path)
	order_parameter_list_path = sn_dir + opList_filename
	sn_order_parameter_distribution_path = sn_dir + opDist_filename
	save_order_parameter_list(simulation_scalar, order_parameter_list_path)
	save_order_parameter_distribution(order_parameter_list_path, sn_order_parameter_distribution_path)
	
	#################################################################################################################

	plot_comparison(lit_dir_loc, vn_order_parameter_distribution_path, sn_order_parameter_distribution_path)

	#################################################################################################################

if __name__ == "__main__":
	system('clear')
	main()

