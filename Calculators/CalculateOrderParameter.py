
import math
from math import cos, sin, atan2, sqrt, pow, ceil
import matplotlib.pyplot as plt

class CalculateOrderParameter:

	def __init__(self):
		pass


	''' Method:     "get_order_parameter_time_evolution"
		input(s):    A list (run) of gsd snapshots 
		output(s):   A list (order_parameter_list) containing the instantaneous polar
				     order parameter of the system at each time step
		description: Given the entire trajectory of the system, this method retrieves
					 the orientation of each particle in the system at each time-step 
					 (dt), obtains the 2D-polar angle from the quaternion orientation, 
					 and then calculates the order parameter at the time-step (dt)
	'''
	def get_order_parameter_time_evolution(self, run):
		print('\n\nCalculating the time-evolution of the order parameter...')
		t = len(run)
		order_parameter_list = []
		for dt in range(0, t):
			orientations = run[dt].particles.orientation
			angles = self.convert_quaternions_to_angles(orientations)
			order_parameter = self.calculate_order_parameter(angles)
			order_parameter_list.append(order_parameter)
		return order_parameter_list

	def convert_quaternions_to_angles(self, quaternion_list):
		angles = []
		for quaternion in quaternion_list:
			angle = quaternion[0]
			angles.append(angle)
		return angles


	''' Method:     "calculate_order_parameter"
		input(s):    A list (angles) containing the 2D-polar angle [-pi,pi] of each
					 particle in the system at a given time-step
		output(s):   The global polar order parameter of the system at the time-step
		description: Given the 2D-polar angle of a particle 'i' in the system, a cartesian 
					 unit vector can be defined as
					 		vec{r}_i = [cos(theta_i)  sin(theta_i)]
					 the polar order of the system can then be found by calculating the 
					 ensemble average vector
					 		vec{R} = 1/N * sum_i vec{r}_i
					 and then taking its magnitude
					 		phi = sqrt[vec{R} * vec{R}]
	'''
	def calculate_order_parameter(self, angles):
		N = len(angles)
		x, y = 0, 0
		for i in range(0, N):
			x+=cos(angles[i])
			y+=sin(angles[i])
		order_parameter = sqrt(x*x + y*y) / N
		return order_parameter

	'''
		Method "calculate_global_avg_orientation":
		input(s):    A list (angles) containing the 2D-polar angle [-pi,pi] of each
					 particle in the system at a given time-step
		output(s):   The orientation of the system's velocity vector
		description: In a Vicsek model simulation, the system orientation unit vector is 
					 defined by the ensemble average
							hat{phi} = 1/N*sum_i [ x_i*hat{x} + y_i*hat{y} ] 
									 = x_phi*hat{x} + y_phi*hat{y}
					 This vector may be used to define new (orthogonal) basis vectors
	'''
	def calculate_global_avg_orientation(self, angles):
		N = len(angles)
		x, y = 0, 0
		for i in range(0, N):
			x+=cos(angles[i])
			y+=sin(angles[i])
		theta_global = atan2(y, x)
		return theta_global


	''' Method:     "calculate_order_parameter_steady_state"
		input(s):    (1) A list (order_parameter_list) of the system's global order parameter
					 at each time step in the simulation (2) the time-step at which to stop
					 averaging
		output(s):   The time-averaged global order parameter
		description: This method simply performs a time-averaging of the global order 
					 parameter phi(t) as
					 			<phi>_t = 1/t * sum_{dt} phi(dt)
	'''
	def calculate_time_averaged_order_parameter(self, order_parameter_list, t):
		sum = 0
		for dt in range(0, t):
			sum+=order_parameter_list[dt]
		time_averaged_order_parameter = sum/t
		return time_averaged_order_parameter


	''' Method:     "calculate_order_parameter_steady_state"
		input(s):    (1) A list (order_parameter_list)of the system's global order parameter
					 at each time step in the simulation (2) a scalar (epsilon) which defines
					 a range
					 		<phi>(infty)-epsilon < <phi>(t_ss) < <phi>(infty)+epsilon
					 in which <phi>(t_ss) is considered to be approximately equal to <phi>(infty)
		output(s):   (1) The time (t_ss) at which the instantaneous average polar order parameter
					 is approximately equal to <phi>(infty) (2) the value of <phi>(t_ss)
		description: In order to avoid errors from early fluctuations in <phi>(t), the minimum 
					 steady-state time to be 10% of the total run time. From this time-step, the 
					 instantaneous order parameter is compared to that at infinity; if it has not
					 converged, t_ss is incremented
	'''
	def calculate_order_parameter_steady_state(self, order_parameter_list, epsilon):
		# for large numbers of time-steps/particles this can be expensive: print status message
		print('\n\nCalculating the order parameter steady-state time...')
		taop_infinity = self.calculate_time_averaged_order_parameter(order_parameter_list, len(order_parameter_list))
		print("<phi>(infinity) = %2.8f" %(taop_infinity))
		taop_ss = -1
		# take the minimum steady-state time to be 10% of the total integration time
		t_ss = ceil(0.1*len(order_parameter_list))
		converged = False
		while(not converged):
			t_ss+=1
			taop_ss = self.calculate_time_averaged_order_parameter(order_parameter_list, t_ss)
			converged = (taop_infinity-epsilon)<taop_ss and taop_ss<(taop_infinity+epsilon)
		print("<phi>(t_ss) = %2.8f" %(taop_ss))
		print("t_ss = %d/%d\n" %(t_ss, len(order_parameter_list)))
		return t_ss, taop_ss

	def plot_order_parameter_time_evolution(self, order_parameter_list):
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		ax.plot(order_parameter_list, color='k', ls='--')
		ax.set_xlim(0.00, len(order_parameter_list))
		ax.set_ylim(0.00, 1.10)
		ax.set_xlabel('$t$')
		ax.set_ylabel('$\\phi$')
		plt.show()
