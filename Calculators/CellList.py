import math
from math import cos, sin, atan2, pi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.lines import Line2D
from matplotlib.pyplot import arrow

'''
	Single Cell Data Structure:
		Method "__init__": 
			A sub-volume cell of the simulation box is defined by,
			1) a vertex defining one corner of the cell, (x_origin, y_origin)
			2) the upper-right corner of the cell defined by the horizontal distance x_max
			3) the lower-left corner of the cell defined by the vertical distance y_max
		Method "is_in_cell": 
			This method takes in the (x,y) coordinate of a single particle and returns 
			True (False) if the particle is (not) in the cell. 
		Method "is_in_row":
			This method is used to increased the efficiency of sorting through a cell-list. 
			See the CellList class to see functionality  
'''
class Cell:
	
	def __init__(self, x_origin=0, y_origin=0, x_max=0, y_min=0):
		self.x_origin = x_origin
		self.y_origin = y_origin
		self.x_max = x_max
		self.y_min = y_min
		self.particles = []
		self.num_particles = 0

	def is_in_cell(self, x, y):
		in_x = self.x_origin<x and x<self.x_max
		in_y = self.y_min<y and y<self.y_origin
		if(in_x and in_y):
			return True
		else:
			return False

	def is_in_row(self, y):
		in_row = self.y_min<y and y<self.y_origin
		if in_row:
			return True
		else:
			return False



class CellList:

	'''	Method 		 "__init__":
		Description: A simulation box in Hoomd is defined by a linear length L with Cartesian coordinates 
					 [-L/2, L/2]. The total LxL box can be divided into n_x*n_y cells whose linear 
					 dimensions are l_i = L/n_i. The CellList class uses (l_x,l_y) to map the simulation
					 box into a 2D array of cells with rows spaced by l_y and columns spaced by l_x. 
	'''
	def __init__(self, L=1, n_x=1, n_y=1, system=[], theta=0, plot=False):
		# essential system variables
		self.L = L
		self.l_x = L/n_x
		self.l_y = L/n_y
		self.system = system
		self.cell_list = []

		self.theta = theta
		# rotated particle positions are stored as ordered pairs [x,y] in this list
		self.rotated_particles = []
		if plot:
			fig = plt.figure()
			self.ax = fig.add_subplot(1, 1, 1)

	'''	Method 		 "make_cell_list":
		Description: The 2D array of cells is created by defined by a cell coordinate (y_row, x_column), 
					 beginning in the upper-left corner of the simulation box with coordinate (-L/2,L/2). 
					 The n_x columns of width l_x=L/n_x are then looped through until x_column=L/2; the
					 n_y rows of height l_y=L/n_y are then looped through until y_row=-L/2. 
	'''
	def make_cell_list(self):
		y_row = self.L/2
		while not (y_row <= -self.L/2):
			row_cell_list = []
			x_column = -self.L/2
			while not (x_column >= self.L/2):
				x_max, y_min = (x_column+self.l_x), (y_row-self.l_y)
				cell = Cell(x_origin=x_column, y_origin=y_row, x_max=x_max, y_min=y_min)
				row_cell_list.append(cell)
				x_column += self.l_x
			self.cell_list.append(row_cell_list)
			y_row -= self.l_y

	'''	Method: 	 "populate_cell_list"
		Description: This method uses the positions of the particles of the system, at some time t>t_ss, in order to 
					 count the number of particles within each cell in cell_list. This process is made more efficient 
					 through two conditions 
					 (1) each particle can only be located in one cell, so once a particle is found to belong in a cell 
					     stop searching the cell list and move onto the next particle 
					 (2) the algorithm loops top-down through rows of cells, comparing the particle with each cell in the 
						 row. If a quick check of the particle's y-coordinate shows that it isn't in the current row, then 
						 there is no need to search the row's cells, and we may move onto the next row 
	'''
	def populate_cell_list(self):
		# add the xy-position of each particle in the simulation box
		particles = self.system.particles.position
		for	i in range(0, len(particles)):
			particle = particles[i]
			r = [particle[0], particle[1]]
			if (self.theta != 0):
				r_prime = self.rotation_matrix(r)
				self.rotated_particles.append(r_prime)
				r = r_prime
			# and determine which particles are in which cells
			cell_found = False
			for row in self.cell_list:
				for cell in row:
					if not cell.is_in_row(r[1]):
						break
					if cell.is_in_cell(r[0],r[1]):
						cell.num_particles += 1
						cell.particles.append(i)
						cell_found = True
				if cell_found:
					break


	'''	Method: 	 "rotation_matrix"
		Description: If the cells are to be oriented at some angle (theta) counter-clockwise to the xy-plane, 
					 then "theta" can be initialized. But, rather than rotating the cell list by theta, it is 
					 simpler to rotate the position vector of each particle according to
								[x']   [cos(theta)  -sin(theta)] [x]
								[y'] = [sin(theta)   cos(theta)] [y]
					 Period boundary conditions are then applied to the primed vector to keep the particle within the 
					 simulation box. 
	'''
	def rotation_matrix(self, r):
		# performing a counter-clockwise rotation about the z-axis
		x, y = r[0], r[1]
		x_prime = x*cos(self.theta)-y*sin(self.theta)
		y_prime = x*sin(self.theta)+y*cos(self.theta) 
		x_prime, y_prime = self.apply_boundary_conditions(x_prime), self.apply_boundary_conditions(y_prime)
		return [x_prime, y_prime]

	def apply_boundary_conditions(self, q):
		# assuming that the simulation box is symmetric, q=x,y should have same b.c.
		if q>self.L/2:
			q -=self.L
		if q<-self.L/2:
			q +=self.L
		return q

	'''	Method: 	 "plot_cell_list"
		Description: 
	'''
	def plot_cell_list(self):
		# the plot window is set to the simulation box size
		self.ax.set_xlim(-self.L/2, self.L/2)
		self.ax.set_ylim(-self.L/2, self.L/2)
		self.ax.legend()
		plt.show()

	def add_system_boundaries(self):
		# plot an outline boundary of the simulation box
		boundary_lines = []
		left = [[-self.L/2,-self.L/2], [self.L/2,-self.L/2]]
		right = [[self.L/2,self.L/2], [self.L/2,-self.L/2]]
		top = [[-self.L/2,self.L/2], [self.L/2,self.L/2]]
		bottom = [[-self.L/2,self.L/2], [-self.L/2,-self.L/2]]
		boundary_lines.append(left)
		boundary_lines.append(right)
		boundary_lines.append(bottom)
		boundary_lines.append(top)
		for boundary in boundary_lines:
			self.ax.add_line(Line2D(boundary[0], boundary[1], linewidth=3, color='black'))

	def add_cell_vertices(self):
		for row in self.cell_list:
			for cell in row:
				self.ax.scatter(cell.x_origin, cell.y_origin, color='black', marker='s')

	def add_cell_grid(self):
		for row in self.cell_list:
			for cell in row:
				# in order to plot with Line2D, an array of x's and y's must be made separately
				hor_line_xs = [cell.x_origin, cell.x_max]
				hor_line_ys = [cell.y_origin, cell.y_origin]
				self.ax.add_line(Line2D(hor_line_xs, hor_line_ys, linewidth=1, linestyle='--', color='red'))
				ver_line_xs = [cell.x_origin, cell.x_origin]
				ver_line_ys = [cell.y_origin, cell.y_min]
				self.ax.add_line(Line2D(ver_line_xs, ver_line_ys, linewidth=1, linestyle='--', color='red'))

	def add_axes(self):
		arr_size = 0.2
		l = (self.L/2)*0.90
		self.ax.arrow(0, 0, 0, l, fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('+y', xy=(arr_size,l))
		self.ax.arrow(0, 0, 0, -l, fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('-y', xy=(arr_size,-l))
		self.ax.arrow(0, 0, l, 0, fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('+x', xy=(l,arr_size))
		self.ax.arrow(0, 0, -l, 0, fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('-x', xy=(-l,arr_size))

	def add_particle_positions(self):
		# plot the original position of each particle
		for row in self.cell_list:
			for cell in row:
				i = 0
				for particle in cell.particles:
					print("plotting particle: ", i, "/", len(cell.particles), end="\r")
					position = self.system.particles.position[particle]
					#self.ax.scatter(position[0], position[1], color='k', marker='o')  # for small number of particles
					self.ax.scatter(position[0], position[1], color='k', marker='.', s=0.1)  # for large number of particles
					i += 1

	def add_particle_velocities(self):
		for i in range(0, len(self.system.particles.position)):
			position = self.system.particles.position[i]
			x0, y0 = position[0], position[1]
			#self.ax.scatter(x0, y0, color='k', marker='o')  # for small number of particles
			self.ax.scatter(x0, y0, color='k', marker='.', s=0.1)  # for large number of particles
			orientation = self.system.particles.orientation[i]
			theta = orientation[0]
			self.ax.arrow(x0,y0, cos(theta),sin(theta), fc="k", ec="k", head_width=0.1, head_length=0.1)
		
	def add_rotated_axes(self):
		arr_size = 0.2
		x_hat, y_hat = [1, 0], [0, 1]
		phi_perp = self.rotation_matrix(x_hat)
		phi_par = self.rotation_matrix(y_hat)
		s = '$\\phi_{||}$'
		l = (self.L/2)*0.90
		self.ax.arrow(0, 0, l*phi_perp[0], l*phi_perp[1], fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('+$\\phi_{||}$', xy=(l*phi_perp[0]+arr_size,l*phi_perp[1]))
		self.ax.arrow(0, 0, -l*phi_perp[0], -l*phi_perp[1], fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('-$\\phi_{||}$', xy=(-l*phi_perp[0]+arr_size,-l*phi_perp[1]))
		self.ax.arrow(0, 0, l*phi_par[0], l*phi_par[1], fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('-$\\phi_{\\perp}$', xy=(l*phi_par[0],l*phi_par[1]-arr_size))
		self.ax.arrow(0, 0, -l*phi_par[0], -l*phi_par[1], fc="k", ec="k", head_width=arr_size, head_length=arr_size)
		self.ax.annotate('+$\\phi_{\\perp}$', xy=(-l*phi_par[0],-l*phi_par[1]+arr_size))

	def add_rotated_particle_positions(self):
		# plot the rotated position of each particle
		for particle in self.rotated_particles:
			self.ax.scatter(particle[0], particle[1], color='k', marker='s')

	def add_rotated_velocities(self):
		for i in range(0, len(self.rotated_particles)):
			position = self.rotated_particles[i]
			x0, y0 = position[0], position[1]
			orientation = self.system.particles.orientation[i]
			theta = orientation[0]
			v = self.rotation_matrix([cos(theta),sin(theta)])
			self.ax.arrow(x0,y0, v[0],v[1], fc="k", ec="k", head_width=0.1, head_length=0.1)

	def add_com(self):
		R_x, R_y = 0, 0
		N = len(self.system.particles.position)
		for position in self.system.particles.position:
			R_x += position[0]/N
			R_y += position[1]/N
		self.ax.scatter(R_x, R_y, color='m', marker='s')
		self.ax.arrow(R_x, R_y, cos(self.theta), sin(self.theta), fc="m", ec="m", head_width=0.1, head_length=0.1)		


# test method to ensure that the classes and methods work properly
def main():
	import hoomd, hoomd.md
	from InitializeUnitCell import InitializeUnitCell
	from CalculateOrderParameter import CalculateOrderParameter

	# Choose the number of particles and linear size of the simulation box
	N, L = 10, 10
	initialize = InitializeUnitCell(N, L)

	# Define the number of cells to divide the system into, and calculate their linear dimensions 
	n_x, n_y = 2, 200
	l_x, l_y = L/n_x, L/n_y
	non_random = True

	# Choose positions and orientations for each of the particles
	if non_random:
		pos_1 = [l_x*1, l_y*(1+1/2), 0]
		pos_2 = [l_x*1, -l_y*(1+1/2), 0]
		positions = [pos_1, pos_2]
		ori_1 = [pi/2, 0, 0, 1]
		ori_2 = [0, 0, 0, 1] 
		orientations = [ori_1, ori_2]
		N = len(positions)
	else:
		positions = initialize.randomized_positions()
		orientations = initialize.randomized_orientations()

	# Create the Hoomd simulation box
	hoomd.context.initialize()
	unit_cell = hoomd.lattice.unitcell(
	    N=N, 
	    a1=[L, 0, 0],
	    a2=[0, L, 0],
	    a3=[0, 0, 1],  
	    dimensions=2,
	    position=positions,
	    orientation=orientations)
	system = hoomd.init.create_lattice(unit_cell, n=1)
	snap = system.take_snapshot()

	'''
		For the defined system, take the quaternion orientations of each particle in the
		snapshot and calculate the global average orientation (theta)
	'''
	op_calc = CalculateOrderParameter()
	orientation_angles = op_calc.convert_quaternions_to_angles(snap.particles.orientation)
	theta = op_calc.calculate_global_avg_orientation(orientation_angles)

	# Define a cell list for the simulation box and populate the cells
	cell_list = CellList(L=L, n_x=n_x, n_y=n_y, system=snap, theta=theta, plot=True)
	cell_list.make_cell_list()
	cell_list.populate_cell_list()

	# Plot the cell list
	cell_list.add_system_boundaries()
	cell_list.add_axes()
	#cell_list.add_cell_vertices()
	cell_list.add_cell_grid()
	cell_list.add_particle_positions()
	cell_list.add_particle_velocities()

	cell_list.add_rotated_axes()
	cell_list.add_rotated_particle_positions()
	cell_list.add_com()
	cell_list.plot_cell_list()


if __name__ == "__main__":
	main()





