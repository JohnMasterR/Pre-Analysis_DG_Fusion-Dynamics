from libc.math cimport tan, pi

cpdef double wrapAngle(double ang):
	x = ang % (2*pi)
	return x

# Particle class
cdef class particle:
	''' class particle: holds the position and charge of a particle

			Attributes
			==========
			charge : float
				Charge of the particle
			position : float
				Angular position of the particle
			step : float
				Step to increment the particle position. Computed externally from the dynamics
			
			Methods
			=======
			update_position()
				updates the position by incrementing by step. Resets to 0 the step. '''
	
	cdef public double charge
	cdef public double position
	cdef public double step	
	
	def __init__(self, charge, position):
		self.charge = charge
		self.position = position
		self.step = 0.0
	def __repr__(self):
		return repr((self.charge, self.position))
	def update_position(self):
		self.position = wrapAngle(self.position+self.step)
		self.step = 0.0

cpdef double force(int i, list particles):
	''' Computes the force on the i-th particle
			
			Parameters
			==========
			i : int
				Particle number
			particles : list
				List of particle class
			Returns : float
				Force on the particle 1 '''
	
	cdef double sum_i
	cdef int j
	
	sum_i = 0
	for j in range(len(particles)):
		if i != j:
			sum_i = sum_i + particles[i].charge*particles[j].charge/tan(0.5*(particles[i].position - particles[j].position))
	return 0.5*sum_i

cpdef void stepup(list particles, normalnoise, freediffusion, double dt):
	cdef int i
	cdef double mu_i
	cdef double force_i
	
	for i in range(len(particles)):
		if freediffusion:
			mu_i = 0
		else:
			# Sum over all charges to get electric force (E(theta))
			force_i = force(i, particles)
			mu_i = force_i*dt
		
		# Generate the random step base on the latter values of the positions
		particles[i].step = mu_i + normalnoise[i]
	
	# End (compute steps)



















