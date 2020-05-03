class Config:

	def __init__(self):
 		
 		self.MeanFlow()
 		self.AcousticField()
 		self.Computations()
 		self.Grid()	
 		self.Plots()
 		self.ComplexPlane()

 	def MeanFlow(self):
		
		self.dens0 = 0.9		
		# Mass Flow Rate
		self.F = 0.559
		# Needs changing depending on the problem geometry
		self.E = 2.514
		# Ratio of specific heats for air
		self.gamma = 1.4202              
		self.dir_coeff = 1

	def AcousticField(self):	
	
		self.Z1 = 10^6
		self.Z2 = 2-1j
		self.m = 26
		self.omega = 25

	def Plots(self):

		self.PlotTracedModes = 1
		self.PlotDuct = 1
		self.PlotMeanFlow = 0
		self.PlotAlpha = 0
		self.PlotMu = 0
		self.PlotSigma = 0
		self.PlotAmplitude = 0
		self.PlotContours = 0

	def Grid(self):

		# The number of x stations	
		self.NNxx = 100
		self.xmin = 0
		self.xmax = 2	
	

	def Computations(self):	
	
		self.maxit = 100
		self.tol = 10^(-6)

	def ComplexPlane(self):


		# Defines the max/min of the complex plane
		self.mod_z_max = 3

		# Used for the number of points in the complex plane
		self.n = 40


   	