import math   
import numpy as np

class Config:

	def __init__(self):
 		
 		self.MeanFlow()
 		self.AcousticField()
 		self.Computations()
 		self.Grid()	
 		self.Plots()
 		self.ComplexPlane()
 		self.Preliminaries()
 		self.R1()
 		self.R2()

	def R1(self):
		x = self.x
		R1 = 0.689 - np.power((0.055+1.131*np.power((1-x/2),2)),0.5)
		# # Set negative entries equal to zero.
		R1[ R1 < 0 ] = 0 
		#R1 [ R1 > 0] = 0 

		self.R1 = R1

	def R2(self):
		x = self.x		
		self.R2 = 1.073-0.198*np.power((1-x/2),2) + 0.109*np.exp(-11*x/2)	
 		

 	def MeanFlow(self):
		
		self.dens0 = 0.9		
		
		# Mass Flow Rate
		self.F = 0.559
		#self.F = 0
		
		# Needs changing depending on the problem geometry
		self.E = 2.514
		
		# Ratio of specific heats for air
		self.gamma = 1.4202              
		self.dir_coeff = 1

	def AcousticField(self):	
	
		self.invZ1 = 0
		self.invZ2 = 1/(2-1j)

		self.m = 26
		self.omega = 25

	def Plots(self):

		self.PlotTracedModes = 1
		self.PlotDuct = 0
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
		self.tol = math.pow(10,-12)

	def ComplexPlane(self):


		# Defines the max/min of the complex plane
		self.mod_z_max = 10

		# Used for the number of points in the complex plane
		self.n = 20
		self.xshift = 35
		self.yshift = 0

	def Preliminaries(self):
		xmin = self.xmin
		xmax = self.xmax
		NNxx = self.NNxx

		# The increment/step for the initial modetracer (i.e., starting from F = 0)	
		self.dF = np.sign(self.F)*(0.01)

		# Set up the x stations
		self.x = np.linspace(xmin,xmax,NNxx)
		self.DeltaX = (xmax-xmin)/NNxx		



   	