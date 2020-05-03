import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from Config import Config

class Eigensolver:

	def __init__(self):
		
		self.Configure()
		self.Preliminaries()
		self.R1()
		self.R2()
		self.Plots()

	def Configure(self):
		self.config 			= Config()
		self.dens0 				= self.config.dens0
		self.F 					= self.config.F
		self.E 					= self.config.E
		self.gamma 				= self.config.gamma              
		self.dir_coeff 			= self.config.dir_coeff
		self.Z1 				= self.config.Z1
		self.Z2 				= self.config.Z2
		self.m 					= self.config.m
		self.omega 				= self.config.omega
		self.PlotTracedModes 	= self.config.PlotTracedModes
		self.PlotDuct 			= self.config.PlotDuct
		self.PlotMeanFlow 		= self.config.PlotMeanFlow
		self.PlotAlpha 			= self.config.PlotAlpha
		self.PlotMu 			= self.config.PlotMu
		self.PlotSigma 			= self.config.PlotSigma
		self.PlotAmplitude 		= self.config.PlotAmplitude
		self.PlotContours 		= self.config.PlotContours
		self.NNxx 				= self.config.NNxx
		self.xmin 				= self.config.xmin
		self.xmax 				= self.config.xmax		
		self.maxit 				= self.config.maxit
		self.tol 				= self.config.tol
		self.mod_z_max 			= self.config.mod_z_max
		self.n 					= self.config.n
		self.xmin 				= self.config.xmin
		self.xmax 				= self.config.xmax	
	
	def Preliminaries(self):
		Z1 = self.Z1
		Z2 = self.Z2
		F  = self.F
		xmin = self.xmin
		xmax = self.xmax
		NNxx = self.NNxx

		self.invZ1 = 1/(Z1)
		self.invZ2 = 1/(Z2)

		# The increment/step for the initial modetracer (i.e., starting from F = 0)	
		self.dF = np.sign(F)*(0.01)

		# Set up the x stations
		
		self.x = np.linspace(xmin,xmax,NNxx)
		self.DeltaX = (xmax-xmin)/NNxx

	def R1(self):
		x = self.x
		R1 = 0.689 - np.power((0.055+1.131*np.power((1-x/2),2)),0.5)
		# Set negative entries equal to zero.
		R1[ R1 < 0 ] = 0 
		self.R1 = R1

	def R2(self):
		x = self.x		
		self.R2 = 1.073-0.198*np.power((1-x/2),2) + 0.109*np.exp(-11*x/2)	
	
	def Plots(self):
		x = self.x
		R1 = self.R1
		R2 = self.R2

		#print(XX.shape)

		fig, ax = plt.subplots()
		ax.plot(x, R1,linewidth=3)
		ax.plot(x,R2,linewidth=3)

		ax.set(xlabel='$x$', ylabel='$r$',title='Plot of the Duct')
		ax.xaxis.label.set_size(25)
		ax.yaxis.label.set_size(25)
		ax.grid()

		fig.savefig("duct.png")
		plt.show()	

		