import math
import numpy as np
from Config import Config


class MeanFlow:
	
	def __init__(self):

		self.Configure()

		self.ComputeMeanFlowDensity()
		self.ComputeMeanFlowAxialVelocity()
		self.ComputeMeanFlowSoundSpeed()
		self.ComputeMeanFlowMachNumber()
		self.ComputeMeanFlowPressure()

	def Configure(self):
		self.config = Config()

		self.N 					= self.config.NNxx
		self.F 					= self.config.F
		self.E 					= self.config.E		
		self.R1					= self.config.R1
		self.R2					= self.config.R2
		self.maxit 				= self.config.maxit
		self.tol 				= self.config.tol
		self.dens0 				= self.config.dens0
		self.gamma 				= self.config.gamma  
	
		# This doesn't do anything as far as the MF is concerned! To modify.		            
		self.dir_coeff 			= self.config.dir_coeff 						

	def ComputeMeanFlowDensity(self):
		N = self.N
		F = self.F
		E = self.E
		R1 = self.R1
		R2 = self.R2
		maxit = self.maxit
		tol = self.tol
		dens0 = self.dens0
		gamma = self.gamma	
		D0 = np.zeros(shape=(N,),dtype=float)
		d0 = np.zeros(shape=(maxit,),dtype=float)
		d0[0] = dens0
		
		k = 0
		while (k < N):
			j = 0
			beta = 0.5*math.pow(F,2)/math.pow((math.pow(R2[k],2)-math.pow(R1[k],2)),2);
			while( j < maxit - 1 ):
			
				f = beta*math.pow(d0[j],-2) + (1/(gamma-1))*math.pow(d0[j],gamma-1) - E
				df = -2*beta*math.pow(d0[j],-3)	+ math.pow(d0[j],gamma-2)

				# Newton's Method	
				d0[j+1] = d0[j] - f/df
					
				a = abs(d0[j+1]-d0[j])

				if(a < tol):
					#print("Newton Algorithm Converged")
					D0[k] = d0[j+1]
					d0 = np.zeros(shape=(maxit,),dtype=float)
					d0[0] = D0[k]
					break
				if(j==maxit-2):
					print("Newton Algorithm Failed to Converge")	
				j = j + 1	
			
			k=k+1
		
		self.D0 = D0
		
		

	def ComputeMeanFlowAxialVelocity(self):
		N = self.N
		R1 = self.R1
		R2 = self.R2
		F = self.F
		D0 = self.D0

		U0 = np.zeros(shape=(N,),dtype=float)
		U0denominator = np.power(R2,2) - np.power(R1,2)
		U0denominator = U0denominator*np.multiply(D0,U0denominator)
		U0 = np.divide(F,U0denominator)	
		self.U0 = U0	


	def ComputeMeanFlowSoundSpeed(self):
		N = self.N
		D0 = self.D0
		gamma = self.gamma
		C0 = np.zeros(shape=(N,),dtype=float)
		C0 = np.power(D0,(gamma-1)/2)
		self.C0 = C0
	

	def ComputeMeanFlowMachNumber(self):
		U0 = self.U0
		C0 = self.C0
		Mach0 = np.divide(U0,C0)
		self.Mach0 = Mach0

	
	def ComputeMeanFlowPressure(self):
		gamma = self.gamma
		D0 = self.D0
		P0 = np.multiply(1/gamma,np.power(D0,gamma))
		self.P0 = P0			