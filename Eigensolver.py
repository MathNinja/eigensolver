import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import numpy.matlib
import math, cmath
from scipy import special

# Custom classes
from Config import Config
from MeanFlow import MeanFlow


class Eigensolver:

	def __init__(self):
		
		self.Configure()
		#self.Preliminaries()

		# Mean Flow
		self.ComputeMeanFlow()

		#self.Plots()
		#self.PlotTheMeanFlow()
		self.CreateComplexPlane()
		self.MultRootFindC()
		self.WhichAlpha()
		self.ModeTracer()

	def Configure(self):
		
		self.config 			= Config()

		self.R1					= self.config.R1
		self.R2					= self.config.R2		
		self.dens0 				= self.config.dens0
		self.F 					= self.config.F
		self.E 					= self.config.E
		self.gamma 				= self.config.gamma              
		self.dir_coeff 			= self.config.dir_coeff
		self.invZ1 				= self.config.invZ1
		self.invZ2 				= self.config.invZ2
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
		self.xshift				= self.config.xshift
		self.yshift				= self.config.yshift
		self.xmin 				= self.config.xmin
		self.xmax 				= self.config.xmax
		self.dF					= self.config.dF
		self.x					= self.config.x
		self.DeltaX				= self.config.DeltaX

	

	def ComputeMeanFlow(self):

		MF = MeanFlow()
		self.D0 = MF.D0
		self.U0 = MF.U0
		self.C0 = MF.C0
		self.P0 = MF.P0
		self.Mach0 = MF.Mach0



			
	

	# def ComputeMeanFlowDensity(self):
	# 	N = self.NNxx
	# 	F = self.F
	# 	E = self.E
	# 	R1 = self.R1
	# 	R2 = self.R2
	# 	maxit = self.maxit
	# 	tol = self.tol
	# 	dens0 = self.dens0
	# 	gamma = self.gamma	
	# 	D0 = np.zeros(shape=(N,),dtype=float)
	# 	d0 = np.zeros(shape=(maxit,),dtype=float)
	# 	d0[0] = dens0
		
	# 	k = 0
	# 	while (k < N):
	# 		j = 0
	# 		beta = 0.5*math.pow(F,2)/math.pow((math.pow(R2[k],2)-math.pow(R1[k],2)),2);
	# 		while( j < maxit - 1 ):
			
	# 			f = beta*math.pow(d0[j],-2) + (1/(gamma-1))*math.pow(d0[j],gamma-1) - E
	# 			df = -2*beta*math.pow(d0[j],-3)	+ math.pow(d0[j],gamma-2)

	# 			# Newton's Method	
	# 			d0[j+1] = d0[j] - f/df
					
	# 			a = abs(d0[j+1]-d0[j])

	# 			if(a < tol):
	# 				#print("Newton Algorithm Converged")
	# 				D0[k] = d0[j+1]
	# 				d0 = np.zeros(shape=(maxit,),dtype=float)
	# 				d0[0] = D0[k]
	# 				break
	# 			if(j==maxit-2):
	# 				print("Newton Algorithm Failed to Converge")	
	# 			j = j + 1	
			
	# 		k=k+1
		
	# 	self.D0 = D0
		
		

	# def ComputeMeanFlowAxialVelocity(self):
	# 	N = self.NNxx
	# 	R1 = self.R1
	# 	R2 = self.R2
	# 	F = self.F
	# 	D0 = self.D0
	# 	U0 = np.zeros(shape=(N,),dtype=float)
	# 	U0denominator = np.power(R2,2) - np.power(R1,2)
	# 	U0denominator = U0denominator*np.multiply(D0,U0denominator)
	# 	U0 = np.divide(F,U0denominator)	
	# 	self.U0 = U0	


	# def ComputeMeanFlowSoundSpeed(self):
	# 	N = self.NNxx
	# 	D0 = self.D0
	# 	gamma = self.gamma
	# 	C0 = np.zeros(shape=(N,),dtype=float)
	# 	C0 = np.power(D0,(gamma-1)/2)
	# 	self.C0 = C0
	

	# def ComputeMeanFlowMachNumber(self):
	# 	U0 = self.U0
	# 	C0 = self.C0
	# 	Mach0 = np.divide(U0,C0)
	# 	self.Mach0 = Mach0

	
	# def ComputeMeanFlowPressure(self):
	# 	gamma = self.gamma
	# 	D0 = self.D0
	# 	P0 = np.multiply(1/gamma,np.power(D0,gamma))
	# 	self.P0 = P0	


	def CreateComplexPlane(self):
		mod_z_max = self.mod_z_max
		n = self.n	
		xshift = self.xshift
		yshift = self.yshift

		x = np.linspace(-mod_z_max,mod_z_max,n)	+ xshift
		real = np.matlib.repmat(x,n,1)
		
		y = np.linspace(mod_z_max,-mod_z_max,n)	+ yshift
		imag = np.matlib.repmat(y,n,1)

		imag = np.transpose(imag)
		imag = np.multiply(1j,imag)

		
		comp = np.zeros((n,n),dtype=complex)
		comp = real + imag 

		self.ComplexPlane = comp

	def Newton(self,z,R1,R2,U0,C0):
		m = self.m
		omega = self.omega
		D0 = self.D0[0]
		invZ1 = self.invZ1
		invZ2 = self.invZ2
		F = self.F
		dir_coeff = self.dir_coeff
		#maxit = self.maxit
		maxit = 10
		tol = self.tol
	

		j=0
		while(j<maxit-1):

			f = self.BesEqn(D0,R1,R2,U0,C0,z,dir_coeff)
			df = self.dBesEqn(D0,R1,R2,U0,C0,z,dir_coeff)
			z = z - f/df	

			a = abs(f)
			if(a < tol):
				#print("Newton algorithm converged\n")
				result = [1,z]
				break									
			j=j+1

			if(j == maxit-1):
				#print("Newton algorithm didn't converge\n")
				result = [0,0]
		return result		

		#exit()		

	def ModeTracer(self):
		#MultRootFindc computes the whereabouts of some the coalesed modes (F=0). We take one of these and trace its progress along the inlet for increasing F. This ensures that we're looking at the same mode (one forward running and one backward running) for F \neq 0.   

		# Sets initial value of alpha equal to the initial modal value (i.e. the calculated mode for F=0)

		alpha_init = self.alpha_init 
		dF = self.dF

		alpha = np.zeros(0,dtype=complex)
		alpha =  np.append(alpha, alpha_init) 
		
		# Set the initial mass flow condition	
		F=0 
		j = 1
		
		#count=1;  We have started our array of results alpha with the initial calculated value, and so our 1st point of interest is the second element of this array
		while(j):
			F=F+dF
			break
		# 	if(j==1): 
		# 		z1=alpha(j);
		# 		alpha(j+1)=alpha(j) + inc;
		# 		z2=alpha(j+1);
		# 		test=riddlers(z1, z2, maxIt, tol, R1, R2, m,omega,D0,invZ1,invZ2,F,dir_coeff);
		
		# if(test(2)==1) % Test for convergence
		# 	alpha(count+1)=test(1);
		# 	count=count+1;
		# end
	
	#else
			
# 			test=Newton(z,R1,R2,U0,C0)
	
# 		if(test(2)==1) % Test for convergence
# 			alpha(count+1)=test(1);
# 			count=count+1;
# 		end
		
# 	end
	
# 	if(F>maxF) %Breaks from loop if F lies outside our region of interest
# 		break;
# 	end
# j=j+1;
# end

# if (isempty(alpha))
#     fprintf('no roots found');
#     [alpha]=[];
# end




	def WhichAlpha(self):
		# Specifies which eigenvalue in the alpha_array to track
		alpha = 1
		self.alpha_init = self.alpha_array[alpha]
		print("Tracking Eigenvalue: " + str(self.alpha_init))


	def MultRootFindC(self):
		D0 = self.D0[0]
		R1 = self.R1[0]
		R2 = self.R2[0]
		U0 = 0 # This is to ensure that we first find the F=0 mode
		C0 = self.C0[0]
		
		z = self.ComplexPlane
		dir_coeff = self.dir_coeff

		B = self.BesEqn(D0,R1,R2,U0,C0,z,dir_coeff)

		Br = np.real(B)
		Bi = np.imag(B)

		N = self.n
		alpha = np.zeros(0,dtype=complex)

		count = 0;
		j = 0;

		while( j < N-1 ):     #Loop through rows
			k = 0;
			while( k < N-1 ):  #Loop through columns
				if( Br[j,k] < 0 and Br[j,k+1] > 0 ) or ( Br[j,k] > 0 and Br[j,k+1] < 0 ): 
					z1=z[j,k];
					test = self.Newton(z1,R1,R2,U0,C0)	
					if(test[0] == 1):
						alpha = np.append(alpha,test[1])
				elif( Br[j,k] > 0 and Br[j+1,k] < 0 ) or ( Br[j,k] < 0 and Br[j+1,k] > 0 ): 
					z1=z[j,k];
					test = self.Newton(z1,R1,R2,U0,C0)	
					if(test[0] == 1):
						alpha = np.append(alpha,test[1])
				elif( Bi[j,k] < 0 and Bi[j,k+1] > 0 ) or ( Bi[j,k] > 0 and Bi[j,k+1] < 0 ): 
					z1=z[j,k];
					test = self.Newton(z1,R1,R2,U0,C0)	
					if(test[0] == 1):
						alpha = np.append(alpha,test[1])
				elif( Bi[j,k] > 0 and Bi[j+1,k] < 0 ) or ( Bi[j,k] < 0 and Bi[j+1,k] > 0 ): 
					z1=z[j,k];
					test = self.Newton(z1,R1,R2,U0,C0)	
					if(test[0] == 1):
						alpha = np.append(alpha,test[1])

				k=k+1
			j=j+1

		#print alpha

		alpha = self.RootFilter(alpha)
		print "Cross-Sectional Eigenvalues Found"
		print alpha
		self.alpha_array = alpha

		#k = self.BesEqn(D0,R1,R2,U0,C0,alpha,dir_coeff)
		#print k			

		# if (isempty(alpha)):
		# 	fprintf('no roots found');
		# 	alpha=[];

	def RootFilter(self,alpha):
		tol = self.tol
		N = np.size(alpha)
		j=0
		while ( j < N ):
			k = j+1
			while(k<N):
				diff = abs(alpha[j] - alpha[k])
				if(diff < tol):
					alpha[k] = 0
				k = k+1	
			j = j+1	
		
		t = np.where(alpha != 0)[0]
		alpha1 = alpha[t]
		return alpha1

	def Riddlers(self,z1,z2,R1,R2,U0,C0):
		m = self.m
		omega = self.omega
		D0 = self.D0[0]
		invZ1 = self.invZ1
		invZ2 = self.invZ2
		F = self.F
		dir_coeff = self.dir_coeff
		maxit = self.maxit
		tol = self.tol
		#print(maxit)
		#test = []


		j=1
		while(j<=maxit):
			s = 0
			f1 = self.BesEqn(D0,R1,R2,U0,C0,z1,dir_coeff)
			f2 = self.BesEqn(D0,R1,R2,U0,C0,z2,dir_coeff)

			if(np.real(f1-f2)<0):
				s=-1
			elif(np.real(f1-f2)==0):
				s=0
			else:
				s=1

			z3 = (z1+z2)/2
			f3 = self.BesEqn(D0,R1,R2,U0,C0,z3,dir_coeff)
			# Riddler's Algorithm
			z4 = z3 + ((z3-z1)*s*f3)/(cmath.sqrt(f3*f3-f1*f2)); 

			f4 = self.BesEqn(D0,R1,R2,U0,C0,z4,dir_coeff)
			a =  abs(f4)
			#exit()
			print j
			print a		

			if(a < tol):
				alpha = z4				
				result = 1
				print("Riddler's Algorithm Converged")
				exit()
				break	
				
			if((np.real(f1) < 0 and np.real(f4) < 0) or (np.real(f1) > 0 and np.real(f4) > 0) ):
				z1 = z4
				z2 = z2
			else:
				z1 = z1
				z2 = z4

			if(j == maxit):
				#print("Riddlers did not converge")
				alpha = z4				
				result = 0
				break
				#test = [alpha,result]	

			j = j+1
		return [alpha,result]



	def BesEqn(self,D0,R1,R2,U0,C0,z,dir_coeff):


		maxit = self.maxit
		tol = self.tol
		m = self.m

		# J1 = special.jv(m, z*R1)
		# J1p1 = special.jv(m+1, z*R1)
		# J2 = special.jv(m, z*R2);
		# J2p1 = special.jv(m+1, z*R2)

		# if(R1 != 0):
		# 	Y1 = special.yv(m, z*R1)
		# 	Y1p1 = special.yv(m+1, z*R1)
		
		# Y2 = special.yv(m, z*R2)
		# Y2p1 = special.yv(m+1, z*R2)

		
		# zeta1 = self.zeta1(D0,R1,U0,C0,z,dir_coeff)
		zeta2 = self.zeta2(D0,R2,U0,C0,z,dir_coeff)
		
		if(R1 == 0):
			#y = -z*R2*J2p1 + (m-zeta2)*J2
			y = z*R2*special.jvp(m,z*R2) - zeta2*special.jv(m,z*R2)

		else:
			ratio_top = (-z*R2*Y2p1 + (m-zeta2)*Y2)
			ratio_bottom = (-z*R1*Y1p1 + (m+zeta1)*Y1)  # Need to double-check this line
			ratio=ratio_top/ratio_bottom
			y = -z*R2*J2p1 + (m-zeta2)*J2 - ratio*(-z*R1*J1p1 + (m+zeta1)*J1)

		return y	

	def dBesEqn(self,D0,R1,R2,U0,C0,z,dir_coeff):
		maxit = self.maxit
		tol = self.tol
		#gamma = self.gamma
		m = self.m		

		J1 = special.jv(m, z*R1)
		J1p1 = special.jv(m+1, z*R1)
		J1p2 = special.jv(m+2, z*R1)

		J2 = special.jv(m, z*R2)
		J2p1 = special.jv(m+1, z*R2)
		J2p2 = special.jv(m+2, z*R2)

		Y1 = special.yv(m, z*R1)
		Y1p1 = special.yv(m+1,z*R1)
		Y1p2 = special.yv(m+2,z*R1)

		Y2 = special.yv(m, z*R2)
		Y2p1 = special.yv(m+1, z*R2)
		Y2p2 = special.yv(m+2,z*R2)

		if( R1 !=0 ):
			dJ1 = (m/(z*R1))*J1 - J1p1
			dY1 = (m/(z*R1))*Y1 - Y1p1
			dJ1p1 = ((m+1)/(z*R1))*J1_p1 - J1p2
			dY1p1 = ((m+1)/(z*R1))*Y1p1 - Y1p2
			zeta1 = zeta1(D0,R1,U0,C0,z,dir_coeff)
			dzeta1 = dZeta1(D0,R1,U0,C0,z,dir_coeff)


		dJ2 = (m/(z*R2))*J2 - J2p1
		dY2 = (m/(z*R2))*Y2 - Y2p1		
		dJ2p1 = ((m+1)/(z*R2))*J2p1 - J2p2
		dY2p1 = ((m+1)/(z*R2))*Y2p1 - Y2p2

		zeta2 = self.zeta2(D0,R2,U0,C0,z,dir_coeff)
		dzeta2 = self.dZeta2(D0,R2,U0,C0,z,dir_coeff)

		if( R1==0 ):

			j2 = (-z*R2*J2p1 + (m-zeta2)*J2)
			dj2 = -R2*J2p1-z*math.pow(R2,2)*dJ2p1 + (m-zeta2)*R2*dJ2-dzeta2*J2
			dy = dj2; 

		else:

			j1 = (-z*R1*J1p1 + (m+zeta1)*J1);  
			j2 = (-z*R2*J2_p1 + (m-zeta2)*J2);
			dj1 = -R1*J1p1-z*math.pow(R1,2)*dJ1p1 + (m+zeta1)*R1*dJ1+dzeta1*J1
			dj2 = -R2*J2p1-z*math.pow(R2,2)*dJ2p1 + (m-zeta2)*R2*dJ2-dzeta2*J2

			y1 = (-z*R1*Y1p1 + (m+zeta1)*Y1) 
			y2 = (-z*R2*Y2_p1 + (m-zeta2)*Y2)
			dy1 = -R1*Y1p1-z*math.pow(R1,2)*dY1_p1 + (m+zeta1)*R1*dY1+dzeta1*Y1
			dy2 = -R2*Y2p1-z*math.pow(R2,2)*dY2_p1 + (m-zeta2)*R2*dY2-dzeta2*Y2

			a = y2/y1
			b = dy2/y1
			c = y2*dy1/(y1*y1);
			dy = (dj2-b*j1 + c*j1-a*dj1)

		return dy




			

	# Need to check that this all works when dir_coeff = -1
	
	def sigma(self,U0,C0,alpha,dir_coeff):
		omega = self.omega
		alpha_squared = np.power(alpha,2)
		k = (math.pow(C0,2) - math.pow(U0,2))/math.pow(omega,2)
		k_alpha = k*alpha_squared
		sigma = dir_coeff*np.sqrt(1-k_alpha)
		return sigma

	def mu(self,U0,C0,alpha,dir_coeff):
		omega=self.omega
		sigma = self.sigma(U0,C0,alpha,dir_coeff)
		k = omega/(math.pow(C0,2) - math.pow(U0,2))
		mu = k*(C0*sigma - U0)	
		return mu

	def zeta1(self,D0,R1,U0,C0,alpha,dir_coeff):
		omega = self.omega
		invZ1 = self.invZ1
		mu = self.mu(U0,C0,alpha,dir_coeff)
		zeta1 = (-1j*D0*R1*invZ1/omega)*np.power((omega - mu*U0),2)
		return zeta1

	def zeta2(self,D0,R2,U0,C0,alpha,dir_coeff):
		omega = self.omega
		invZ2 = self.invZ2
		mu = self.mu(U0,C0,alpha,dir_coeff)
		zeta2 = (-1j*D0*R2*invZ2/omega)*np.power((omega - mu*U0),2)
		return zeta2

	def dZeta1(self,D0,R1,U0,C0,z,dir_coeff):
		omega = self.omega
		invZ1 = self.invZ1
		numerator = -1j*C0*D0*R1*U0*z*invZ1
		denominator = omega*cmath.sqrt(-C0*C0*z*z + U0*U0*z*z + omega*omega)
		dZeta1 = dir_coeff*numerator/denominator
		return dZeta1

	def dZeta2(self,D0,R2,U0,C0,z,dir_coeff):
		omega = self.omega
		invZ2 = self.invZ2
		numerator = -1j*C0*D0*R2*U0*z*invZ2
		denominator = omega*cmath.sqrt(-C0*C0*z*z + U0*U0*z*z + omega*omega)
		dZeta2 = dir_coeff*numerator/denominator
		return dZeta2			



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

	def PlotTheMeanFlow(self):
		x = self.x
		U0 = self.U0
		C0 = self.C0
		Mach0 = self.Mach0
		fig, ax = plt.subplots()
		ax.plot(x, U0,linewidth=3)	
		ax.plot(x, C0,linewidth=3)	
		ax.plot(x, Mach0,linewidth=3)
		ax.set(xlabel='$x$',title='Mean Flow')	

		fig.savefig("MeanFlow.png")
		plt.show()

	def PlotContour(self):
		x = np.real(z[0,:])
		y = np.imag(z[:,0])
		X, Y = np.meshgrid(x,y)
		absB = abs(B)
		#levels = [0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1]	
		levels = np.linspace(0,1,10)

		plt.contour(X,Y,absB,levels)
		plt.show()


		