#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys
sys.path.append("/home/tg4/workspace/PyOrbit/ext/laserstripping/PyModules/")
import math

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

import os
import orbit_mpi
import random
import part_generator

from orbit_mpi import mpi_comm,mpi_datatype,mpi_op

rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)

if(rank == 0): print "Start."

random.seed((rank+1)*1257)

z0=2.0e-3
Bx=0.1/10000.
n_states=3
TK = 1.0
#fx=fy=-0.1
la=355.0e-9

#wx=59.0e-6 
wy=850.0e-6
#power=5000000.
n_step = 50000
N_part = 1




alphaX = -3.294               # [rad]
betaX = 26.7                  # [m]
emtX = 0.225e-6               # [m*rad]

alphaY = 1.979            # [rad]
betaY = 3.906             # [m]
emtY = 0.380e-6           # [m*rad]

relativeSpread = 1.0e-4

dispD = 0.    # [m]
dispDP = 2.58 # [rad]



cutOffX = math.sqrt(emtX*betaX)*3.0
trGenX = part_generator.TransverseCoordGen(alphaX,betaX,emtX,cutOffX)

cutOffY = math.sqrt(emtY*betaY)*3.0
trGenY = part_generator.TransverseCoordGen(alphaY,betaY,emtY,cutOffY)

if(rank == 0): print "cutOffX= %6.2f cutOffY= %6.2f "%(cutOffX*1.0e+3,cutOffY*1.0e+3)

pGen = part_generator.EnergyGen(TK,relativeSpread)

partGen = part_generator.ParticlesGen(dispD,dispDP,trGenX,trGenY,pGen)


bunch_target = Bunch()
bunch_target.charge(0)
bunch_target.mass(0.938256 + 0.000511)
bunch_target.addPartAttr("Amplitudes")

E = bunch_target.mass() + TK
P = math.sqrt(E*E - bunch_target.mass()*bunch_target.mass())
vz=299792458*P/E
time_step = (2*z0/vz)/n_step




av=0;
nstat=100000
for i in range(nstat):
    (x,px,y,py,z,pz) = partGen.getCoords()
    av+=pz
    
av=av/nstat


av2=0
for i in range(nstat):
    (x,px,y,py,z,pz) = partGen.getCoords()
    av2+=px*px
    
av2=math.sqrt(av2/nstat)

print av2/av
sys.exit(1)




Stark=HydrogenStarkParam("/home/tg4/workspace/PyOrbit/ext/laserstripping/transitions/",n_states)

#print alpha
#print kz
#sys.exit(1)

if(rank == 0):
	f_out = open("/home/tg4/workspace/PyOrbit/result_ls.dat","w")
	f_out.write("cpu_time   W [MW]     wx [um]   dist[cm]   Population \n")
	f_out.close()
	print "======START LOOP=========="


time_start = orbit_mpi.MPI_Wtime()

for power,wx,fxy in [(power,wx,fxy) 
                    for power in [1.0e6,1.5e6,2.0e6] 
                    for wx in [60.0e-6 ,80.0e-6,100.0e-6,120.0e-6,140.0e-6,160.0e-6,180.0e-6,200.0e-6] 
                    for fxy in [-0.00,-0.02,-0.04,-0.06,-0.08,-0.1,-0.12,-0.14,-0.16,-0.18,-0.2]]:
	population = 0.
	for i in range(N_part):
		bunch_target = Bunch()
		bunch_target.charge(0)
		bunch_target.mass(0.938256 + 0.000511)
		bunch_target.addPartAttr("Amplitudes")		
		(x,px,y,py,z,pz) = partGen.getCoords()
		bunch_target.addParticle(0.,0.,0.,0.,-z0,pGen.getP0())
		#bunch_target.addParticle(x,px,y,py,z,pz)
		bunch_target.partAttrValue("Amplitudes",0,1,1.0)

		la0= 2*math.pi*5.291772108e-11/7.297352570e-3/(Stark.getStarkEnergy(bunch_target.mass(),   0,2,0,  0.,0.,0.,  Bx,0.,0.,   0.,0.,P) -Stark.getStarkEnergy(bunch_target.mass(),   0,0,0,  0.,0.,0.,  Bx,0.,0.,   0.,0.,P))
	
		kz=-1/math.sqrt(math.pow(P/(bunch_target.mass()*(la/la0-1)-TK),2)-1)
    
		#       alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi
		#       print alpha
		fx=fy=fxy
		LFS=HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,la) 
		LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   1.,0.,1./kz,  0.,1.,0.)

		tracker = RungeKuttaTracker(1000.0)
		First = LasStripExternalEffects(LFS,Stark,10000.)
		fS=LSFieldSource(0.,0.,0.,Bx,0.,0.)
		print "debug  ===track ========== rank=",rank," starkEnrg=",Stark.getStarkEnergy(bunch_target.mass(),   0,2,0,  0.,0.,0.,  Bx,0.,0.,   0.,0.,P)
		tracker.track(bunch_target,0,time_step*n_step, time_step,fS,First)
		print "debug  amp=",(1-bunch_target.partAttrValue("Amplitudes",0,1))
		population += (1-bunch_target.partAttrValue("Amplitudes",0,1))

	population = population/N_part 	
	mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
	op = mpi_op.MPI_SUM
	data_type = mpi_datatype.MPI_DOUBLE
	population = orbit_mpi.MPI_Allreduce(population,data_type,op,mpi_comm.MPI_COMM_WORLD)	
	population = population/mpi_size
	time_live = orbit_mpi.MPI_Wtime() - time_start
	res = " %6.0f  %4.1f  %4.1f  %4.1f    %7.3f "%(time_live,power/1.0e+6,1.0e+6*wx,fxy*100,population)
	if(rank == 0): print  "W [MW]= %4.1f  wx [um] = %4.1f  dist[cm]= %4.1f  Population: %7.3f "%(power/1.0e+6,1.0e+6*wx,fxy*100,population)
	if(rank == 0):
		f_out = open("/home/tg4/workspace/PyOrbit/result_ls.dat","a")
		f_out.write(res + "\n")
		f_out.close()

