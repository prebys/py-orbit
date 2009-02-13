#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------
addr="./ext/laserstripping/"
import sys
sys.path.append(addr+"/PyModules/")
import math

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

import os
import orbit_mpi
import random
import part_generator

print "Start."





z0=2.0e-3
Bx=1./10000.
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
emtX = 0.225e-6/(0.736*0.736) # [m*rad]

alphaY = 1.979            # [rad]
betaY = 3.906             # [m]
emtY = 3.803e-6/(3.5*3.5) # [m*rad]

relativeSpread = 1.0e-4

dispD = 0.    # [m]
dispDP = 2.58 # [rad]



cutOffX = math.sqrt(emtX*betaX)*3.0
trGenX = part_generator.TransverseCoordGen(alphaX,betaX,emtX,cutOffX)

cutOffY = math.sqrt(emtY*betaY)*3.0
trGenY = part_generator.TransverseCoordGen(alphaY,betaY,emtY,cutOffY)

#print "cutOffX= %6.2f cutOffY= %6.2f "%(cutOffX*1.0e+3,cutOffY*1.0e+3)

pGen = part_generator.EnergyGen(TK,relativeSpread)


partGen = part_generator.ParticlesGen(dispD,dispDP,trGenX,trGenY,pGen)


bunch = Bunch()
bunch_target = Bunch()
bunch.charge(0)
bunch_target.charge(0)
bunch.addPartAttr("Amplitudes")
for i in range(N_part):
    (x,px,y,py,z,pz) = partGen.getCoords()
    bunch.addParticle(x,px,y,py,-z0,pz)
    bunch.partAttrValue("Amplitudes",i,1,1.0)

#two bunches should have the same particle attributes to copy operation


E = bunch.mass() + TK
P = math.sqrt(E*E - bunch.mass()*bunch.mass())
vz=299792458*P/E
time_step = (2*z0/vz)/n_step



Stark=HydrogenStarkParam(addr+"/transitions/",n_states)

#print alpha
#print kz
#sys.exit(1)


for power,wx,fxy in [(power,wx,fxy) 
                    for power in [1.0e6,1.5e6,2.0e6] 
                    for wx in [60.0e-6 ,80.0e-6,100.0e-6,120.0e-6,140.0e-6,160.0e-6,180.0e-6,200.0e-6] 
                    for fxy in [-0.00,-0.02,-0.04,-0.06,-0.08,-0.1,-0.12,-0.14,-0.16,-0.18,-0.2]]:

	
	bunch_target.deleteAllParticles()
	bunch.copyBunchTo(bunch_target)

	la0= 2*math.pi*5.291772108e-11/7.297352570e-3/(
		Stark.getStarkEnergy(bunch_target.mass(),   0,2,0,  0.,0.,0.,  Bx,0.,0.,   0.,0.,P)
		-Stark.getStarkEnergy(bunch_target.mass(),   0,0,0,  0.,0.,0.,  Bx,0.,0.,   0.,0.,P))
	
	kz=-1/math.sqrt(math.pow(P/(bunch_target.mass()*(la/la0-1)-TK),2)-1)
    
#    alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi
#    print alpha

	fx=fy=fxy
	LFS=HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,la) 
	LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   1.,0.,1./kz,  0.,1.,0.)

	tracker = RungeKuttaTracker(1000.0)
	First = LasStripExternalEffects(LFS,Stark,10000.)
	fS=LSFieldSource(0.,0.,0.,Bx,0.,0.)
	tracker.track(bunch_target,0,time_step*n_step, time_step,fS,First)
	
	population = 0.
	for i in range(bunch_target.getSize()):
		population += (1-bunch_target.partAttrValue("Amplitudes",i,1))
	population = population/bunch_target.getSize()
	#file=open("020_tuning.txt",'a')
	#file.write('%f'%Bx+"\t")
	#file.write('%f'%population+"\n")
	#file.close() 
	bunch_target.deleteAllParticles()
	print  "W = %4.1f  wx [um] = %4.1f  dist[cm]= %4.1f  Population: %7.3f "%(power,1.0e+6*wx,fxy*100,population)


