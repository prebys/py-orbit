#-----------------------------------------------------
#Creates Grid2D for charge density and Poisson Solver
#for a charged string (2D case) 
#Copmares with the exact result.
#-----------------------------------------------------
#The one very useful property of the FFT Poisson solver - scalability 
#you can change the absolute value of steps X and Y and you will get 
#the right results as soon you keep stepX/stepY constant.
#

import sys
import math

import orbit_mpi

from spacecharge import Grid2D
from spacecharge import PoissonSolverFFT2D

print "Start."

sizeX = 200
sizeY = 200
xMin = -5.0
xMax = +5.0
yMin = -5.0
yMax = +5.0

scale_coeff =3.0
solver = PoissonSolverFFT2D(sizeX,sizeY,xMin/scale_coeff,xMax/scale_coeff,yMin/scale_coeff,yMax/scale_coeff)

gridRho = Grid2D(sizeX,sizeY,xMin,xMax,yMin,yMax)
gridPhi = Grid2D(sizeX,sizeY,xMin,xMax,yMin,yMax)

chrage_pos_x = 2.5
chrage_pos_y = 0.0
charge = 1.0
gridRho.binValue(charge,chrage_pos_x,chrage_pos_y)

solver.findPotential(gridRho,gridPhi)

r_test = 4.0
n_angle_steps = 10
angle_step = 360./(n_angle_steps - 1)
print "  i    x       y         phi         phi_theory    ratio phi/theory  "
for i in xrange(n_angle_steps):
	angle = math.pi*i*angle_step/180.
	x = r_test*math.cos(angle)
	y = r_test*math.sin(angle)
	phi = gridPhi.getValue(x,y)
	dist = (chrage_pos_x - x)*(chrage_pos_x - x) + (chrage_pos_y - y)*(chrage_pos_y - y)
	dist = math.sqrt(dist)
	phi_th = -math.log(dist)
	print "",i," %7.4f  %7.4f  %12.5g  %12.5g  %12.5g  "%(x,y,phi,phi_th,(phi/phi_th)) 

print "Stop."

