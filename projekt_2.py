import numpy as np
import scipy
import math as m
import numpy as np
from scipy import integrate
import random
import concurrent.futures
import functools
import multiprocessing
import matplotlib
import matplotlib.pyplot as plt

#input parameters
 
'''number_of_asteroids = 1000

#standard gravitational parameter
sun_sgp = 13.34*pow(10,10)
#[km^3/s^3]

earth_sgp = 4*pow(10,5)

moon_sqp = 5*pow(10,3) 

sun_coordinates = [0,0]

earth_sun_distance = 1.5*pow(10,8)
#[km]
earth_moon_distance = 4*pow(10,5)

step= 3600
#[s]
liminal_time = 15768*pow(10,4)
#[s]
#estimated 5 years

radial_earth_speed = 2*pow(10,-7)
#[radiants/s]

radial_moon_speed = 2.6*pow(10,-6)

moon_radius = 1.7* pow(10,3)

earth_radius = 6.3*pow(10,3)

sun_radius = 7*pow(10,5)

asteroid_sun_distance = 2.28*pow(10,8)

mars_orbiting_speed = 24
#[km/s]'''
number_of_asteroids = 100
sun_sgp=4*10**3
earth_sgp=3
moon_sqp=2
sun_coordinates=[0,0]
earth_sun_distance=3
earth_moon_distance=0.5
step=0.01
liminal_time=25.0
radial_earth_speed=1
radial_moon_speed=12
moon_radius=0.052
earth_radius = 0.07
sun_radius = 0.1
asteroid_sun_distance=5
mars_orbiting_speed=24

#defining useful vector functions
def angle_between(x1, y1, x2, y2):
    dot = x1*x2 + y1*y2
    det = x1*y2 - y1*x2
    return m.atan2(det, dot)

def vec_function(x,y):
	return np.sqrt(np.power(x,2)-np.power(y,2))

def cylindrical_to_cartesian_transformation(r,theta,Z):
	Z[:,0]=r*np.cos(theta)
	Z[:,2]=r*np.sin(theta)
	return Z

Y = np.zeros(shape=(number_of_asteroids,4))

#random thetas --> angles between vector tangent to Mars orbit and asteroid velocity vector
starting_thetas=[random.uniform(0,2*m.pi) for __ in range(number_of_asteroids)]

starting_velocities_r = [-mars_orbiting_speed*m.sin(random.uniform(0,m.pi/2)) for __ in range(number_of_asteroids)]
starting_velocities_theta = vec_function(mars_orbiting_speed ,starting_velocities_r)

#put starting random asteroid data to array Y
def make_starting_array(r,theta,funct):
	Y[:,1] = starting_velocities_r
	Y[:,3] = starting_velocities_theta
	return funct(r,theta,Y)
make_starting_array(asteroid_sun_distance,starting_thetas,cylindrical_to_cartesian_transformation)

#funtion for earth motion
def earth_coordinates(theta,funct):
	return funct(earth_sun_distance,theta,B)

#put earth data in next time steps to array
B = C = np.zeros(shape=(int(liminal_time/step),3))
B[:,0] = C[:,0] = np.arange(0,liminal_time,step)

#x/y earth coordinates in next time steps
B[:,1]= earth_sun_distance*np.cos(radial_earth_speed*B[:,0])
B[:,2]= earth_sun_distance*np.sin(radial_earth_speed*B[:,0])

#earth x/y coordinates- mooon x/y coordinates
C[:,1]=  earth_moon_distance*np.cos(radial_moon_speed*B[:,0])
C[:,2]=  earth_moon_distance*np.sin(radial_moon_speed*B[:,0])

#justification: memory used for storing B and C is small, they are useful later when we estimate when asteroid crushes into moon

'''file = open("dane.txt","w+")

file.truncate(0)
file.close()'''

Ax_earth_speed = -radial_earth_speed*earth_sun_distance*np.sin(radial_earth_speed*B[:,0])
Ay_earth_speed = radial_earth_speed*earth_sun_distance*np.cos(radial_earth_speed*B[:,0])
#earth velocity

def solvr(P,t):
    return[P[1],-sun_sgp/(np.hypot(P[0], P[2])**3)*P[0]-earth_sgp/(np.hypot(P[0]-P[4], P[2]-P[5])**3)*(P[0]-P[4]) ,P[3],-sun_sgp/(np.hypot(P[0], P[2])**3)*P[2]-earth_sgp/(np.hypot(P[0]-P[4], P[2]-P[5])**3)*(P[2]-P[5]),-radial_earth_speed*P[5],radial_earth_speed*P[4]]
#differential equation for asteroid motion, with gravity force from earth and sun
#justification: we ignored force from moon for code clarity 


def main(i):
    a_t = np.arange(0, liminal_time, step)
    asol = integrate.odeint(solvr, [Y[i,0],Y[i,1],Y[i,2],Y[i,3],earth_sun_distance,0], a_t)
    astack = np.c_[a_t, asol[:,0], asol[:,2]]
    for index in range(0, int(liminal_time/step)):
        if m.hypot(asol[index,0]-C[index,1]-B[index,1],asol[index,2]-C[index,2]-B[index,2])< moon_radius:
            angle = angle_between(-C[index,1],-C[index,2], asol[index,0]-C[index,1]-B[index,1],asol[index,2]-C[index,2]-B[index,2])
            return angle
            
# returns angle between vectors: from earth centre to sun centre and from moon centre to asteroid collision point in range (-pi,pi)
        if m.hypot(asol[index,0]-B[index,1],asol[index,2]-B[index,2]) < earth_radius or m.hypot(asol[index,0],asol[index,2])< sun_radius:
            return None

def plotData(plotInput):
    # Data for plotting

    zakresy = (-3.14,3.14)
    numbin = 20
    plt.hist(plotInput, numbin , zakresy)
    plt.xlabel("Angle in radians")
    plt.ylabel("Quantity")
    plt.show()


try:
    # collisions to graph
    plotInput = []
    bindex = np.arange(0,number_of_asteroids, 1)
    with concurrent.futures.ProcessPoolExecutor(multiprocessing.cpu_count()) as executor:
        for number, angle in zip(bindex, executor.map(main, bindex)):
            print('asteroid %d angle %s' % (number, angle))
            if angle != None:
                plotInput.append(angle)
    plotData(plotInput)

except ODEintWarning as error:
    pass
