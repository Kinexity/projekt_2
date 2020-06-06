import numpy as np
import scipy
import math as m
import numpy as np
from scipy import integrate
import random
import multiprocessing

number_of_asteroids = 3

Y = np.zeros(shape=(number_of_asteroids,4))

asteroid_sun_distance=5
mars_orbiting_speed=24

#24 km/s
print(Y)

def vec_function(x,y):
	return np.sqrt(np.power(x,2)-np.power(y,2))

def cylindrical_to_cartesian_transformation(r,theta,Z):
	Z[:,0]=r*np.cos(theta)
	Z[:,2]=r*np.sin(theta)
	return Z

starting_thetas=[random.uniform(0,2*m.pi) for __ in range(number_of_asteroids)]
'''print(cylindrical_to_cartesian_transformation(asteroid_sun_distance,starting_thetas))'''

#liczymy macierz obrotu od wektora mars_orbiting_speed (wikipedia)
starting_velocities_r = [-mars_orbiting_speed*m.sin(random.uniform(0,m.pi/2)) for __ in range(number_of_asteroids)]
starting_velocities_theta = vec_function(mars_orbiting_speed ,starting_velocities_r)
'''print(starting_velocities_theta)'''

def make_starting_array(r,theta,funct):
	Y[:,1] = starting_velocities_r
	Y[:,3] = starting_velocities_theta
	return funct(r,theta,Y)
print(make_starting_array(asteroid_sun_distance,starting_thetas,cylindrical_to_cartesian_transformation))
print(Y[0,1])

sun_sgp=4
earth_sgp=3
moon_sqp=2
sun_coordinates=[0,0]
earth_sun_distance=3
step=0.01
liminal_time=25.0
radial_earth_speed=1
def earth_coordinates(theta,funct):
	return funct(earth_sun_distance,theta,B)

B = np.zeros(shape=(int(liminal_time/step),3))
B[:,0]= np.arange(0,liminal_time,step)
#wspolrzedne x-owe
'''B[:,1]= earth_sun_distance*np.cos(radial_speed*B[:,0])
B[:,2]= earth_sun_distance*np.sin(radial_speed*B[:,0])'''
print(B)



'''class asteroid_data:
	def '''

file = open("dane.txt","w+")
'''data.csv'''
file.truncate(0)
file.close()


def solvr(Y,t):
    return[-Y[1],sun_sgp/(np.hypot(Y[0], Y[2])**3)*Y[0],-Y[3],sun_sgp/(np.hypot(Y[0], Y[2])**3)*Y[2]]

def main(i):
    a_t = np.arange(0, liminal_time, step)


    asol = integrate.odeint(solvr, [Y[i,0],Y[i,1],Y[i,2],Y[i,3]], a_t)
    astack = np.c_[a_t, asol[:,0], asol[:, 1],asol[:,2], asol[:, 3]]
    np.savetxt('dane.txt', astack, delimiter=',', header='t, dx,vx,dy,vy', comments='')

processes = []
i=0
for _ in range(number_of_asteroids):
	p= multiprocessing.Process(target=main(i))
	p.start()
	i+=1
	processes.append(p)
