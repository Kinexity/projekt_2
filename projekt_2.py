import numpy as np
import scipy
import math as m
import numpy as np
from scipy import integrate
import random
import multiprocessing
import functools

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
#losowo wybrane wspolrzedne asteroid na orbicie marsa

#liczymy macierz obrotu od wektora mars_orbiting_speed
starting_velocities_r = [-mars_orbiting_speed*m.sin(random.uniform(0,m.pi/2)) for __ in range(number_of_asteroids)]
starting_velocities_theta = vec_function(mars_orbiting_speed ,starting_velocities_r)
'''print(starting_velocities_theta)'''

#wektor predkosci=mars_orbiting_speed odchylam o randomowy kat
def make_starting_array(r,theta,funct):
	Y[:,1] = starting_velocities_r
	Y[:,3] = starting_velocities_theta
	return funct(r,theta,Y)
print(make_starting_array(asteroid_sun_distance,starting_thetas,cylindrical_to_cartesian_transformation))
print(Y[0,1])

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
def earth_coordinates(theta,funct):
	return funct(earth_sun_distance,theta,B)

#tu będą licznoe wspolrzedne ziemi i ksiezyca
B = C = np.zeros(shape=(int(liminal_time/step),3))
B[:,0] = C[:,0] = np.arange(0,liminal_time,step)
#wsp x-owe ziemi
B[:,1]= earth_sun_distance*np.cos(radial_earth_speed*B[:,0])
B[:,2]= earth_sun_distance*np.sin(radial_earth_speed*B[:,0])
print(B)
#wsp x-owe ks
C[:,1]= B[:,1] + earth_moon_distance*np.cos(radial_moon_speed*B[:,0])
C[:,2]= B[:,2] + earth_moon_distance*np.sin(radial_moon_speed*B[:,0])
print(C)

'''def solve(Y,t):
    return[Y[1],-earth_sgp/(np.hypot(Y[0], Y[2])**3)*Y[0],Y[3],-earth_sgp/(np.hypot(Y[0], Y[2])**3)*Y[2]]'''

'''class asteroid_data:
	def '''

file = open("data.csv","w+")

file.truncate(0)
file.close()
Ax_earth_speed = -radial_earth_speed*earth_sun_distance*np.sin(radial_earth_speed*B[:,0])
Ay_earth_speed = radial_earth_speed*earth_sun_distance*np.cos(radial_earth_speed*B[:,0])
#rownanie rozniczkowe na przyciaganie od  slonca i ziemi- tor ruchu asteroidy, gdzieś w nim jest blad bo nie liczy dobrze
def solvr(P,t):
    return[P[1],-sun_sgp/(np.hypot(P[0], P[2])**3)*P[0]-earth_sgp/(np.hypot(P[0]-P[4], P[2]-P[5])**3)*(P[0]-P[4]) ,P[3],-sun_sgp/(np.hypot(P[0], P[2])**3)*P[2]-earth_sgp/(np.hypot(P[0]-P[4], P[2]-P[5])**3)*(P[2]-P[5]),-radial_earth_speed*P[5],radial_earth_speed*P[4]]


def main(i):
    a_t = np.arange(0, liminal_time, step)
    asol = integrate.odeint(solvr, [Y[i,0],Y[i,1],Y[i,2],Y[i,3],earth_sun_distance,0], a_t)
    astack = np.c_[a_t, asol[:,0], asol[:, 1],asol[:,2], asol[:, 3]]
    np.savetxt('data.csv', astack, delimiter=',', header='time, asteroid_x_from_sun,asteroid_velocity_x,asteroid_y_from_sun,asteroid_velocity_y', comments='')
    

'''def main(i):
    a_t = np.arange(0, liminal_time, step)


    asol = integrate.odeint(solvr, [Y[i,0],Y[i,1],Y[i,2],Y[i,3]], a_t)
    astack = np.c_[a_t, asol[:,0], asol[:, 1],asol[:,2], asol[:, 3]]
    np.savetxt('dane.txt', astack, delimiter=',', header='t, dx,vx,dy,vy', comments='')'''

processes = []
i=0
#rozwiązuje w pętli z multiprocessig

    
    
'''for _ in range(number_of_asteroids):
	p= multiprocessing.Process(target=main(i))
	p.start()
	i+=1
	processes.append(p)'''
main(1)



#ostatecznie zapisuje do pliku main(1) nie wiedziałem jak zrobic zeby dopisywal w porzadku po multiprocessing i zeby dopisywal arraye a nie tworzyl caly plik od nowa po kazdej serii

