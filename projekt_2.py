import numpy as np
import scipy
import math as m
import numpy as np
from scipy import integrate
import random
import concurrent.futures
import functools

number_of_asteroids = 100


Y = np.zeros(shape=(number_of_asteroids,4))

asteroid_sun_distance=5
mars_orbiting_speed=24

#24 km/s
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

starting_thetas=[random.uniform(0,2*m.pi) for __ in range(number_of_asteroids)]

#losowo wybrane wspolrzedne asteroid na orbicie marsa

#liczymy macierz obrotu od wektora mars_orbiting_speed
starting_velocities_r = [-mars_orbiting_speed*m.sin(random.uniform(0,m.pi/2)) for __ in range(number_of_asteroids)]
starting_velocities_theta = vec_function(mars_orbiting_speed ,starting_velocities_r)


#wektor predkosci=mars_orbiting_speed odchylam o randomowy kat
def make_starting_array(r,theta,funct):
	Y[:,1] = starting_velocities_r
	Y[:,3] = starting_velocities_theta
	return funct(r,theta,Y)
print(make_starting_array(asteroid_sun_distance,starting_thetas,cylindrical_to_cartesian_transformation))


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
def earth_coordinates(theta,funct):
	return funct(earth_sun_distance,theta,B)

#tu będą licznoe wspolrzedne ziemi i ksiezyca
B = C = np.zeros(shape=(int(liminal_time/step),3))
B[:,0] = C[:,0] = np.arange(0,liminal_time,step)
#wsp x-owe ziemi
B[:,1]= earth_sun_distance*np.cos(radial_earth_speed*B[:,0])
B[:,2]= earth_sun_distance*np.sin(radial_earth_speed*B[:,0])

#wsp x/y-owe ks - wsp x/y-owe ziemi
C[:,1]=  earth_moon_distance*np.cos(radial_moon_speed*B[:,0])
C[:,2]=  earth_moon_distance*np.sin(radial_moon_speed*B[:,0])


file = open("dane.txt","w+")

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
    astack = np.c_[a_t, asol[:,0], asol[:,2]]
    for index in range(0, int(liminal_time/step)):
        if m.hypot(asol[index,0]-C[index,1]-B[index,1],asol[index,2]-C[index,2]-B[index,2])< moon_radius:
            angle = angle_between(-C[index,1],-C[index,2], asol[index,0]-C[index,1]-B[index,1],asol[index,2]-C[index,2]-B[index,2])
            print(angle)
# wywala kąt między wektorem łąćzącym  środek ks i punkt zderzenia oraz wektorem od sr ks do ziemi
        if m.hypot(asol[index,0]-B[index,1],asol[index,2]-B[index,2]) < earth_radius or m.hypot(asol[index,0],asol[index,2])< sun_radius:
            break

#mozna tez wrzucic map dla lepszego porzadku ale to tylko zwalnia nie wplywa na wyniki

try:
    with concurrent.futures.ProcessPoolExecutor() as executor:
        bindex = np.arange(0,number_of_asteroids, 1)
        results = [executor.submit(main,number) for number in bindex]
except ODEintWarning as error:
    pass
     

