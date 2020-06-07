import numpy as np
import scipy
import math as m
from scipy import integrate
import random
import multiprocessing
import functools
import operator
from threading import Lock
import astropy.constants as astro_consts
from concurrent.futures import ThreadPoolExecutor

earth_sgp = astro_consts.GM_earth
earth_orbital_period = 31557600. #wystarczająco dobre przybliżenie
earth_orbit_radius = 1.496e+8
earth_angular_velocity = 2 * m.pi / earth_orbital_period
moon_sgp = 4.9048696e+12
moon_orbital_period = 27.321661 * 86400
moon_orbit_radius = 1.496e+8
moon_angular_velocity = 2 * m.pi / earth_orbital_period
sun_sgp = astro_consts.GM_sun

def cyl_to_cart(ro, fi):
	return [ro * np.cos(fi), ro * np.sin(fi)]

def earth_pos(t):
	return cyl_to_cart(earth_orbit_radius, earth_angular_velocity * t)

def moon_pos_angle_randomize():
	fi = 2 * m.pi * random.random()
	return lambda t: list(map(operator.add, cyl_to_cart(moon_orbit_radius, moon_angular_velocity * t + fi), earth_pos(t)))

number_of_asteroids = 1000

def angle_between(x1, y1, x2, y2):
	dot = x1 * x2 + y1 * y2      
	det = x1 * y2 - y1 * x2     
	return m.atan2(det, dot)

#rownanie rozniczkowe na przyciaganie od slonca i ziemi- tor ruchu asteroidy,
#gdzieś w nim jest blad bo nie liczy dobrze
def asteroid_position(P,t):
	earth_current_pos = earth_pos(t)
	sun_distance_3 = np.hypot(P[0], P[2]) ** 3
	earth_distance_3 = np.hypot(P[0] - earth_current_pos[0], P[2] - earth_current_pos[1]) ** 3
	return[P[1],
		   -sun_sgp / sun_distance_3 * P[0] - earth_sgp / earth_distance_3 * (P[0] - earth_current_pos[0]),
		   P[3],
		   -sun_sgp / sun_distance_3 * P[2] - earth_sgp / earth_distance_3 * (P[2] - earth_current_pos[1])]

def place_asteroid():
	return [] # będzie zwracać listę z wektorami prędkości i położenia

angles_list = []
lock = Lock()

def asteroid_batch(thr_index):
	while (len(angles_list) < 10000):
		local_angles_list = []
		for i in range(number_of_asteroids):
			solver = lambda t : integrate.odeint(asteroid_position, place_asteroid(), t)
		lock.acquire()
		print("<---------->")
		print("Wątek: ", thr_index, "	Wykryte uderzenia:", len(local_angles_list))
		angles_list.extend(local_angles_list)
		lock.release()

#mozna tez wrzucic map dla lepszego porzadku ale to tylko zwalnia nie wplywa na
#wyniki
def main():
	with ThreadPoolExecutor() as executor:
		executor.map(asteroid_batch, list(range(multiprocessing.cpu_count())))

main()