
import scipy
import scipy
import math
import math


number_of_asteroids = 1000

class vector: # wektor euklidesowy dla oddziaływań Słońce-Ziemia-Księżyc
class vector: # wektor euklidesowy dla oddziaływań Słońce-Ziemia-Księżyc
	x = float()
	x = float()
	y = float()
	y = float()


	def __init__(x_arg: float, y_arg: float):
	def __init__(x_arg: float, y_arg: float):
		x, y = x_arg, y_arg
		x, y = x_arg, y_arg


	def __init__():
		x, y = 0, 0

	def lenght():
	def lenght():
		return math.hypot(x,y)
		return math.hypot(x,y)


	def __add__(self, vec2: vector):
	def __add__(self, vec2):
		return vector(self.x + vec2.x, self.y + vec2.y)
		return vector(self.x + vec2.x, self.y + vec2.y)


	def __sub__(self, vec2: vector):
	def __sub__(self, vec2):
		return vector(self.x - vec2.x, self.y - vec2.y)
		return vector(self.x - vec2.x, self.y - vec2.y)


	def __iadd__(self, vec2: vector):
	def __iadd__(self, vec2):
		temp = self + vec2
		temp = self + vec2
		x = temp.x
		x = temp.x
		y = temp.y
		y = temp.y
		return temp
		return temp


	def __isub__(self, vec2: vector):
	def __isub__(self, vec2):
		temp = self - vec2
		temp = self - vec2
		x = temp.x
		x = temp.x
		y = temp.y
		y = temp.y
	def __div__(self, divisor: float):




class celestial_body: # ciało Słońce/Ziemia/Księżyc
class celestial_body: # ciało Słońce/Ziemia/Księżyc
	position = vector()
	position : vector
	velocity = vector()
	velocity : vector
	acceleration = vector()
	acceleration : vector
	sgp = 0. # standardowy parametr grawitacyjny
	sgp = 0. # standardowy parametr grawitacyjny
	radius = 0.
	radius = 0.


	def interact(self, cb2: celestial_body):
	def interact(self, cb2):
		vec_to_cb2 = cb2.position - self.position
		vec_to_cb2 = cb2.position - self.position
		distance = (self.position - cb2.position).lenght()
		distance = (self.position - cb2.position).lenght()
		acceleration += vec_to_cb2 * sgp / distance ** 3
		acceleration += vec_to_cb2 * sgp / distance ** 3


	def update():
	def update():
		velocity += acceleration
		velocity += acceleration
		position += velocity
		position += velocity
		acceleration = vector(0.,0.)
		acceleration = vector(0.,0.)



class simualtion:
	asteroid_x = np.zeros(shape=(number_of_asteroids)) # numpy array'e asteroidów wrzuci się do oddzielnej klasy
	asteroid_y = np.zeros(shape=(number_of_asteroids))
	asteroid_velocity_x = np.zeros(shape=(number_of_asteroids))
	asteroid_velocity_y = np.zeros(shape=(number_of_asteroids))
	asteroid_acceleration_x = np.zeros(shape=(number_of_asteroids))
	asteroid_acceleration_y = np.zeros(shape=(number_of_asteroids))
	Earth : celestial_body
	Sun : celestial_body
	Moon : celestial_body
	collision_angles = []
	asteroid_collision_mask = np.ones(shape=(number_of_asteroids)) # maska do oznaczania, czy asteroida uderzyła w ciało

	def __init__():
		# tu będzie implementacja rozmieszczania ciał niebieskich i asteroid
		return

	def simulate(): # główna metoda symulacji
		while True: #warunek przzerwania symulacji zostanie dodany potem
			#interakcje grawitacyjne ciał
			Earth.interact(Sun)
			Earth.interact(Moon)
			Sun.interact(Earth) # można by przykleić Słońce w punkcie (0,0), do rozważenia
			Sun.interact(Moon)
			Moon.interact(Earth)
			Moon.interact(Sun)
			for body in [Sun, Earth, Moon]:
				asteroid_vec_to_body_x = numpy.full((number_of_asteroids), body.position.x) - asteroid_x
				asteroid_vec_to_body_y = numpy.full((number_of_asteroids), body.position.y) - asteroid_y
				asteroid_distance = np.hypot(asteroid_vec_to_body_x, asteroid_vec_to_body_y)
				asteroid_acceleration_x += asteroid_vec_to_body_x * body.sgp / asteroid_distance ** 3
				asteroid_acceleration_y += asteroid_vec_to_body_y * body.sgp / asteroid_distance ** 3
			#aktualizacje pozycji
			Earth.update()
			Moon.update()
			Sun.update()
			asteroid_velocity_x += asteroid_acceleration_x
			asteroid_velocity_y += asteroid_acceleration_y
			asteroid_x += asteroid_velocity_x
			asteroid_y += asteroid_velocity_y
			asteroid_acceleration_x = np.zeros(shape=(number_of_asteroids)) #czyszczenie przypieszeń po kroku
			asteroid_acceleration_y = np.zeros(shape=(number_of_asteroids))
			#kolizje ogarnie się potem 
