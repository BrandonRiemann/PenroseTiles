# Penrose tiling

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.colors import to_rgba

class Tile:
	def __init__(self, xy, label, color='blue', alpha=0.5):
		self.patch = Polygon(xy, facecolor=tuple(np.array([i for i in to_rgba(color)])-[0,0,0,1-alpha]), edgecolor='k', linewidth=0.3)
		self.label = label
		self.xy = xy

	def get_patch(self):
		return self.patch

	def get_label(self):
		return self.label

	def get_xy(self):
		return self.xy

	def set_patch(self, xy, color='blue'):
		self.patch = Polygon(xy, color=color)

# Prototiles: tL, TL, tR, TR
thin_height = 1.90211303/2
thin_width = 0.61803399

thick_height = 1.17557050/2
thick_width = 1.61803399

tL = Tile(np.array([[0,0], [thin_width,0], [thin_width/2, thin_height]]), 'tL', color='blue')
TL = Tile(np.array([[0,0], [thick_width,0], [thick_width/2, thick_height]]), 'TL', color='red')
tR = Tile(np.array([[0,0], [thin_width,0], [thin_width/2, thin_height]]), 'tR', color='blue')
TR = Tile(np.array([[0,0], [thick_width,0], [thick_width/2, thick_height]]), 'TR', color='red')

def subdivide(tile):
	"""
		tile : Tile object
		returns a new Tile object with subdivision
	"""
	label = tile.get_label()
	xy = tile.get_xy()

	base = np.sqrt((xy[0][0]-xy[1][0])**2+(xy[0][1]-xy[1][1])**2)
	orient = np.sign(xy[1][0]-xy[0][0])
	rotate = np.arctan((xy[1][1]-xy[0][1])/(xy[1][0]-xy[0][0]))+(orient-1)/2*np.pi
	shift_x, shift_y = xy[0][0], xy[0][1]

	# Perform subdivision on prototile, then transform to match vertices of given tile object
	if label == 'tL':
		scale = base/np.sqrt((base/thin_width*thin_height)**2+(base/2)**2)
		h = 2*np.sqrt(base**2-(scale*base/thin_width*thin_height)**2)
		alpha = np.arctan(thin_height/(0.5*thin_width))

		p1 = Tile(np.array([[shift_x+h*np.cos(alpha+rotate), shift_y+h*np.sin(alpha+rotate)],
							[shift_x, shift_y],
							[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)]]),
							'tL', color='blue')

		dist = np.sqrt((base/2)**2+(thin_height*base/thin_width)**2)
		p2 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
							[shift_x+dist*np.cos(alpha+rotate), shift_y+dist*np.sin(alpha+rotate)],
							[shift_x+h*np.cos(alpha+rotate), shift_y+h*np.sin(alpha+rotate)]]),
							'TL', color='red')

		return [p1, p2]

	if label == 'tR':
		scale = base/np.sqrt((base/thin_width*thin_height)**2+(base/2)**2)
		h = 2*np.sqrt(base**2-(scale*base/thin_width*thin_height)**2)
		alpha = np.arctan(thin_height/(0.5*thin_width))
		
		radius = np.sqrt((h*np.sin(alpha))**2+(base-h*np.cos(alpha))**2)
		theta = np.arctan(h*np.sin(alpha)/(base-h*np.cos(alpha)))+rotate

		p1 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
							[shift_x+radius*np.cos(theta), shift_y+radius*np.sin(theta)],
							[shift_x, shift_y]]),
							'tR', color='blue')

		dist = np.sqrt((base/2)**2+(thin_height*base/thin_width)**2)
		p2 = Tile(np.array([[shift_x+dist*np.cos(alpha+rotate), shift_y+dist*np.sin(alpha+rotate)],
							[shift_x, shift_y],
							[shift_x+radius*np.cos(theta), shift_y+radius*np.sin(theta)]]),
							'TR', color='red')

		return [p1, p2]

	if label == 'TL':
		scale = base/thick_width
		height = scale*thick_height
		x = np.sqrt((base/2)**2+height**2)
		scale2 = x/thick_width
		h = thick_height*scale2
		l2 = np.sqrt(1/4*x**2+h**2)
		l1 = base-l2

		p1 = Tile(np.array([[shift_x, shift_y],
							[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)],
							[shift_x+l1/2*np.cos(rotate)-thick_height*l1/thick_width*np.sin(rotate), shift_y+thick_height*l1/thick_width*np.cos(rotate)+l1/2*np.sin(rotate)]]),
							'TR', color='red')

		p2 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
							[shift_x+base/2*np.cos(rotate)-height*np.sin(rotate), shift_y+base/2*np.sin(rotate)+height*np.cos(rotate)],
							[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
							'TL', color='red')

		p3 = Tile(np.array([[shift_x+base/2*np.cos(rotate)-height*np.sin(rotate), shift_y+base/2*np.sin(rotate)+height*np.cos(rotate)],
							[shift_x+l1/2*np.cos(rotate)-thick_height*l1/thick_width*np.sin(rotate), shift_y+l1/2*np.sin(rotate)+thick_height*l1/thick_width*np.cos(rotate)],
							[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
							'tR', color='blue')
	
		return [p1, p2, p3]

	if label == 'TR':
		x = np.sqrt((base/2)**2+(base/thick_width*thick_height)**2)
		h = x/thick_width*thick_height
		l1 = np.sqrt((x/2)**2+h**2)
		l2 = base - l1
		p1 = Tile(np.array([[shift_x+base/2*np.cos(rotate)-base/thick_width*thick_height*np.sin(rotate), shift_y+base/thick_width*thick_height*np.cos(rotate)+base/2*np.sin(rotate)],
							[shift_x, shift_y],	
							[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
							'TR', color='red')

		p2 = Tile(np.array([[shift_x+(l1+l2/2)*np.cos(rotate)-l2/thick_width*thick_height*np.sin(rotate), shift_y+(l1+l2/2)*np.sin(rotate)+l2/thick_width*thick_height*np.cos(rotate)],
							[shift_x+base/2*np.cos(rotate)-base/thick_width*thick_height*np.sin(rotate), shift_y+base/2*np.sin(rotate)+base/thick_width*thick_height*np.cos(rotate)],
							[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
							'tL', color='blue')

		p3 = Tile(np.array([[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)],
							[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
							[shift_x+(l1+l2/2)*np.cos(rotate)-l2/thick_width*thick_height*np.sin(rotate), shift_y+(l1+l2/2)*np.sin(rotate)+l2/thick_width*thick_height*np.cos(rotate)]]),
							'TL', color='red')
		
		return [p1, p2, p3]

	return []

fig, ax = plt.subplots(figsize=(800/100., 800/100.), dpi=100)

ax.set_aspect('equal')
ax.set_xlim(0,thick_width)
ax.set_ylim(-thick_height,thick_height)
plt.axis('off')

angle = np.radians(0)
base = thick_width
dist = np.sqrt((base/2)**2+(thick_height)**2)
test_TL = Tile(np.array([[0,0],
						 [base*np.cos(angle), base*np.sin(angle)],
						 [dist*np.cos(angle+np.arctan(thick_height/(0.5*base))), dist*np.sin(angle+np.arctan(thick_height/(0.5*base)))]]),
						 'TL', color='green')

angle = np.radians(180)
test_TR = Tile(np.array([[thick_width,0],
						 [thick_width+base*np.cos(angle), base*np.sin(angle)],
						 [thick_width+dist*np.cos(angle+np.arctan(thick_height/(0.5*base))), dist*np.sin(angle+np.arctan(thick_height/(0.5*base)))]]),
						 'TR', color='green')

num_iterations = 6
initial_tiles = subdivide(test_TL)
tiles = []

for i in range(num_iterations):
	for tile in initial_tiles:
		t = subdivide(tile)
		for a in t:
			tiles.append(a)
	initial_tiles = tiles
	tiles = []

for tile in initial_tiles:
	ax.add_patch(tile.get_patch())

initial_tiles = subdivide(test_TR)

for i in range(num_iterations):
	for tile in initial_tiles:
		t = subdivide(tile)
		for a in t:
			tiles.append(a)
	initial_tiles = tiles
	tiles = []

for tile in initial_tiles:
	ax.add_patch(tile.get_patch())

#fig.savefig('subdiv'+str(num_iterations)+'.png', quality=100, dpi=100)
plt.show()