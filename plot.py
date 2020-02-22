# Penrose tiling - performs deflation to generate the pattern

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.colors import to_rgba

class Tile:
	def __init__(self, xy, label, color='blue', alpha=0.5):
		self.patch = Polygon(xy,
			facecolor=tuple(np.array([i for i in to_rgba(color)])-[0,0,0,1-alpha]),
			edgecolor='k',
			linewidth=0.3)

		self.label = label
		self.xy = xy

	def get_patch(self):
		return self.patch

	def get_label(self):
		return self.label

	def get_xy(self):
		return self.xy

	def set_patch(self, xy, color='blue', alpha=0.5):
		self.patch = Polygon(xy,
			facecolor=tuple(np.array([i for i in to_rgba(color)])-[0,0,0,1-alpha]),
			edgecolor='k',
			linewidth=0.3)

	def subdivide(self):
		"""
			returns a list of tiles according to subdivision rules
		"""

		label = self.label
		xy = self.xy

		# Base length of the given tile (triangle) - assumes first two vertices determine the base
		base = np.sqrt((xy[0][0]-xy[1][0])**2+(xy[0][1]-xy[1][1])**2)

		# Determine the angle of rotation of the given tile about its first vertex
		# The ordering of the vertices are important: The base vertices should be ordered from L->R
		orient = np.sign(xy[1][0]-xy[0][0])
		rotate = np.arctan((xy[1][1]-xy[0][1])/(xy[1][0]-xy[0][0]))+(orient-1)/2*np.pi

		# Determine the translation to apply to the subdivisions
		shift_x, shift_y = xy[0][0], xy[0][1]

		# The given tile will be labeled either thinLeft, thickLeft, thinRight, or thickRight, from which we determine how to subdivide
		if label == 'thinLeft':
			scale = base/np.sqrt((base/thin_width*thin_height)**2+(base/2)**2)
			h = 2*np.sqrt(base**2-(scale*base/thin_width*thin_height)**2)
			alpha = np.arctan(thin_height/(0.5*thin_width))

			p1 = Tile(np.array([[shift_x+h*np.cos(alpha+rotate), shift_y+h*np.sin(alpha+rotate)],
				[shift_x, shift_y],
				[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)]]),
				'thinLeft', color='blue')

			dist = np.sqrt((base/2)**2+(thin_height*base/thin_width)**2)
			p2 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
				[shift_x+dist*np.cos(alpha+rotate), shift_y+dist*np.sin(alpha+rotate)],
				[shift_x+h*np.cos(alpha+rotate), shift_y+h*np.sin(alpha+rotate)]]),
				'thickLeft', color='red')

			return [p1, p2]

		if label == 'thinRight':
			scale = base/np.sqrt((base/thin_width*thin_height)**2+(base/2)**2)
			h = 2*np.sqrt(base**2-(scale*base/thin_width*thin_height)**2)
			alpha = np.arctan(thin_height/(0.5*thin_width))
			
			radius = np.sqrt((h*np.sin(alpha))**2+(base-h*np.cos(alpha))**2)
			theta = np.arctan(h*np.sin(alpha)/(base-h*np.cos(alpha)))+rotate

			p1 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
				[shift_x+radius*np.cos(theta), shift_y+radius*np.sin(theta)],
				[shift_x, shift_y]]),
				'thinRight', color='blue')

			dist = np.sqrt((base/2)**2+(thin_height*base/thin_width)**2)
			p2 = Tile(np.array([[shift_x+dist*np.cos(alpha+rotate), shift_y+dist*np.sin(alpha+rotate)],
				[shift_x, shift_y],
				[shift_x+radius*np.cos(theta), shift_y+radius*np.sin(theta)]]),
				'thickRight', color='red')

			return [p1, p2]

		if label == 'thickLeft':
			bh = base/thick_width*thick_height
			x = np.sqrt((base/2)**2+bh**2)
			h = x/thick_width*thick_height
			l2 = np.sqrt((x/2)**2+h**2)
			l1 = base - l2
			l1h = thick_height*l1/thick_width

			p1 = Tile(np.array([[shift_x, shift_y],
				[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)],
				[shift_x+l1/2*np.cos(rotate)-l1h*np.sin(rotate),
					shift_y+l1h*np.cos(rotate)+l1/2*np.sin(rotate)]]),
				'thickRight', color='red')

			p2 = Tile(np.array([[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
				[shift_x+base/2*np.cos(rotate)-bh*np.sin(rotate),
					shift_y+base/2*np.sin(rotate)+bh*np.cos(rotate)],
				[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
				'thickLeft', color='red')

			p3 = Tile(np.array([[shift_x+base/2*np.cos(rotate)-bh*np.sin(rotate),
					shift_y+base/2*np.sin(rotate)+bh*np.cos(rotate)],
				[shift_x+l1/2*np.cos(rotate)-l1h*np.sin(rotate),
					shift_y+l1/2*np.sin(rotate)+l1h*np.cos(rotate)],
				[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
				'thinRight', color='blue')
		
			return [p1, p2, p3]

		if label == 'thickRight':
			bh = base/thick_width*thick_height
			x = np.sqrt((base/2)**2+bh**2)
			h = x/thick_width*thick_height
			l1 = np.sqrt((x/2)**2+h**2)
			l2 = base - l1
			l2h = thick_height*l2/thick_width

			p1 = Tile(np.array([[shift_x+base/2*np.cos(rotate)-bh*np.sin(rotate),
					shift_y+bh*np.cos(rotate)+base/2*np.sin(rotate)],
				[shift_x, shift_y],	
				[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
				'thickRight', color='red')

			p2 = Tile(np.array([[shift_x+(l1+l2/2)*np.cos(rotate)-l2h*np.sin(rotate),
					shift_y+(l1+l2/2)*np.sin(rotate)+l2h*np.cos(rotate)],
				[shift_x+base/2*np.cos(rotate)-bh*np.sin(rotate),
					shift_y+base/2*np.sin(rotate)+bh*np.cos(rotate)],
				[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)]]),
				'thinLeft', color='blue')

			p3 = Tile(np.array([[shift_x+l1*np.cos(rotate), shift_y+l1*np.sin(rotate)],
				[shift_x+base*np.cos(rotate), shift_y+base*np.sin(rotate)],
				[shift_x+(l1+l2/2)*np.cos(rotate)-l2h*np.sin(rotate),
					shift_y+(l1+l2/2)*np.sin(rotate)+l2h*np.cos(rotate)]]),
				'thickLeft', color='red')
			
			return [p1, p2, p3]

		# Return empty list if label is not one of the four prototile labels
		return []

thin_height = np.sin(np.radians(72))
thin_width = (np.sqrt(5)-1)/2

thick_height = np.sin(np.radians(36))
thick_width = (1+np.sqrt(5))/2

# Prototiles: thinLeft, thickLeft, thinRight, thickRight
thinLeft = Tile(np.array([[0,0], [thin_width,0], [thin_width/2, thin_height]]), 'thinLeft', color='blue')
thickLeft = Tile(np.array([[0,0], [thick_width,0], [thick_width/2, thick_height]]), 'thickLeft', color='red')
thinRight = Tile(np.array([[0,0], [thin_width,0], [thin_width/2, thin_height]]), 'thinRight', color='blue')
thickRight = Tile(np.array([[0,0], [thick_width,0], [thick_width/2, thick_height]]), 'thickRight', color='red')

fig, ax = plt.subplots(figsize=(800/100., 800/100.), dpi=100)

# Set up the axes
ax.set_aspect('equal')
ax.set_xlim(0,thick_width)
ax.set_ylim(-thick_height,thick_height)
plt.axis('off')

# Initial triangles/tiles to start with (same dimensions as prototiles)
angle = np.radians(0) # you can change this angle
base = thick_width    # you can specify a different base length

dist = np.sqrt((base/2)**2+(thick_height)**2)
test_thickLeft = Tile(np.array([[0,0],
	[base*np.cos(angle), base*np.sin(angle)],
	[dist*np.cos(angle+np.arctan(thick_height/(0.5*base))),
	dist*np.sin(angle+np.arctan(thick_height/(0.5*base)))]]),
	'thickLeft', color='green')

angle = np.radians(180) # you can change this angle
test_thickRight = Tile(np.array([[base*-np.cos(angle),base*-np.sin(angle)],
	[base*-np.cos(angle)+base*np.cos(angle), base*np.sin(angle)+base*-np.sin(angle)],
	[base*-np.cos(angle)+dist*np.cos(angle+np.arctan(thick_height/(0.5*base))),
	dist*np.sin(angle+np.arctan(thick_height/(0.5*base))+base*-np.sin(angle))]]),
	'thickRight', color='green')

# The number of times we wish to subdivide the tiles
num_iterations = 6

# Perform the subdivision
initial_tiles = test_thickLeft.subdivide()
tiles = []

for i in range(num_iterations):
	for tile in initial_tiles:
		t = tile.subdivide()
		for a in t:
			tiles.append(a)
	initial_tiles = tiles
	tiles = []

for tile in initial_tiles:
	ax.add_patch(tile.get_patch())

initial_tiles = test_thickRight.subdivide()

for i in range(num_iterations):
	for tile in initial_tiles:
		t = tile.subdivide()
		for a in t:
			tiles.append(a)
	initial_tiles = tiles
	tiles = []

for tile in initial_tiles:
	ax.add_patch(tile.get_patch())

#fig.savefig('subdiv'+str(num_iterations)+'.png', quality=100, dpi=100)
plt.show()