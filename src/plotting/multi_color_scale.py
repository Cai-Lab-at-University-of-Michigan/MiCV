import numpy as np
from colour import Color

class MultiColorScale():

	m = 1
	gamma = 2.4

	def __init__(self, n_levels=None):
		if (n_levels is None):
			n_levels = 32
		self.scale_dict = self.generate_scale_dict(n_levels)
		self.value_dict = self.generate_value_dict()
		self.scale_list = self.generate_scale_list()
		self.value_list = self.generate_value_list()

	# provides color hex key for each value
	def generate_scale_dict(self, n_levels):
	    x = 0
	    interp = np.linspace(0, self.m, n_levels)
	    #print("[DEBUG] interp: " + str(interp))
	    s = {}
	    for i in interp:
	        for j in interp:
	            for k in interp:
	                x = str((i+j+k) / (3.0*self.m))
	                s[x] = Color(rgb=(1 - i**self.gamma, 1 - j**self.gamma, 1 - k**self.gamma)).hex_l
	    return s

	# provides value for each color hex key
	def generate_value_dict(self):
		return {v: k for k, v in self.scale_dict.items()}

	def generate_scale_list(self):
		return list(map(list, self.scale_dict.items()))

	def generate_value_list(self):
		return list(map(list, self.value_dict.items()))

	def get_scale_list(self):
		return self.scale_list

	def get_scale_dict(self):
		return self.scale_dict

	def get_value_list(self):
		return self.value_list

	def get_value_dict(self):
		return self.value_dict

	def calculate_hex_color(self, RGB_color):
		if (len(RGB_color) == 3):
			r,g,b = RGB_color
			hex_color = Color(rgb=(1 - r**self.gamma, 1 - g**self.gamma, 1 - b**self.gamma)).hex_l
			return hex_color
		else:
			print("[ERROR] RGB color is not len==3")
			return "#000000"

	def calculate_value(self, RGB_color):
		hex_color = self.calculate_hex_color(RGB_color)
		return self.value_dict[hex_color]

if __name__ == "__main__":
	c = MultiColorScale(81)
	print((c.get_scale_list())[0:10])