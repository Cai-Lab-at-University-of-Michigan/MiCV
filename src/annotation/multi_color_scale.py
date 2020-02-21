import numpy as np

class MultiColorScale():

	def __init__(self, n_levels):
		if (n_levels is None):
			n_levels = 17
		self.scale_dict = self.generate_scale_dict(n_levels)
		self.value_dict = self.generate_value_dict()
		self.scale_list = self.generate_scale_list()
		self.value_list = self.generate_value_list()

	# provides color hex key for each value
	def generate_scale_dict(self, n_levels):
	    x = 0
	    interp = np.linspace(0, 255, n_levels)
	    s = {}
	    for i in interp:
	        for j in interp:
	            for k in interp:
	                s[x] = '#%02x%02x%02x' % (int(i), int(j), int(k))
	                x += 1
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
			hex_color = '#%02x%02x%02x' % (int(RGB_color[0]), int(RGB_color[1]), int(RGB_color[2]))
			return hex_color
		else:
			return "#000000"

	def calculate_value(self, RGB_color):
		hex_color = self.calculate_hex_color(RGB_color)	
		return self.value_dict[hex_color]

if __name__ == "__main__":
	c = MultiColorScale(81)
	print((c.get_scale_list())[0:10])