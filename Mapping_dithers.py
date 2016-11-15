from numpy import *
from string import *
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from shapely.ops import cascaded_union

#map_file = 'IFUcen_HETDEX.txt'		#VIRUS unit
map_file = '/Users/Briana/Documents/cure/virusp1/config/IFUcen_VP2_27m.txt'	#VIRUS-P
#map_file = 'LRS2_B_UV_mapping.txt'	#UltraViolet
#map_file = 'LRS2_B_OR_mapping.txt'	#Orange
#map_file = 'LRS2_R_NR_mapping.txt'	#(Near) Red
#map_file = 'LRS2_R_FR_mapping.txt'	#Far Red

dith_file = 'dith2.txt'

fib_lis = open(map_file).read().splitlines()
fib_inf = split(fib_lis[1],' ')

fibsize = float(fib_inf[0])
fibsep = float(fib_inf[-1])

print 'Fiber Size: '+str(fibsize)
print 'Fiber Separation: '+str(fibsep)

fib_num, fib_x, fib_y = loadtxt(map_file, usecols = (0,1,2), skiprows=3, unpack=True)
dith_x, dith_y = loadtxt(dith_file, usecols = (0,1), unpack=True)

plt.axes()

dith_color = ['blue', 'green', 'purple', 'orange', 'red', 'black']

for d in range(len(dith_x)):
	for i in range(len(fib_num)):
		x = fib_x[i]+dith_x[d]
		y = fib_y[i]+dith_y[d]
		s = str(int(fib_num[i]))
		circle = plt.Circle((x, y), radius=(fibsize/2), ec=dith_color[d], fc = 'none')
		plt.gca().add_patch(circle)
		#plt.text(x+0.2, y-0.1, s, fontsize=8)

R0=float64(15) # radius in arcsec arround fiber0 to calculate barycenter
point = plt.Circle((-5.22883324015,-3.97419367307),radius=1.5, ec = 'red',fc='red')
plt.gca().add_patch(point)
point = plt.Circle((-5.22883324015,-3.97419367307),radius=R0, ec = 'red',fc='none')
plt.gca().add_patch(point)

z_line = zeros(10)
x_line = linspace(max(fib_x)+fibsep,min(fib_x)-fibsep,10)
y_line = linspace(max(fib_y)+fibsep,min(fib_y)-fibsep,10)

plt.plot(x_line, z_line, color='red')
plt.plot(z_line, y_line, color='red')

plt.title('6 Dither Fiber Mapping for '+str(map_file))

plt.gca().invert_xaxis()
plt.axis('scaled')
plt.show()
