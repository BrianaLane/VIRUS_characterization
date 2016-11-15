from numpy import * 
from string import *
import matplotlib.pyplot as plt
import matplotlib.lines as lines

#map_file = 'IFUcen_HETDEX.txt'		#VIRUS unit
#map_file = 'IFUcen_VP2_27m.txt'	#VIRUS-P
#map_file = 'LRS2_B_UV_mapping.txt'	#UltraViolet
#map_file = 'LRS2_B_OR_mapping.txt'	#Orange
map_file = 'LRS2_R_NR_mapping.txt'	#(Near) Red
#map_file = 'LRS2_R_FR_mapping.txt'	#Far Red
#map_file = '/Users/Briana/Documents/cure_v2/cure/config/IFUcen_cam020.txt'

fib_lis = open(map_file).read().splitlines()
#fib_inf = split(fib_lis[14],' ')
fib_inf = split(fib_lis[1],' ')
print fib_inf

fibsize = float(fib_inf[0])
fibsep = float(fib_inf[-1])

print 'Fiber Size: '+str(fibsize)
print 'Fiber Separation: '+str(fibsep)

#fib_num, fib_x, fib_y = loadtxt(map_file, usecols = (0,1,2), skiprows=25, unpack=True)
fib_num, fib_x, fib_y = loadtxt(map_file, usecols = (0,1,2), skiprows=5, unpack=True)

plt.axes()

for i in range(len(fib_num)):
	x = fib_x[i]
	y = fib_y[i]
	s = str(int(fib_num[i]))
	circle = plt.Circle((x, y), radius=(fibsize/2), ec='b', fc = 'none')
	plt.gca().add_patch(circle)
	plt.text(x+0.2, y-0.1, s, fontsize=8)

z_line = zeros(10)
x_line = linspace(max(fib_x)+fibsep,min(fib_x)-fibsep,10)
y_line = linspace(max(fib_y)+fibsep,min(fib_y)-fibsep,10)

plt.plot(x_line, z_line, color='red')
plt.plot(z_line, y_line, color='red')


plt.title('Fiber Mapping for '+str(map_file))

plt.gca().invert_xaxis()
plt.axis('scaled')
plt.show()