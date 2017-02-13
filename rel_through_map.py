from numpy import *
from string import *
import pyfits
import subprocess
import matplotlib.pyplot as plt
import matplotlib.lines as lines

#$CUREBIN/fiberextract -d mastertrace_004_L.dist -f mastertrace_004_L.fmod -c -n 2048 mastertrace_004_L.fits

map_file = '/Users/Briana/Documents/cure_v2/cure/config/IFUcen_HETDEX.txt'

specID = ['020']
pix_min = 200
pix_max = 1800

imL = pyfits.open('Femastertrace_'+specID[0]+'_L.fits')
dL =  imL[0].data
hL = imL[0].header

imR = pyfits.open('Femastertrace_'+specID[0]+'_R.fits')
dR =  imR[0].data
hR = imR[0].header

Left = dL[:,pix_min:pix_max]
Right = dR[:,pix_min:pix_max]

totfluxL = divide(sum(Left, axis = 1),pix_max-pix_min)
totfluxR = divide(sum(Right, axis = 1),pix_max-pix_min)

print len(totfluxL),len(totfluxR)
print average(totfluxL),std(totfluxL)
print average(totfluxR),std(totfluxR)

#fiber mapping
fib_lis = open(map_file).read().splitlines()
fib_inf = split(fib_lis[14],' ')

fibsize = float(fib_inf[0])
fibsep = float(fib_inf[-1])

print 'Fiber Size: '+str(fibsize)
print 'Fiber Separation: '+str(fibsep)

fib_num, fib_x, fib_y = loadtxt(map_file, usecols = (0,1,2), skiprows=25, unpack=True)

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


