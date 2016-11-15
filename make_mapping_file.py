from numpy import * 

#quandrant for the numbering to start in
quad = 4

#chan = 'LRS2_B_UV'  #UltraViolet
#chan = 'LRS2_B_OR'  #Orange
#chan = 'LRS2_R_NR'  #(Near) Red
chan = 'LRS2_R_FR'  #Far Red

side = 'R'

#diameter of the fiber
fibsize = 0.59
#separation of the fibers
fibsep = 0.59 
#number of fiber columns 
numcols = 13
#number of fibers in the longest row. Assumes every other row is 1 shorter.
longrow = 22

f1 = open(str(chan)+'_mapping.txt','w')

#map_file = 'mapping_cols_LRS2.txt'
#fib_num, fib_col = loadtxt(map_file, usecols = (0,1), unpack=True)

#find starting point for fiber 1 
y = float((numcols-1)/2)*fibsep
x = ((longrow*fibsep)/2)-(0.5*fibsep)
print 'X: '+str(x), 'Y: '+str(y)


x_lis = []
y_lis = []
for c in range(numcols):
	if c%2 ==0:
		numrows = longrow
		xOff = 0
	else:
		numrows = longrow-1
		xOff = 0.5*fibsep

	for r in range(numrows):
		x_n = x - (r*fibsep) - xOff
		if -0.01<x_n<0.01:
			x_n = 0
		y_n = y - (c*fibsep)

		x_lis.append(x_n)
		y_lis.append(y_n)

if quad == 2:
	x_lis = multiply(x_lis, -1)
if quad == 3:
	y_lis = multiply(y_lis, -1)
if quad == 4:
	x_lis = multiply(x_lis, -1)
	y_lis = multiply(y_lis, -1)

f1.write('# FIBERD   FIBERSEP'+'\n')
f1.write(str(fibsize)+'  '+str(fibsep)+'\n')
f1.write('# NFIBX	NFIBY'+'\n')
f1.write(str(longrow)+'  '+str(numcols)+'\n')
f1.write('# ID    FIBER_X   FIBER_Y'+'\n')
for i in range(len(x_lis)):
	f1.write(str(i+1)+' '+str(x_lis[i])+' '+str(y_lis[i])+'	'+side+'	'+str(i+1)+'	'+'1.000'+'\n')

f1.close()
print 'File Done'




