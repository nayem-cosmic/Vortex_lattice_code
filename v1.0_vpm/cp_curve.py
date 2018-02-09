import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
#from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/cp_curve.txt', 'r')

x_list = []
y_list = []

    
heading = f.readline()

for l in f:
    lstr = l.split()
    x = float(lstr[0])
    x_list.append(x)
    y = float(lstr[1])
    y_list.append(y)

#finter = interp1d(x_list, y_list, kind='cubic')
#xnew_list = np.linspace(min(x_list), max(x_list), num = 50)
#plt.plot(xnew_list, finter(xnew_list), color='black', linewidth=.5, linestyle='-', marker='o', markerfacecolor='r', label = 'Line Graph')
plt.plot(x_list, y_list, color='black', linewidth=1, linestyle='-', marker='o', markerfacecolor='red')

plt.title("Pressure Distribution on Upper and Lower Surface")
plt.xlabel("Chordwise Length, x")
plt.ylabel("Pressure Coefficient, Cp")
plt.grid()

#savefig('figures/alphacl.png')

plt.show()

f.close()
