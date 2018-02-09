import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
#from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('alphacm.txt', 'r')

x_list = []
y_list = []
x_ref = []
y_ref = []

data = []
for i in range(3):
    data.append(f.readline())
    
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

for i in np.arange(1,24,.01):
    x_ref.append(i)
    y_ref.append(2.743*i*3.1416/180)

plt.plot(x_list, y_list, color='red', linewidth=1, linestyle='-',marker='s',label='Calculated Data')

plt.plot(x_ref, y_ref, color='green', linewidth=1, linestyle='-', label='Standard Data')



plt.title("Moment Coefficient vs. Angle of attack\nComparison with Warren 12 Planform Standard Data")
plt.xlabel("Angle of Attack")
plt.ylabel("Moment Coefficient")
plt.grid()
plt.legend(loc=4)

savefig('alphaml.png')

plt.show()

f.close()
