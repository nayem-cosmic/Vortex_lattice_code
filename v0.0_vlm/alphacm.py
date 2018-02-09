import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
#from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/alphacm.txt', 'r')

x_list = []
y_list = []

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
plt.plot(x_list, y_list, color='black', linewidth=1, linestyle='-', marker='s', markerfacecolor='green')

plt.title("Moment Coefficient vs. Angle of Attack")
plt.text(max(x_list)*0.05,max(y_list)*0.95,data[0],fontsize=10)
plt.text(max(x_list)*0.05,max(y_list)*0.91,data[1],fontsize=10)
plt.xlabel("Angle of Attack")
plt.ylabel("Moment Coefficient, CM")
plt.grid()

savefig('figures/alphacm.png')

plt.show()

f.close()
