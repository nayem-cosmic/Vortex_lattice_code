import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
#from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/chcd.txt', 'r')

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
plt.plot(x_list, y_list, color='black', linewidth=1, linestyle='-', marker='p', markerfacecolor='blue')

plt.title("Drag Coefficient vs. Height Above Ground")
plt.text(max(x_list)*0.5,max(y_list)*0.9,data[0],fontsize=10)
plt.text(max(x_list)*0.5,max(y_list)*0.88,data[1],fontsize=10)
plt.xlabel("Height Above Ground")
plt.ylabel("Drag Coefficient, CD")
plt.grid()

savefig('figures/chcd.png')

plt.show()

f.close()
