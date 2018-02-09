import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)


fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/d_drag_chord.txt', 'r')

x_list = []
y_list = []
xline_list = []

data = []
for i in range(3):
    data.append(f.readline())
    
dystr = data[2].split()
dx = float(dystr[1])
    
heading = f.readline()

for l in f:
    lstr = l.split()
    x = dx*(float(lstr[0])-1)
    x_list.append(x)
    xline = dx*(float(lstr[0])-0.5)
    xline_list.append(xline)
    y = float(lstr[1])
    y_list.append(y)

#bar graph    
#ax.bar(x_list, y_list, dx, color='#ebe8f9', label='Bar Graph')

finter = interp1d(xline_list, y_list, kind='cubic')
xnew_list = np.linspace(min(xline_list), max(xline_list), num = 50)
plt.plot(xnew_list, finter(xnew_list), color='black', linewidth=.5, linestyle='-', marker='p', markerfacecolor='b', label = 'Line Graph')

plt.title("Chordwise Drag Distribution (Drag per Length)")
plt.text(max(x_list)*0.05,max(y_list)*0.95,data[0],fontsize=10)
plt.text(max(x_list)*0.05,max(y_list)*0.91,data[1],fontsize=10)
plt.xlabel('Average-Chordwise Length')
plt.ylabel('Drag per Chord Length')
plt.grid()
plt.legend(loc=4)

savefig('figures/d_drag_chord.png')

plt.show()

f.close()
