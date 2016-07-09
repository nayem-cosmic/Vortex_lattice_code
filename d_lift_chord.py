import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/d_lift_chord.txt', 'r')

x_list = []
y_list = []
xline_list = []

data = []
for i in range(4):
    data.append(f.readline())
    
dystr = data[3].split()
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
    
ax.bar(x_list, y_list, dx, color='#fcc2c2', label = 'Bar Graph')

finter = interp1d(xline_list, y_list, kind='cubic')
xnew_list = np.linspace(min(xline_list), max(xline_list), num = 40)
plt.plot(xnew_list, finter(xnew_list), color='black', linewidth=.5, linestyle='-', marker='o', markerfacecolor='r', label = 'Line Graph')

plt.title("Chordwise Lift Distribution (Drag per Length)")
plt.text(max(x_list)*0.7,max(y_list)*0.95,data[0],fontsize=8)
plt.text(max(x_list)*0.7,max(y_list)*0.91,data[1],fontsize=8)
plt.text(max(x_list)*0.7,max(y_list)*0.87,data[2],fontsize=8)
plt.xlabel('Average Chordwise Length')
plt.ylabel('Lift per Chord Length')
plt.grid()
plt.legend(loc=3)

plt.show()

f.close()
