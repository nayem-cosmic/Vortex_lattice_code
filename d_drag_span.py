import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/d_drag_span.txt', 'r')

x_list = []
y_list = []
x_list_neg = []
xline_list = []
xline_list_neg = []

data = []
for i in range(4):
    data.append(f.readline())
    
dystr = data[3].split()
dy = float(dystr[1])
    
heading = f.readline()

for l in f:
    lstr = l.split()
    x = dy*(float(lstr[0])-1)
    x_list.append(x)
    x = -x-dy
    x_list_neg.append(x)
    xline = dy*(float(lstr[0])-0.5)
    xline_list.append(xline)
    xline = -xline
    xline_list_neg.append(xline)
    y = float(lstr[1])
    y_list.append(y)
    
ax.bar(x_list, y_list, dy, color='#ebe8f9', label='Bar Graph')
ax.bar(x_list_neg, y_list, dy, color='#ebe8f9')

finter = interp1d(xline_list, y_list, kind='cubic')
xnew_list = np.linspace(min(xline_list), max(xline_list), num = 40)
plt.plot(xnew_list, finter(xnew_list), color='black', linewidth=.5, linestyle='-', marker='v', markerfacecolor='b', label = 'Line Graph')
finter_neg = interp1d(xline_list_neg, y_list, kind='cubic')
xnew_list_neg = np.linspace(min(xline_list_neg), max(xline_list_neg), num = 40)
plt.plot(xnew_list_neg, finter_neg(xnew_list_neg), color='black', linewidth=.5, linestyle='-', marker='v', markerfacecolor='b')

plt.title("Spanwise Drag Distribution (Drag per Length)")
plt.text(min(x_list_neg),max(y_list)*0.95,data[0],fontsize=8)
plt.text(min(x_list_neg),max(y_list)*0.91,data[1],fontsize=8)
plt.text(min(x_list_neg),max(y_list)*0.87,data[2],fontsize=8)
plt.xlabel('Spanwise Length')
plt.ylabel('Drag per Span Length')
plt.grid()
plt.legend()

plt.show()

f.close()
