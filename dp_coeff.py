import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(121)

f = open('outputdata/dp_coeff.txt', 'r')

x_list = []
y_list = []
z_list = []
y_list_neg = []

data = f.readline()
heading = f.readline()

for l in f:
    lstr = l.split()
    x = float(lstr[0])
    y = float(lstr[1])
    z = float(lstr[2])
    x_list.append(x)
    y_list.append(y)
    z_list.append(z)

for i in y_list:
    i=-i
    y_list_neg.append(i)
    
# pressure coefficient contour graph for right side of wing    
xc = np.linspace(min(x_list), max(x_list), num=500)
yc = np.linspace(min(y_list), max(y_list), num = 500)
x_grid, y_grid = np.meshgrid(xc,yc)
# zip gives paired value
zc = griddata(zip(x_list, y_list), z_list, (x_grid, y_grid), method='cubic')
cp = ax.contourf(xc, yc, zc)

# for left side of wing
ycL = np.linspace(min(y_list_neg), max(y_list_neg), num=500)
x_grid, y_grid_neg = np.meshgrid(xc, ycL)
zcL = griddata(zip(x_list, y_list_neg), z_list, (x_grid, y_grid_neg), method='cubic')
ax.contourf(xc, ycL, zcL)

plt.colorbar(cp)
plt.title("Pressure Coefficient at Different Regions of the Wing")
plt.text(max(x_list)*0.15,max(y_list)*0.85,data,fontsize=8)
plt.xlabel('Chordwise Length')
plt.ylabel('Spanwise Length')

plt.show()

f.close()
