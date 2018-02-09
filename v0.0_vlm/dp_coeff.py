import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
import numpy as np
from pylab import savefig

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/dp_coeff.txt', 'r')

x_list = []
y_list = []
z_list = []
x_list_neg = []
y_list_rev = []

data = f.readline()
heading = f.readline()

for l in f:
    lstr = l.split()
    x_list.append(float(lstr[1]))
    x_list_neg.append(-float(lstr[1]))
    y_list.append(float(lstr[0]))
    z_list.append(float(lstr[2]))
    
y_max = max(y_list)

for i in range(len(y_list)):
    y_list_rev.append(y_max-y_list[i])
    
# pressure coefficient contour graph for right side of wing    
xc = np.linspace(min(x_list), max(x_list), num=500)
yc = np.linspace(min(y_list_rev), max(y_list_rev), num = 500)
x_grid, y_grid = np.meshgrid(xc,yc)
# zip gives paired value
zc = griddata(zip(x_list, y_list_rev), z_list, (x_grid, y_grid), method='cubic')
cp = ax.contourf(xc, yc, zc)

# for left side of wing
xcL = np.linspace(min(x_list_neg), max(x_list_neg), num=500)
x_grid_neg, y_grid = np.meshgrid(xcL, yc)
zcL = griddata(zip(x_list_neg, y_list_rev), z_list, (x_grid_neg, y_grid), method='cubic')
ax.contourf(xcL, yc, zcL)

plt.colorbar(cp)
plt.title("Pressure Coefficient at Different Regions of the Wing")
plt.text(min(x_list_neg)*0.95,max(y_list)*0.85,data,fontsize=10)
plt.xlabel('Spanwise Length')
plt.ylabel('Chordwise Length')

savefig('figures/dp_coeff.png')

plt.show()

f.close()
