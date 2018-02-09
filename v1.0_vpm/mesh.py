import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import savefig

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

f = open('outputdata/grid.txt', 'r')

x_list = []
y_list = []
z_list = []

#first_line = f.readline()

for l in f:
    lstr = l.split()
    x = float(lstr[0])
    y = float(lstr[1])
    z = float(lstr[2])
    x_list.append(x)
    y_list.append(y)
    z_list.append(z)


ax.plot_wireframe(x_list, y_list, z_list,color='purple')
# rstride : array row stride (step size), defaults to 1
# cstride : array column stride (step size), defaults to 1

#ax.auto_scale_xyz([0,4],[0,13],[-.5,0])
#ax.set_aspect('equal'): this only works for 2D
plt.xlabel('Chordwise Length')
plt.ylabel('Spanwise Length')

#savefig('figures/vlm_mesh.png')

plt.show()

f.close()
