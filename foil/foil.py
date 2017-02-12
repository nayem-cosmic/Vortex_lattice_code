import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from pylab import savefig

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
f = open("foil_mesh.scr", 'w')


c = .2
s = 1.
m = .04
p = .4
t = .12
ib = 50
jb = 10

dx = c/ib
dy = s/jb

xco = []
yco = []
zco = []

for j in range(jb+1):
    y = dy*j
    for i in range(ib+1):
            x = dx*i
            if x<p*c:
                grad = 2*m*(p-x/c)/(p**2)
                zc=m*x*(2*p-x/c)/(p**2)

            else:
                grad = 2*m*(p-x/c)/((1-p)**2)
                zc=m*(c-x)*(1+x/c-2*p)/((1-p)**2)


            zt = 5*t*(.2969*np.sqrt(x/c)-.1260*(x/c)-.3516*np.power(x/c,2)+.2843*np.power(x/c,3)-.1015*np.power(x/c,4))
            xu = x-zt*np.sin(np.arctan(grad))
            zu = zc+zt*np.cos(np.arctan(grad))
            xco.append(xu)
            yco.append(y)
            zco.append(zu)
    for i in range(ib,-1,-1):
             x = dx*i
             if x<p*c:
                 grad = 2*m*(p-x/c)/p**2
                 zc=m*x*(2*p-x/c)/p**2
 
             else:
                 grad = 2*m*(p-x/c)/np.power((1-p), 2)
                 zc=m*(c-x)*(1+x/c-2*p)/(1-p)**2
             
             zt = 5*t*(.2969*np.sqrt(x/c)-.1260*(x/c)-.3516*np.power(x/c,2)+.2843*np.power(x/c,3)-.1015*np.power(x/c,4))
             xl = x+zt*np.sin(np.arctan(grad))
             zl = zc-zt*np.cos(np.arctan(grad))
             xco.append(xl)
             yco.append(y)
             zco.append(zl) 
for i in range(ib+1):
    x = dx*i
    if x<p*c:
        grad = 2*m*(p-x/c)/p**2
        zc=m*x*(2*p-x/c)/p**2
    else:
        grad = 2*m*(p-x/c)/(1-p)**2
        zc=m*(c-x)*(1+x/c-2*p)/(1-p)**2
 
    zt = 5*t*(.2969*np.sqrt(x/c)-.1260*(x/c)-.3516*np.power(x/c,2)+.2843*np.power(x/c,3)-.1015*np.power(x/c,4))
    xu = x-zt*np.sin(np.arctan(grad))
    zu = zc+zt*np.cos(np.arctan(grad))
    if i%2==0:
        for j in range(jb+1):
            y = dy*j
            xco.append(xu)
            yco.append(y)
            zco.append(zu)
    else:
        for j in range(jb,-1,-1):
            y = dy*j
            xco.append(xu)
            yco.append(y)
            zco.append(zu)
for i in range(ib,-1,-1):
    x = dx*i
    if x<p*c:
        grad = 2*m*(p-x/c)/p**2
        zc=m*x*(2*p-x/c)/p**2
    else:
        grad = 2*m*(p-x/c)/(1-p)**2
        zc=m*(c-x)*(1+x/c-2*p)/(1-p)**2
 
    zt = 5*t*(.2969*np.sqrt(x/c)-.1260*(x/c)-.3516*np.power(x/c,2)+.2843*np.power(x/c,3)-.1015*np.power(x/c,4))
    xl = x+zt*np.sin(np.arctan(grad))

    zl = zc-zt*np.cos(np.arctan(grad))

    if i%2==0:
        for j in range(jb+1):
            y = dy*j
            xco.append(xl)
            yco.append(y)
            zco.append(zl)
    else:
        for j in range(jb,-1,-1):
            y = dy*j
            xco.append(xl)
            yco.append(y)
            zco.append(zl)
 
ax.plot_wireframe(xco, yco, zco,color='purple')
# rstride : array row stride (step size), defaults to 1
# cstride : array column stride (step size), defaults to 1
#ax.auto_scale_xyz([0,4],[0,13],[-.5,0])
#ax.set_aspect('equal'): this only works for 2D
plt.xlabel('Chordwise Length')
plt.ylabel('Spanwise Length')

savefig('foil_mesh.png')

f.truncate()
f.write("PFACE\n")
i=0
for data in xco:
    f.write(str(xco[i]))
    f.write(",")
    f.write(str(yco[i]))
    f.write(",")
    f.write(str(zco[i]))
    f.write("\n")
    i=i+1
f.close()

plt.show()
