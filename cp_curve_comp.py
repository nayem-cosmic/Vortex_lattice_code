import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig
from scipy.interpolate import interp1d

plt.rc('font',family='serif')
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/cp_curve.txt', 'r')
g = open('naca0012_0.txt','r')

x_list = []
y_list = []
x_list_st = []
y_list_st = []

    
heading_f = f.readline()
heading_g = g.readline()

for l in f:
    lstr = l.split()
    x = float(lstr[0])
    x_list.append(x)
    y = float(lstr[1])
    y_list.append(y)

for l in g:
    lstr = l.split()
    x = float(lstr[0])/100.
    x_list_st.append(x)
    y = float(lstr[1])
    y_list_st.append(y)





#finter = interp1d(x_list, y_list, kind='cubic')
#xnew_list = np.linspace(min(x_list), max(x_list), num = 50)
#plt.plot(xnew_list, finter(xnew_list), color='black', linewidth=.5, linestyle='-', marker='o', markerfacecolor='r', label = 'Line Graph')
plt.plot(x_list, y_list, color='black', linewidth=1, linestyle='-', marker='o', markerfacecolor='red', label='Result')

# standard cp curve
#finter = interp1d(x_list_st, y_list_st, kind='cubic')
#xnew_list = np.linspace(min(x_list_st), max(x_list_st), num = 50)
#plt.plot(xnew_list, finter(xnew_list), color='blue', linewidth=2, linestyle='',marker='o' ,label = 'Standard Data')
plt.plot(x_list_st, y_list_st, color='black', linestyle='', marker='o', markerfacecolor='blue')


plt.title("Pressure Distribution on Upper and Lower Surface")
plt.xlabel("Chordwise Length, x")
plt.ylabel("Pressure Coefficient, Cp")
plt.grid()

#savefig('figures/alphacl.png')

plt.show()

f.close()
