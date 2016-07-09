import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import interp1d

#plt.rc('font',family='serif')
plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

f = open('outputdata/lambdacm.txt', 'r')

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

plt.title("Sweep Angle vs. CM")
plt.text(-max(x_list)*0.95,max(y_list)*0.9,data[0],fontsize=8)
plt.text(-max(x_list)*0.95,max(y_list)*0.85,data[1],fontsize=8)
plt.xlabel("Lambda")
plt.ylabel("CM")
plt.grid()

plt.show()

f.close()
