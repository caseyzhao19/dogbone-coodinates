import matplotlib.pyplot as plt
import numpy as np

lines = open( 'plotdata.txt', "r" ).readlines()
data = []
for line in lines:
    l = line[1:-1].split()
    l = [float(i[:-1]) for i in l]
    data.append(l)
print(data)

for i in range(138):
    figure, axes = plt.subplots()
    plt.figure(figsize=(4, 7))
    axes.set_ylabel('distance')
    axes.set_ylabel('iteration')
    plt.title("side length = " + "{:.2f}".format(0.02+0.01*i))
    plt.axhline(0, color='lightgray')
    plt.xlim([0.5, 12.5])
    plt.ylim([-0.9, 0.8])
    plt.xticks(range(1, 13))
    plt.yticks(np.arange(-0.9, 0.9, 0.1))
    plt.plot(range(1, 13), data[i])
    plt.savefig('./plots/plot' + str(i) + '.jpg', dpi=600)