import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as anim


dir = ""
loop = False
time = 0.1
width = 100
height = 100

j = 1
while j < len(sys.argv):
    if sys.argv[j][:2] == "-d":
        j += 1
        dir = sys.argv[j]
    elif sys.argv[j][:2] == "-l":
        loop = True
    elif sys.argv[j][:2] == "-t":
        j += 1
        time = float(sys.argv[j])
    elif sys.argv[j][:2] == "-w":
        j += 1
        width = float(sys.argv[j])
    elif sys.argv[j][:2] == "-h":
        j += 1
        height = float(sys.argv[j])

    j += 1

if len(dir) == 0:
    print("Missing input directory")

particles_list = []

for file in sorted(os.scandir(dir), key=lambda entry: int((os.path.basename(entry.name)[6:])[:-4])):
    if file.path.endswith(".txt"):
        particles_list.append(pd.read_csv(file.path, header=None, names=['mass', 'x', 'y']))

fig = plt.figure()


def init():
    plt.xlim(-500, 500)
    plt.ylim(-500, 500)
    plt.scatter(particles_list[0].x, particles_list[0].y, s=particles_list[0].mass)


def update(i):
    plt.cla()
    plt.scatter(particles_list[i].x, particles_list[i].y, s=particles_list[i].mass)


ani = anim.FuncAnimation(fig, update, frames=len(particles_list), interval=time*1000)
plt.show()
