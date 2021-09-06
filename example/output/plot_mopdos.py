import numpy as np 
import matplotlib.pyplot as plt

def read_mopdos(filename='MOPDOS.out'):  # read data from MOPDOS.out
    f = open(filename,'r')
    lines = f.readlines()  # read in all lines of file
    f.close()
    data = {} 
    for line in lines:
        if line.find('spin') >=0:
            words = line.split()
            spin = words[1]+words[2]
            data[spin] ={}
            count = 0
        elif line.find('bandnumber') >=0:
            words = line.split()
            band = int(words[-1])
            data[spin][band] ={'energy':[],'mopdos':[]}
            count += 1
        elif len(line) > 1:
            words = line.split()
            data[spin][band]['energy'].append(float(words[0]))
            data[spin][band]['mopdos'].append(float(words[1]))

    return data

def plot_mopdos(ax, data):
    for band in data:
        energy = data[band]['energy']
        mopdos = data[band]['mopdos']
        ax.plot(energy, mopdos, label=band)

    # format axes
    ax.set_xlabel('$E-E_F$ (eV)')
    ax.set_ylabel('MOPDOS (arb. units)')
    plt.legend()

# main program
data = read_mopdos()
fig, ax = plt.subplots()
plot_mopdos(ax, data['spin1'])
plt.savefig('MOPDOS.png')
plt.show()


