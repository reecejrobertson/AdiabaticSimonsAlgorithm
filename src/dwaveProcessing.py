#---------------------------------------------------------------------------
#                             IMPORTS
#---------------------------------------------------------------------------
import os
from matplotlib import pyplot as plt

#---------------------------------------------------------------------------
#                             CONSTANTS
#---------------------------------------------------------------------------
START = 2
END = 41
SHOTS = 512
PLOT_LAYOUT = (8, 5)
FIGSIZE = (15, 24)
DPI = 300


#---------------------------------------------------------------------------
#                             FUNCIONS
#---------------------------------------------------------------------------
def plotFigure(dictionary, filename):
    fig, axs = plt.subplots(
        PLOT_LAYOUT[0],
        PLOT_LAYOUT[1],
        tight_layout = True,
        sharey = True,
        figsize = FIGSIZE
    )
    for n in range(START, END):
        ax = axs[int((n-1)/PLOT_LAYOUT[1]), int((n-1)%PLOT_LAYOUT[1])]
        penalty = float(dictionary[n]["penalty"])
        time = dictionary[n]["time"]
        bars = ax.bar(
            ['$|0⟩$', '$|n-1⟩$', 'Others'],
            [
                dictionary[n]['zeros'] / SHOTS,
                dictionary[n]['ones'] / SHOTS,
                dictionary[n]['others'] / SHOTS
            ],
            color = ['b', 'g', 'k']
        )
        for bar in bars:
            yval = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                yval + 0.01,
                f'{yval*100:.1f}%',
                ha='center',
                va='bottom'
            )
        ax.set_title(f'$n = {n}, N = {n + 2*(n-1)}, p = {penalty}, t = {time}$')
        ax.set_ylabel('Frequency (%)')
    plt.savefig(f'../{filename}.png', dpi=DPI)
    plt.close()


#---------------------------------------------------------------------------
#                             EXECUTION
#---------------------------------------------------------------------------
optimalZeros = dict()
optimalOnes = dict()
for n in range(START, END):
    optimalZeros[n] = {'zeros': 0, 'ones': 0, 'others': 0, 'penalty': 0, 'time':0}
    optimalOnes[n] = {'zeros': 0, 'ones': 0, 'others': 0, 'penalty': 0, 'time':0}

for filename in os.listdir('../data/'):
    try:
        device, size, penalty, time, protocol = filename.split('-')
    except:
        continue
    n = int(size)
    if device == 'zephyer':
        with open(f'../data/{filename}') as file:
            data = file.readlines()
        zeros = 0
        ones = 0
        others = 0
        for item in data:
            item = item.strip()
            bitstring, value = item.split(':')
            count = int(value)
            inputs, outputs, ancillas = bitstring.split(' ')
            if (
                inputs == '0' * n and
                outputs == '0' * (n - 1) and
                ancillas == '0' * (n-1)
            ):
                zeros += count
            elif (
                inputs == '1' * n and
                outputs == '0' * (n - 1) and
                ancillas == '1' * (n-1)
            ):
                ones += count
            else:
                others += count
        if zeros > optimalZeros[n]['zeros']:
            optimalZeros[n]['zeros'] = zeros
            optimalZeros[n]['ones'] = ones
            optimalZeros[n]['others'] = others
            optimalZeros[n]['penalty'] = penalty
            optimalZeros[n]['time'] = time
        if ones > optimalOnes[n]['ones']:
            optimalOnes[n]['zeros'] = zeros
            optimalOnes[n]['ones'] = ones
            optimalOnes[n]['others'] = others
            optimalOnes[n]['penalty'] = penalty
            optimalOnes[n]['time'] = time

plotFigure(optimalZeros, 'optimalZeros')
plotFigure(optimalOnes, 'optimalOnes')
