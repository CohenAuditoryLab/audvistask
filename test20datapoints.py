import matplotlib.pyplot as pyplot
from scipy.interpolate import spline
import numpy as np

def main(file_string):

    data = []  # create empty list to append values to

    input = open(file_string, 'rU')  # open and read csv file
    input.readline()  # skip headers

    for line in input:
        line = line.strip('\n').split(',')  # split lines at commas
        data.append(line)  # add split data to list

    input.close()  # close file

    #split data by coherence & percentage of times for each coherence that subject answered high
    coherence = [0, .25, .5, .75, 1]

    high0 = 0
    high25 = 0
    high50 = 0
    high75 = 0
    high100 = 0

    coh0 = 0
    coh25 = 0
    coh50 = 0
    coh75 = 0
    coh100 = 0

    for entry in data:
        if entry[2] == '0':
            coh0 += 1
            if entry[7] == '2':
                high0 += 1
        if entry[2] == '0.25':
            coh25 += 1
            if entry[7] == '2':
                high25 += 1
        if entry[2] == '0.5':
            coh50 += 1
            if entry[7] == '2':
                high50 += 1
        if entry[2] == '0.75':
            coh75 += 1
            if entry[7] == '2':
                high75 += 1
        if entry[2] == '1':
            coh100 += 1
            if entry[7] == '2':
                high100 += 1

    high0 = float(high0) / coh0 * 100
    high25 = float(high25) / coh25 * 100
    high50 = float(high50) / coh50 * 100
    high75 = float(high75) / coh75 * 100
    high100 = float(high100) / coh100 * 100

    perhigh = [high0, high25, high50, high75, high100]

    # calculate polynomial
    z1 = np.polyfit(coherence, perhigh, 3)
    f1 = np.poly1d(z1)

    # calculate new x's and y's
    x_new = np.linspace(0, 1, 50)
    y_new = f1(x_new)

    pyplot.plot(coherence, perhigh, 'o', x_new, y_new)
    pyplot.ylabel('% Time Subject Chose High')
    pyplot.xlabel('Desired Coherence')
    pyplot.title('Coherence vs. High Selection Rate')
    pyplot.show()

    #create plot reaction time vs. coherence
    coh = []
    rt = []

    for entry in data:
        coh.append(float(entry[2]))
        rt.append(float(entry[12]))

    # calculate polynomial
    z = np.polyfit(coh, rt, 2)
    f = np.poly1d(z)

    # calculate new x's and y's
    x_new = np.linspace(0, 1, 50)
    y_new = f(x_new)

    pyplot.plot(coh, rt, 'o', x_new, y_new)
    pyplot.xlim([-.1, 1.1])
    pyplot.xlabel('Desired Coherence')
    pyplot.ylabel('Reaction Time (ms)')
    pyplot.title('Coherence vs. Reaction Time')
    pyplot.show()

main('AudVisTask_v1_Diana_170524_1154.csv')