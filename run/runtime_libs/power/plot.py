#!/usr/bin/python

import sys
import getopt
import numpy as np

from matplotlib import pyplot as plt

def main(argv):
    cores = 0
    bench = ''

    fname = 'mic_power.csv'
    
    # Passing one argument: name of the file
    try:
        opts, args = getopt.getopt(argv, "hc:b:")
    except getopt.GetoptError:
        print 'plot.py -f <file_name>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'plot.py -f <file_name>'
            sys.exit()
        elif opt == '-f':
            fname = str(arg)



    time = []
    power = []
    with open(fname, 'r') as f:
        # skip the header
        f.readline()
        
        for line in f:
            columns = line.split(',')
            time.append(float(columns[0]))
            power.append(float(columns[1]))

    fig, ax = plt.subplots()
    line1, = ax.plot(time, power, '-', color = 'blue')

    ax.grid()
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Power (W)')
    ax.set_xlim(min(time), max(time))

    fig.savefig('power.png', bbox_inches='tight')
        

if __name__ == "__main__":
    main(sys.argv[1:])
