#!/usr/bin/env python3
'''Make a histogram for a metric_csv_file
Usage:
    ./make_histogram.py metric_csv_file
'''

import sys

import numpy as np
import matplotlib.pyplot as plt


def load_csv_data(csv_file):
    '''Load CSV data into a list of lists.
    Each inner list is a column of data. 
    '''
    with open(csv_file, 'r') as f:
        lines = f.readlines()

        data = []
        for d in lines[0].split(','):
            data.append([])

        for line in lines:
            for j, d in enumerate(line.split(',')):

                data[j].append(d)

    return data

def make_histogram(data, column_id=1):
    '''Make histogram of a column of a data matrix.'''
    data_col = [float(x) for x in data[column_id]]
   
    print('min = {0}'.format(min(data_col)))
    print('max = {0}'.format(max(data_col)))
    print('mean = {0}'.format(np.mean(data_col)))
    print('standard deviation = {0}'.format(np.std(data_col)))
    print('median = {0}'.format(np.median(data_col)))
    print('3rd quartile - 1st quartile = {0}'.format(np.percentile(data_col, 75) - np.percentile(data_col, 25))) 

    median = np.median(data_col)
    percentile_low = np.percentile(data_col, 10)
    percentile_up = np.percentile(data_col, 90)
    upper_cut = median + 1.5 * (percentile_up - percentile_low)
    lower_cut = median - 1.5 * (percentile_up - percentile_low)

    num_bins = 100
    bin_width = (upper_cut - lower_cut) / num_bins

    bins = [lower_cut + bin_width * i for i in range(num_bins)]

    hist, bin_edges = np.histogram(data_col, bins=bins)

    plt.bar(bin_edges[0:-1], hist, width=bin_width)
    plt.show()

if __name__ == '__main__':
    metric_csv_file = sys.argv[1]

   
    make_histogram(load_csv_data(metric_csv_file))
