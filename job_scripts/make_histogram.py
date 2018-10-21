#!/usr/bin/env python3
'''Make a histogram for a metric_csv_file
Usage:
    ./make_histogram.py [output_file] metric_csv_file1 [metric_csv_file2 ...]
'''

from optparse import OptionParser

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

def make_histogram(data_list, output_file=None, column_id=1, normalize=False, data_labels=None, xlabel=None):
    '''Make histogram of a column of a data matrix.'''
    
    for k, data in enumerate(data_list):
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

        if normalize:
            hist = [h / max(hist) for h in hist]

        if data_labels:
            plt.bar(bin_edges[0:-1], hist, width=bin_width, alpha=1.0/len(data_list), label=data_labels[k])
        else:    
            plt.bar(bin_edges[0:-1], hist, width=bin_width, alpha=1.0/len(data_list))
   
    plt.legend()

    if xlabel:
        plt.xlabel(xlabel)

    if normalize:
        plt.ylabel('Normalized count')
    else:
        plt.ylabel('Count')

    if output_file: 
        plt.title(output_file.split('.')[0])
        plt.savefig(output_file)
    else:
        plt.show()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="output")
    parser.add_option("-n", action="store_true", dest="normalize")
    parser.add_option("-c", "--column", default=1, dest="column")
    parser.add_option("-x", "--xlabel", default=None, dest="xlabel")
    (options, args) = parser.parse_args()

    metric_csv_files = args
    output_file = options.output
  
    data_list = [load_csv_data(f) for f in metric_csv_files]

    make_histogram(data_list, output_file, normalize=options.normalize, data_labels=metric_csv_files, column_id=int(options.column), xlabel=options.xlabel)
