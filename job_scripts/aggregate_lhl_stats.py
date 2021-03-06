#!/usr/bin/env python3

import os
import sys
import json



if __name__ == '__main__':
    
    data_path = sys.argv[1]
    input_path = sys.argv[2]
   
    median_number_of_lhl_units = []
    lhl_lengths = []
    helix_rmsds = []
    helix_rmsds_within_5res_loop = []
    loop_lengths = []

    for f in os.listdir(input_path):

        if f.endswith('.json'):
            with open(os.path.join(input_path, f), 'r') as jf:
                lhl_info = json.load(jf)

                median_number_of_lhl_units.append(lhl_info['median_number_of_lhl_units'])
                lhl_lengths += lhl_info['lhl_lengths']
                helix_rmsds += lhl_info['helix_rmsds']
                helix_rmsds_within_5res_loop += lhl_info['helix_rmsds_within_5res_loop']
                loop_lengths += lhl_info['loop_lengths']

    with open(os.path.join(data_path, 'median_number_of_lhl_units.json'), 'w') as f:
        json.dump(median_number_of_lhl_units, f)

    with open(os.path.join(data_path, 'lhl_lengths.json'), 'w') as f:
        json.dump(lhl_lengths, f)

    with open(os.path.join(data_path, 'helix_rmsds.json'), 'w') as f:
        json.dump(helix_rmsds, f)

    with open(os.path.join(data_path, 'helix_rmsds_within_5res_loop.json'), 'w') as f:
        json.dump(helix_rmsds_within_5res_loop, f)

    with open(os.path.join(data_path, 'loop_lengths.json'), 'w') as f:
        json.dump(loop_lengths, f)
