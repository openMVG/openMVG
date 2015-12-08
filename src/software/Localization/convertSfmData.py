#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 11:27:30 2015

@author: sgaspari
"""

import json
import argparse

if __name__ == '__main__':  
    
    parser = argparse.ArgumentParser(description='Convert a sfm_data into the new format with polymorphic View')
    parser.add_argument('-i', '--inputsfm', required=True, help='The sfm data file to convert')
    parser.add_argument('-o', '--outputsfm', required=True, help='The name of the converted sfm data file')
    args = parser.parse_args()

   
    with open(args.inputsfm) as data_file:    
        data = json.load(data_file)
        
    numView = len(data['views']) 
    
    print('Found ' + str(numView) + ' Views')
    
    for i in range(numView):
        data['views'][i]['value']["polymorphic_id"] = 1073741824;
    
    with open(args.outputsfm, 'w') as outfile:
        json.dump(data, outfile, indent=2)