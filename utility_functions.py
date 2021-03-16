#! /opt/local/bin/python

import os
import string

import numpy as np

"""
This is just a bunch of functions that I use all the time
"""

def read_file(filename, delimiter=None, startline=0):
    """General function to read text file into a 2D list."""

    data_list = []
    ifile = open(filename,'rU')

    for line in ifile:

        if delimiter:
            data = line.split(delimiter)
        else:
            data = line.split()
        data_list.append(data)

    ifile.close()
    return data_list[startline:]

def get_dict_list(data_array):
    """Returns a list of dictionaries based on the column headers
    (the 0th line in the column headers)

    """

    key_list = data_array[0]
    dict_list = []
    for index, line in enumerate(data_array[1:]):
        params = {}
        for i in range(len(key_list)):
            # try:
            #     params[key_list[i]] = float(line[i])
            # except ValueError:
                params[key_list[i]] = line[i]
        dict_list.append(params)

    return dict_list

def make_dir(path):
    """General function for making a new directory without raising errors"""

    if not os.path.isdir(path):
        os.mkdir(path)


#DEFIITION OF FUNCTION TO SORT A LIST IN PLACE USING A KEY THAT'S CONTAINED WITHIN A STRING
def sort_by_key(list, split_char, key_position):

    def find_key(line):
        
        key = int(line.split(split_char)[key_position])
        return key

    list.sort(key=find_key)

    return list

#DEFINITION OF A FUNCTION TO SAVE AN ARRAY AS A JUSTIFIED TEXT FILE
def save_data_array(array, save_path):
    """Function to write a (square) array with justified column widths"""

    #gets column width
    column_width_list = []
    for column in zip(*array):

        column = map(str,column)
        column_width = max(len(x) for x in column) + 2
        column_width_list.append(column_width)

    #writes array to file
    ofile = open(save_path,'w')

    for i in range(len(array)):
        for j in range(len(array[i])):

            element = str(array[i][j]).ljust(column_width_list[j])
            ofile.write(element + '  ')

        ofile.write('\n')

    ofile.close

#DEFINITION OF FIND FIRST HIGHER INDEX
def find_first_higher_index(list,target):
    
    endindex = 0
    for index, x in enumerate(list):

        if x < target:

            endindex = index

    return endindex

#DEFINITION OF FUNCTION TO FIND INDEX OF VALUE IN A LIST NEAREST TO A TARGET
def find_nearest(list, target):
    
    target_index = (np.abs(list - target)).argmin()

    return target_index

#DEFINITION OF FUNCTION TO FIND INDEX OF VALUE IN A LIST FURTHEST FROM TARGET
def find_furthest(list, target):
    
    target_index = (np.abs(list - target)).argmax()

    return target_index

#DEFINITION OF FUNCTION TO COMPRESS MIXED LIST OF ITERABLE AND NON-SEQUENCE TYPES TO A 1-D LIST
#(THIS FUNCTION WORKS FOR TO INFINITE DIMENSIONS)
def flatten_list(old_list):

    repeat = 'yes'
    while repeat == 'yes':

        new_list = []
        repeat = 'no'
        
        for item in old_list:

            try:
                getattr(item,'__iter__')
                new_list.extend(item)
                repeat = 'yes'
                
            except AttributeError:
                new_list.append(item)

        old_list = new_list
            
    return new_list

#DEFINITION OF SMOOTHING FUNCTION
def smooth_moving_window(l, window_len=11, include_edges='Off'):

    if window_len%2==0:
        raise ValueError('>window_len< kwarg in function >smooth_moving_window< must be odd')

    l = np.reshape(map(float,l),len(l))
    w = np.ones(window_len,'d')
    
    if include_edges == 'On':
        edge_list = np.ones(window_len)
        begin_list = [x * l[0] for x in edge_list]
        end_list = [x * l[-1] for x in edge_list]
        
        s = np.r_[begin_list, l, end_list]
        
        y = np.convolve(w/w.sum(), s , mode='same')
        y = y[window_len + 1:-window_len + 1]
        
    elif include_edges == 'Wrap':
        s=np.r_[2 * l[0] - l[window_len-1::-1], l, 2 * l[-1] - l[-1:-window_len:-1]]
        y = np.convolve(w/w.sum(), s , mode='same')
        y = y[window_len:-window_len+1]

    elif include_edges == 'Off':
        y = np.convolve(w/w.sum(), l, mode='valid')

    else:
        raise NameError('Error in >include_edges< kwarg of function >smooth_moving_window<')
    
    return y
