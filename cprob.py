#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Araik Tamazian

import argparse
import numpy as np


def quickSort(alist):
   quickSortHelper(alist,0,len(alist)-1)

def quickSortHelper(alist,first,last):
   if first<last:

       splitpoint = partition(alist,first,last)

       quickSortHelper(alist,first,splitpoint-1)
       quickSortHelper(alist,splitpoint+1,last)


def partition(alist,first,last):
   pivotvalue = alist[first]

   leftmark = first+1
   rightmark = last

   done = False
   while not done:

       while leftmark <= rightmark and alist[leftmark] <= pivotvalue:
           leftmark = leftmark + 1

       while alist[rightmark] >= pivotvalue and rightmark >= leftmark:
           rightmark = rightmark -1

       if rightmark < leftmark:
           done = True
       else:
           temp = alist[leftmark]
           alist[leftmark] = alist[rightmark]
           alist[rightmark] = temp

   temp = alist[first]
   alist[first] = alist[rightmark]
   alist[rightmark] = temp


   return rightmark


def gen_log_space(limit, n):
    """
    # Generates log-spaced integers
    # Taken from https://stackoverflow.com/questions/12418234/logarithmically-spaced-integers
    :param limit: end value of sequence
    :param n: number of elements in sequence
    :return: num logarithmically spaced integers
    """
    ratio = ''
    result = [1]
    if n > 1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result) < n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array(map(lambda x: round(x)-1, result), dtype=np.uint64)


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser('Estimate empirical cumulative '
                                     'probability function (CDF).')
    parser.add_argument('input', help='an input file')
    parser.add_argument('output', help='an output file')
    parser.add_argument('-c', '--compl', action='store_true',
                        default=False, help='estimate complementary CDF')
    parser.add_argument('-n', '--num', default=100,
                        help='number of points')
    parser.add_argument('-l', '--logbin', action='store_true',
                        default=False,
                        help='use logarithmic binning')
    return parser.parse_args()


def estimate_cprob(input_handle, output_handle, compl, num, logbin):
    """
    Given handles of input and output files, read the connection
    statistics from the input file and output the connections grouped
    by their destination to the output handle.

    :type input_handle: file
    :param input_handle: an input handle where source data
        are stored
    :type output_handle: file
    :param output_handle: an output handle where CDF
        is written to
    :param compl: estimate complementary CDF
    :type compl: boolean
    :param logbin: use logarithmic binning
    :type logbin: boolean
    :param num: number of CDF points
    :type num: int
    :return: the tuple of two numbers: the values of random variable
        and the corresponding cumulative probabilities to write
        to the specified output file
    :rtype: tuple
    """
    var = np.loadtxt(input_handle)

    # Get random variable size
    vsize = var.size

    var_sorted = quickSort(var)
    prob = np.linspace(0, 1, vsize)
    if compl:
        prob = 1 - prob
    if logbin:
        ind = gen_log_space(vsize, num)
        var_sorted = var_sorted[ind]
        prob = prob[ind]
    else:
        step = vsize/num
        var_sorted = var_sorted[::step]
        prob = prob[::step]

    np.savetxt(output_handle, np.transpose((var_sorted, prob)), fmt="%f\t%f")

if __name__ == '__main__':
    args = parse_args()
    with open(args.input) as input_file:
        with open(args.output, 'w') as output_file:
            estimate_cprob(input_file, output_file, args.compl, args.num, args.logbin)
