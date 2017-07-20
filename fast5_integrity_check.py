#!/usr/bin/env python3

import os
import argparse
import h5py
import sys


def main():
    args = get_arguments()

    print('\nLooking for fast5 files in: ' + args.dir, file=sys.stderr)
    fast5_files = find_all_fast5s(args.dir)
    print('  Found ' + int_to_str(len(fast5_files)) + ' reads\n', file=sys.stderr)
    if not fast5_files:
        sys.exit()

    print('Checking file integrity', file=sys.stderr, end='')
    good_fast5_count = 0
    bad_fast5_count = 0
    last_read_good = True
    for fast5_file in fast5_files:
        try:
            hdf5_file = h5py.File(fast5_file, 'r')
            get_hdf5_names(hdf5_file)
            good_fast5_count += 1
            last_read_good = True
            if good_fast5_count % 100 == 0:
                print('.', file=sys.stderr, end='', flush=True)
        except (IOError, RuntimeError):
            bad_fast5_count += 1
            if last_read_good:
                print('', file=sys.stderr, end='', flush=True)
            print(fast5_file)
            last_read_good = False

    print('\nResults:', file=sys.stderr)
    print('  ' + int_to_str(good_fast5_count) +
          ' good fast5 file' + ('' if good_fast5_count == 1 else 's'), file=sys.stderr)
    print('  ' + int_to_str(bad_fast5_count) +
          ' bad fast5 file' + ('' if good_fast5_count == 1 else 's'), file=sys.stderr)


def get_arguments():
    parser = argparse.ArgumentParser(description='FAST5 integrity check '
                                                 '(prints bad fast5 files to stdout)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dir', type=str,
                        help='directory of FAST5 reads to check (will be searched recursively)')
    args = parser.parse_args()
    args.dir = os.path.abspath(args.dir)
    return args


def find_all_fast5s(directory):
    fast5s = []
    for dir_name, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fast5'):
                fast5s.append(os.path.join(dir_name, filename))
    return fast5s


def get_hdf5_names(hdf5_file):
    names = []
    hdf5_file.visit(names.append)
    return names


def int_to_str(num):
    return '{:,}'.format(num)


if __name__ == '__main__':
    main()
