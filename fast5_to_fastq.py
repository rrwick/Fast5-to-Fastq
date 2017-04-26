#!/usr/bin/env python3

import os
import argparse
import h5py
import statistics
import sys


def main():
    args = get_arguments()

    print('\nLooking for fast5 files in: ' + args.dir, file=sys.stderr)
    fast5_files = find_all_fast5s(args.dir)
    print('  Found ' + int_to_str(len(fast5_files)) + ' reads\n', file=sys.stderr)
    if not fast5_files:
        sys.exit()

    filters = ['file integrity']
    if args.min_length:
        filters.append('length >= ' + int_to_str(args.min_length) + ' bp')
    if args.min_mean_qual:
        filters.append('mean quality >= ' + str(args.min_mean_qual))
    if args.min_qual_window:
        filters.append('min window quality >= ' + str(args.min_qual_window))
    message = 'Filtering reads based on: ' + ', '.join(filters)
    print(message, file=sys.stderr)
    filtered_fast5_files = []
    total_bases = 0
    for fast5_file in fast5_files:
        passes, length = check_filters(fast5_file, args.min_length, args.min_mean_qual,
                                       args.min_qual_window, args.window_size)
        if passes:
            filtered_fast5_files.append(fast5_file)
            total_bases += length
    fast5_files = filtered_fast5_files
    print('  ' + int_to_str(len(fast5_files)) + ' reads remain after filtering (' +
          int_to_str(total_bases) + ' bp)\n', file=sys.stderr)
    if not fast5_files:
        sys.exit()

    good_fast5_files = set()
    if args.target_bases:
        print('Automatically setting a minimum window quality threshold in order to reach a target '
              'of ' + int_to_str(args.target_bases) + ' bp', file=sys.stderr)
        if total_bases < args.target_bases:
            print('Not enough total bases to reach target\n', file=sys.stderr)
        else:
            min_window_quals_and_lengths = sorted([min_window_qual_and_length(f, args.window_size)
                                                   for f in fast5_files], reverse=True)
            total_bases = 0
            min_window_qual_threshold = 0.0
            for min_window_qual, length, fast5_file in min_window_quals_and_lengths:
                total_bases += length
                min_window_qual_threshold = min_window_qual
                good_fast5_files.add(fast5_file)
                if total_bases > args.target_bases:
                    break
            fast5_files = [f for f in fast5_files if f in good_fast5_files]
            print('  min window quality threshold = ' + '%.2f' % min_window_qual_threshold,
                  file=sys.stderr)
            print('  ' + int_to_str(len(fast5_files)) + ' reads remain (' +
                  int_to_str(total_bases) + ' bp)\n', file=sys.stderr)

    print('Printing FASTQs to stdout', file=sys.stderr)
    for fast5_file in fast5_files:
        try:
            hdf5_file = h5py.File(fast5_file, 'r')
            names = get_hdf5_names(hdf5_file)
            basecall_location = get_best_fastq_hdf5_location(hdf5_file, names)
            if basecall_location:
                print(hdf5_file[basecall_location].value.decode(), end='', flush=True)
        except IOError:
            pass
    print('  Done!\n', file=sys.stderr)


def get_arguments():
    parser = argparse.ArgumentParser(description='FAST5 to FASTQ',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dir', type=str,
                        help='directory of FAST5 reads to extract (will be searched recursively)')
    parser.add_argument('--min_length', type=int, default=0,
                        help='Exclude reads shorter than this length (in bp)')
    parser.add_argument('--min_mean_qual', type=float, default=0.0,
                        help='Exclude reads with a mean qscore less than this value')
    parser.add_argument('--min_qual_window', type=float, default=0.0,
                        help='Exclude reads where their mean qscore in a sliding window drops '
                             'below this value')
    parser.add_argument('--window_size', type=int, default=50,
                        help='The size of the sliding window used for --min_qual_window')
    parser.add_argument('--target_bases', type=int, default=None,
                        help='If set, exclude the worst reads (as judged by their minimum qscore '
                             'in a sliding window) such that only this many bases remain')
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


def get_mean_score(hdf5_file, basecall_location):
    q = hdf5_file[basecall_location].value.decode().split('\n')[3]
    return statistics.mean([ord(c) - 33 for c in q])


def get_best_fastq_hdf5_location(hdf5_file, names):
    """
    This function returns the path in the FAST5 file to the best FASTQ. If there are multiple
    basecall locations, it returns the last one (hopefully from the most recent basecalling).
    """
    basecall_locations = sorted([x for x in names if x.upper().endswith('FASTQ')])
    two_d_locations = [x for x in basecall_locations if 'BASECALLED_2D' in x.upper()]
    template_locations = [x for x in basecall_locations if 'TEMPLATE' in x.upper()]
    complement_locations = [x for x in basecall_locations if 'COMPLEMENT' in x.upper()]

    # If the read has 2D basecalling, then that's what we use.
    if two_d_locations:
        return two_d_locations[-1]

    # If the read has both template and complement basecalling, then we choose the best based on
    # mean qscore.
    elif template_locations and complement_locations:
        template_location = template_locations[-1]
        complement_location = complement_locations[-1]
        mean_template_qscore = get_mean_score(hdf5_file, template_location)
        mean_complement_qscore = get_mean_score(hdf5_file, complement_location)
        if mean_template_qscore >= mean_complement_qscore:
            return template_location
        else:
            return complement_location

    # If the read has only template basecalling (normal for 1D) or only complement, then that's
    # what we use.
    elif template_locations:
        return template_locations[-1]
    elif complement_locations:
        return complement_locations[-1]

    # If the read has none of the above, but still has a fastq value in its hdf5, that's weird, but
    # we'll consider it a 1d read and use it.
    elif basecall_locations:
        return basecall_locations[-1]

    return None


def get_mean_qscore(quals):
    """
    Returns the mean qscore over the entire length of the qscore string.
    """
    try:
        return sum([ord(q) - 33 for q in quals]) / len(quals)
    except ZeroDivisionError:
        return 0.0


def get_min_window_qscore(quals, window_size):
    """
    Returns the minimum mean qscore over a sliding window.
    """
    quals = [ord(q) - 33 for q in quals]  # covert to numbers
    current_window_qscore = statistics.mean(quals[:window_size])
    shift_count = len(quals) - window_size
    if shift_count < 1:
        return current_window_qscore
    min_window_qscore = current_window_qscore
    for i in range(shift_count):
        leaving_window = quals[i]
        entering_window = quals[i + window_size]
        current_window_qscore += (entering_window - leaving_window) / window_size
        min_window_qscore = min(min_window_qscore, current_window_qscore)
    return min_window_qscore


def check_filters(fast5_file, min_length, min_mean_qual, min_qual_window, window_size):
    try:
        hdf5_file = h5py.File(fast5_file, 'r')
        names = get_hdf5_names(hdf5_file)
        basecall_location = get_best_fastq_hdf5_location(hdf5_file, names)
        if basecall_location:
            fastq_str = hdf5_file[basecall_location].value.decode()
            try:
                parts = fastq_str.split('\n')
                seq, quals = parts[1], parts[3]
            except IndexError:
                fastq_str, seq, quals = '', '', ''
            if not fastq_str or not seq:
                return False, 0
            if min_mean_qual and get_mean_qscore(quals) < min_mean_qual:
                return False, 0
            if min_length and len(seq) < min_length:
                return False, 0
            if min_qual_window and get_min_window_qscore(quals, window_size) < min_qual_window:
                return False, 0
            return True, len(seq)
    except (IOError, RuntimeError):
        pass
    return False, 0


def min_window_qual_and_length(fast5_file, window_size):
    try:
        hdf5_file = h5py.File(fast5_file, 'r')
        names = get_hdf5_names(hdf5_file)
        basecall_location = get_best_fastq_hdf5_location(hdf5_file, names)
        if basecall_location:
            fastq_str = hdf5_file[basecall_location].value.decode()
            try:
                parts = fastq_str.split('\n')
                seq, quals = parts[1], parts[3]
                return get_min_window_qscore(quals, window_size), len(seq), fast5_file
            except IndexError:
                pass
    except (IOError, RuntimeError):
        pass
    return 0.0, 0, fast5_file


def int_to_str(num):
    return '{:,}'.format(num)


if __name__ == '__main__':
    main()
