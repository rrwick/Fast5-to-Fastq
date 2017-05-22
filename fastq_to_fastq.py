#!/usr/bin/env python3

import os
import argparse
import statistics
import sys
import gzip


def main():
    args = get_arguments()

    print('\nLoading reads from: ' + args.input_fastq, file=sys.stderr)
    input_reads = load_fastq(args.input_fastq)
    print('  Found ' + int_to_str(len(input_reads)) + ' reads (' +
          int_to_str(sum(len(seq) for _, seq, _ in input_reads)) + ' bp)\n', file=sys.stderr)
    if not input_reads:
        sys.exit()

    filters = []
    if args.min_length:
        filters.append('length >= ' + int_to_str(args.min_length) + ' bp')
    if args.min_mean_qual:
        filters.append('mean quality >= ' + str(args.min_mean_qual))
    if args.min_qual_window:
        filters.append('min window quality >= ' + str(args.min_qual_window))
    message = 'Filtering reads based on: ' + ', '.join(filters)
    print(message, file=sys.stderr)
    filtered_reads = []
    total_bases = 0
    for name, seq, quals in input_reads:
        passes, length = check_filters(seq, quals, args.min_length, args.min_mean_qual,
                                       args.min_qual_window, args.window_size)
        if passes:
            filtered_reads.append((name, seq, quals))
            total_bases += length
    print('  ' + int_to_str(len(filtered_reads)) + ' reads remain after filtering (' +
          int_to_str(total_bases) + ' bp)\n', file=sys.stderr)
    if not filtered_reads:
        sys.exit()

    good_read_names = set()
    if args.target_bases:
        print('Automatically setting a minimum window quality threshold in order to reach a target '
              'of ' + int_to_str(args.target_bases) + ' bp', file=sys.stderr)
        if total_bases < args.target_bases:
            print('  not enough total bases to reach target\n', file=sys.stderr)
        else:
            min_window_quals_and_lengths = \
                sorted([min_window_qual_and_length(name, quals, args.window_size)
                        for name, _, quals in filtered_reads], reverse=True)
            total_bases = 0
            min_window_qual_threshold = 0.0
            for min_window_qual, length, read_name in min_window_quals_and_lengths:
                total_bases += length
                min_window_qual_threshold = min_window_qual
                good_read_names.add(read_name)
                if total_bases > args.target_bases:
                    break
            filtered_reads = [(name, seq, quals) for name, seq, quals in filtered_reads
                              if name in good_read_names]
            print('  min window quality threshold = ' + '%.2f' % min_window_qual_threshold,
                  file=sys.stderr)
            print('  ' + int_to_str(len(filtered_reads)) + ' reads remain (' +
                  int_to_str(total_bases) + ' bp)\n', file=sys.stderr)

    print('Printing FASTQs to stdout', file=sys.stderr)
    for name, seq, quals in filtered_reads:
        print(b'@' + name)
        print(seq)
        print(b'+')
        print(quals, flush=True)
    print('  Done!\n', file=sys.stderr)


def get_arguments():
    parser = argparse.ArgumentParser(description='FASTQ filter tool',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_fastq', type=str,
                        help='FASTQ file of reads to be filtered (can be gzipped)')
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

    args.input_fastq = os.path.abspath(args.input_fastq)
    if not os.path.isfile(args.input_fastq):
        sys.exit('Error: could not find ' + args.input_fastq)

    if args.min_length == 0 and args.min_mean_qual == 0.0 and args.min_qual_window == 0.0 and \
            args.target_bases is None:
        sys.exit('Error: no filters were used so this tool refuses to run (because the output\n'
                 '       FASTQ would be identical to the input FASTQ). Please use one of the\n'
                 '       following filters: --min_length, --min_mean_qual, --min_qual_window\n'
                 '       or --target_bases.')
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
    q = hdf5_file[basecall_location].value.split(b'\n')[3]
    return statistics.mean([c - 33 for c in q])


def get_mean_qscore(quals):
    """
    Returns the mean qscore over the entire length of the qscore string.
    """
    try:
        return sum([q - 33 for q in quals]) / len(quals)
    except ZeroDivisionError:
        return 0.0


def get_min_window_qscore(quals, window_size):
    """
    Returns the minimum mean qscore over a sliding window.
    """
    quals = [q - 33 for q in quals]  # covert to numbers
    current_window_qscore = statistics.mean(quals[:window_size])
    shift_count = len(quals) - window_size
    if shift_count < 1:
        return current_window_qscore
    min_window_qscore = current_window_qscore
    for i in range(shift_count):
        leaving_window = quals[i]
        entering_window = quals[i + window_size]
        current_window_qscore += (entering_window - leaving_window) / window_size
        if current_window_qscore < min_window_qscore:
            min_window_qscore = current_window_qscore
    return min_window_qscore


def check_filters(seq, quals, min_length, min_mean_qual, min_qual_window, window_size):
    if not seq:
        return False, 0
    if min_mean_qual and get_mean_qscore(quals) < min_mean_qual:
        return False, 0
    if min_length and len(seq) < min_length:
        return False, 0
    if min_qual_window and get_min_window_qscore(quals, window_size) < min_qual_window:
        return False, 0
    return True, len(seq)


def min_window_qual_and_length(name, quals, window_size):
    return get_min_window_qscore(quals, window_size), len(quals), name


def int_to_str(num):
    return '{:,}'.format(num)


def load_fastq(filename):
    reads = []
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    with open_func(filename, 'rb') as fastq:
        for line in fastq:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            if not stripped_line.startswith(b'@'):
                continue
            name = stripped_line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, qualities))
    return reads


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


if __name__ == '__main__':
    main()
