# FAST5 to FASTQ

This is just a simple script to extract FASTQ files from FAST5 files.

There are a number of other tools which can do this, including [Poretools](http://poretools.readthedocs.io/), [PoRe](https://sourceforge.net/p/rpore/wiki/Home/), [nanopolish extract](https://github.com/jts/nanopolish) and more.

I made this one for a couple of specific features:
* If there are multiple FASTQ groups in a FAST5 file (i.e. basecalling has been performed more than once), it extracts the last group. It does this on a per-read basis, so it's okay if some reads have only one group and others have more than one.
* Ability to filter using:
  * Read length (`--min_length`)
  * Mean Phred score (`--min_mean_qual`)
  * Minimum mean Phred score in a window, to exclude reads with low quality regions (`--min_qual_window`)
* Ability to automatically set `--min_qual_window` do get a target number of bases (`--target_bases`)


# Requirements

* Python 3.4 or later
* [h5py](https://github.com/h5py/h5py)


# Usage

Extracting all reads to FASTQ:
* `fast5_to_fastq.py path/to/fast5_directory > output.fastq
* This will search through the target directory recursively.

Gzip while you extract:
* `fast5_to_fastq.py path/to/fast5_directory | gzip > output.fastq

Filter based on length:
* `fast5_to_fastq.py --min_length 10000 path/to/fast5_directory | gzip > output.fastq
* Reads must be 10 kbp or longer to be included in the output.

Filter based on mean qscore:
* `fast5_to_fastq.py --min_mean_qual 11.5 path/to/fast5_directory | gzip > output.fastq
* Reads must have a mean qscore of at least 11.5 to be included in the output.

Filter based on min qscore over a sliding window:
* `fast5_to_fastq.py --min_qual_window 10.0 path/to/fast5_directory | gzip > output.fastq
* Reads must have a mean qscore over a sliding window (default 50 bp window, configurable with `--window_size`) that never drops below 10.0.

Aim for a target number of bases:
* `fast5_to_fastq.py --target_bases 100000000 path/to/fast5_directory | gzip > output.fastq
* Only outputs the best 100 Mbp of reads, as judged by their min qscore over a sliding window.
* Effectively sets `--min_qual_window` automatically to get the number of desired bases in the output.

How I (Ryan) like to use it:
* `fast5_to_fastq.py --min_length 2000 --target_bases 500000000 path/to/fast5_directory | gzip > output.fastq
* Since I mainly use Nanopore reads for bacterial isolate assembly, anything much over 100x depth is overkill, so I use `--target_bases` to aim for about 500 Mbp of reads.
* Repeat sequences of about 1 kbp are common in bacterial genomes (e.g. insertion sequences), so I use `--min_length` to exclude anything less than 2 kbp. That's large enough that I'll end with very few reads which are entirely contained within a repeat, but small enough that I'm not excluding small plasmid sequences.


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
