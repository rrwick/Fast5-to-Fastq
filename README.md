# Fast5 to Fastq

This is a simple script to extract FASTQ files from FAST5 files.

There are a number of other tools which can do this, including [Poretools](http://poretools.readthedocs.io/), [PoRe](https://sourceforge.net/p/rpore/wiki/Home/), [nanopolish extract](https://github.com/jts/nanopolish) and more. I made this one for a couple of specific features:
* If there are multiple FASTQ groups in a FAST5 file (i.e. basecalling has been performed more than once), it extracts the last group. It does this on a per-read basis, so it's okay if some reads have only one group and others have more than one.
* Ability to filter using:
  * Read length (`--min_length`)
  * Mean Phred score (`--min_mean_qual`)
  * Minimum mean Phred score in a window, to exclude reads with low quality regions (`--min_qual_window`)
* Ability to automatically set `--min_qual_window` to get a target number of bases (`--target_bases`)


# Requirements

* Python 3.4 or later
* [h5py](https://github.com/h5py/h5py)


# Installation

No installation is required - it's all just in one Python script:
```
git clone https://github.com/rrwick/Fast5-to-Fastq
Fast5-to-Fastq/fast5_to_fastq.py --help
```


# Usage

Extracting all reads from FAST5 to FASTQ:
* `fast5_to_fastq.py path/to/fast5_directory > output.fastq`
* This will search through the target directory recursively.

Gzip while you extract:
* `fast5_to_fastq.py path/to/fast5_directory | gzip > output.fastq.gz`

Filter based on length:
* `fast5_to_fastq.py --min_length 10000 path/to/fast5_directory | gzip > output.fastq.gz`
* To be included in the output, reads must be 10 kbp or longer.

Filter based on mean Phred quality score:
* `fast5_to_fastq.py --min_mean_qual 11.5 path/to/fast5_directory | gzip > output.fastq.gz`
* To be included in the output, reads must have a mean Phred score of at least 11.5.

Filter based on min Phred score over a sliding window:
* `fast5_to_fastq.py --min_qual_window 10.0 path/to/fast5_directory | gzip > output.fastq.gz`
* To be included in the output, reads must have a mean Phred score over a sliding window that never drops below 10.0.
* The default window size is 50 bp, but it's configurable with `--window_size`.

Aim for a target number of bases:
* `fast5_to_fastq.py --target_bases 100000000 path/to/fast5_directory | gzip > output.fastq.gz`
* Only outputs the best 100 Mbp of reads, as judged by their mininum mean Phred score over a sliding window.
* Effectively sets `--min_qual_window` automatically to get the number of desired bases in the output.

How I (Ryan) like to use it:
* `fast5_to_fastq.py --min_length 2000 --target_bases 500000000 path/to/fast5_directory | gzip > output.fastq.gz`
* I mainly use Nanopore reads for bacterial isolate assembly, and anything over 100x depth is probably overkill. So I use `--target_bases` to aim for about 500 Mbp of reads.
* Repeat sequences of about 1 kbp are common in bacterial genomes (e.g. insertion sequences), so I use `--min_length` to exclude anything less than 2 kbp. That's large enough that there should be very few reads which are entirely contained within a repeat (which aren't useful for assembly), but small enough that I'm not excluding small plasmid sequences.
* I'll then pass the output through [Porechop](https://github.com/rrwick/Porechop) to get rid of adapters and split/discard chimeric reads.


# FASTQ filtering

Albacore v1.1 and later can basecall directly to FASTQ, which is very nice! I therefore made a new version of this script (`fastq_to_fastq.py`) which takes a FASTQ as input instead of a FAST5 directory. This allows you to perform the various read filters (described above) when you did straight-to-FASTQ basecalling.

The usage is exactly the same as for `fast5_to_fastq.py`, just replace `path/to/fast5_directory` with `path/to/reads.fastq`:
* `fastq_to_fastq.py --min_length 2000 --target_bases 500000000 path/to/reads.fastq | gzip > output.fastq.gz`

Both `*.fastq` and `*.fastq.gz` should work as input formats.


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
