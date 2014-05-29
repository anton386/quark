# Quark - A method to rank and reveal viral quasispecies within Illumina data

## Prerequisites

Python Libraries
  * numpy
  * scipy
  * matplotlib
  * pysam
  * seaborn
  * networkx

## Usage

1. Top Ten and Top One
```
python quark.py [top_ten|top_one|rank_ten] <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -p --position <pos_of_interest>
  -i --increment <increment - default:10>
  -c --use-consensus-as-ref
  -f --export-consensus-to-fasta <filename>
  -v --export-variants-to-vcf <filename>
  -h --help
```

2. Getting Optimum Window
```
python quark.py windows <options>
  -b --bam
  -r --ref
  -w --window
  -p --position
  -i --increment
  -c --use-consensus-as-ref
  -h --help show this help screen
```

3. Variant Calling
```
python quark.py variant_calling <options>
  -b --bam
  -r --ref
  -p --position
  -i --increment
  -c --use-consensus-as-ref
  -v --export-variants-to-vcf
  -h --help show this help screen
```

4. Build Consensus
```
python quark.py consensus <options>
  -b --bam
  -r --ref
  -f --export-consensus-to-fasta
```

5. Haplotype Calling
```
python quark.py haplotype_calling <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -p --position <pos_of_interest>
  -i --increment <increment - default:10>
  -c --use-consensus-as-ref
  -f --export-consensus-to-fasta <filename>
  -v --export-variants-to-vcf <filename>
  -m --format <haplotype=0, csv=1, raw=2>
  -g --genbank <filename>
  -o --out <out-file>
  -h --help
```

6. Haplotype Entropy
```
python quark.py haplotype_entropy <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -s --start <start>
  -e --end <end>
  -c --use-consensus-as-ref
  -t --threads <threads>
  -o --out <out-file>
  -z --option <fixed_window|moving_window>
  -h --help
```