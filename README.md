<!--
---
output:
  html_document:
    keep_md: yes
---

-->

# *P-test* - a tool for estimating the rate of parallel molecular evolution at the genome-wide level


<!-- ########################################################################################################## -->
## Introduction
<!-- ########################################################################################################## -->

The *P-test* perl script allows users to compare ratios of synonymous and nonsynonymous parallel evolution in groups of four species. The main advantage of the approach is proper normalization of parallel synonymous and parallel nonsynonymous substitution counts. The method works on large samples of protein-coding sequences: it's best on hundreds or thousands of genes. The *P-test* approach is described fully in our paper, available in [Genome Biology and Evolution](https://academic.oup.com/gbe/article/doi/10.1093/gbe/evaa138/5870375?fbclid=IwAR03Zywh-So4xTERD5PAkqkE5YZE0KE3oJKysoz_P8kprNSHg5P20iiFsDQ).


***

<!-- ########################################################################################################## -->
## Main script
<!-- ########################################################################################################## -->

The P_test.pl script requirs perl 5 (tested on perl-5.30.1) and some specific (mainly BioPerl) libraries: Bio::AlignIO, Bio::SimpleAlign, Bio::TreeIO, Bio::SeqIO, List::Util.

***

<!-- ########################################################################################################## -->
## Example files
<!-- ########################################################################################################## -->

Other files are added to the repository for test run:
1) Each line of species_quartets.csv file contains comma delimited list of four species. Last common ancestor (LCA) of first two species and LCA of other two species should be younger, than LCA of all four species. In given tree topology script looks for substitutions, which appeared independently in path I (species 1 and 2) and in path II (species 3 and 4).
2) The tree.newick file contains phylogenetic tree of species, which gonna be analyzed (other species can be included in the tree too). The tree is used to calculate distance between LCAs of path I and path II.
3) The test_alignments.zip archive contains 800 fasta files with in-frame alignments of protein-coding sequences.
 

***

<!-- ########################################################################################################## -->
## Example run
<!-- ########################################################################################################## -->

P_test.pl script needs 6 arguments:

1) Line number in file with species quartets (each run of script handles one group of four species)
2) Name of file with species quartets
3) Name of file with phylogenetic tree
4) Path to gene alignments (that directory should include only nucleotide alignments in fasta format)
5) Path to two input files (i.e. to file with species quartets and to file with tree)
6) Path for output files

So the command should look like this:

perl P_test.pl 1 species_quartets.csv tree.newick /path/to/input/files/ /path/to/input/alignments/ /path/for/output/files


***

<!-- ########################################################################################################## -->
## Output
<!-- ########################################################################################################## -->

First three files contain information about parallel evolution:

P_test_with_threshold.csv

P_test_with_pseudocounts.csv

P_test_raw_counts.csv

Two of them contain P test values, calculated independently for 6 classes of mutations (AC, AT, AG, CT, CG, TG). As sometimes there are very few parallel substitutions (when little sample of genes is used or when species are closely related), P test could invoke division by zero. To avoid such cases, we use two solutions. In first one we add threshold value for number of parallel substitutions. If number of parallel substitutions is less than 3 for particular dinucleotide class, P test returns NA. In second approach we add pseudocounts.
The third file contains raw counts of different mutation patterns. It can be used, if you wish to recalculate some statistics.

Other three files contain analog of dn/ds test for two species from path II (species 3 and 4). Rate of synonymous substitutions is counted at 4-fold degenerate sites, while rate of nonsynonymous substitutions is based on non-degenerate sites. I added dn/ds to the output as a control: it should be always much lower than one, while P test could give even higher than one values.

dn_ds_analog_with_threshold.csv

dn_ds_analog_with_pseudocounts.csv

dn_ds_analog_raw_counts.csv

First two files contain dn/ds analog, calculated independently for 6 classes of mutations (AC, AT, AG, CT, CG, TG). Third file contains raw counts of different mutation patterns in path II.
