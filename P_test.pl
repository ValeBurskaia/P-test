#!/usr/bin/perl
use strict;
use warnings;

#use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::TreeIO;
use Bio::SeqIO;
use List::Util qw(sum);


#FILTRATIONS (could be switched off)

use constant GC_FILTR_ON     => 1;   #   
use constant CONSERV_FILTR_ON     => 1;  #


#PATHS TO FILES:

#1 specs_4_gammarids_sequenced_alone_no_duplications_FOR_FULL.csv misha_tree.newick /home/valya.burskaya/GAMMARUS_M/scripts/input/
#/home/valya.burskaya/GAMMARUS_M/eulimno_misha/cds_algn_4_2019/ /home/valya.burskaya/GAMMARUS_M/scripts/OUTPUT_NONPOP_PARALOGS_FULL_DATA/

#number of line in file with species quartets (one run of this script counts substitutions in one species quartet):
my $arr_id=$ARGV[0];

#name of file with species quartets:
my $species_groups = $ARGV[1];

#name of file with phylogenetic tree
my $tree_name = $ARGV[2];

#path to directory which contains file with species quartets and tree file
my $info_dir=$ARGV[3];

#path to directory which contains fasta files with alignments
my $rawdata_dir=$ARGV[4];

#path to directory for output files
my $output_dir=$ARGV[5];


my $species_names = "<" . $info_dir . $species_groups;
open (FILE_SPNAMES, $species_names) or die ("Could not open file");
my @array=<FILE_SPNAMES>;
close FILE_SPNAMES;
#my $line=$array[rand @array];

my $line=$array[$arr_id];

chomp($line);
#print "Random line is $line\n";
(my $sp_name_1, my $sp_name_2, my $sp_name_3, my $sp_name_4)=(split(",",$line))[0,1,2,3];



#avg_dist function measures distance between LCAs of species pair I and species pair II (and some other distances)

sub avg_dist{
	
	my $spec_1 = shift(@_);
	my $spec_2 = shift(@_);
	my $spec_3 = shift(@_);
	my $spec_4 = shift(@_);	
	
	
    #Reads the tree
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => "<" . $info_dir . $tree_name);
    my $tree = $treeio->next_tree;
	
	#Giving names to inner nodes
	my $num = 0;
	foreach my $node ($tree->get_nodes){
		if (not $node->is_Leaf){
			my $node_name="NODE_" . $num;
			$node->set_tag_value(-id,$node_name);
			$num++;
		}
	}	

	my $node1 = $tree->find_node(-id => $spec_1);
	my $node2 = $tree->find_node(-id => $spec_2);
	my $node3 = $tree->find_node(-id => $spec_3);
	my $node4 = $tree->find_node(-id => $spec_4);
	
	#distances themselves
	#lca
	my $lca_1_2 = $tree->get_lca($node1, $node2);
	my $lca_3_4 = $tree->get_lca($node3, $node4);
	my $distance_lca = $tree->distance($lca_1_2,$lca_3_4); 
	
	#out
	my $distance_1_3 = $tree->distance(-nodes => [$node1,$node3]); 
	my $distance_1_4 = $tree->distance(-nodes => [$node1,$node4]); 
	my $distance_2_3 = $tree->distance(-nodes => [$node2,$node3]); 
	my $distance_2_4 = $tree->distance(-nodes => [$node2,$node4]); 
	my $avg_distance_out = ($distance_1_3 +$distance_1_4+$distance_2_3+$distance_2_4)/4;
	
	#in
    my $distance_1_2 = $tree->distance(-nodes => [$node1,$node2]); #дистанция между парой узлов
	my $distance_3_4 = $tree->distance(-nodes => [$node3,$node4]); #дистанция между парой узлов
	my $avg_distance_in = ($distance_1_2 +$distance_3_4)/2;


	my @arr = ($avg_distance_in, $avg_distance_out, $distance_1_2, $distance_3_4, $distance_lca);
	return @arr;
	
}




#Below there are three functions for codon type classification 

#1)4-fold degenerate at 3rd position:

sub degenerate_codon_3pos{
	
	my %degenerate_3pos;
	
	foreach my $i (qw(A C T G)) {
		foreach my $j (qw(A C T G)) {
			my $codon_1 = Bio::Seq->new(
				-seq      => $i . $j . "A",
				-alphabet => 'dna'
			);
	
			my $aa1 = $codon_1->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_2 = Bio::Seq->new(
				-seq      => $i . $j . "C",
				-alphabet => 'dna'
			);
	
			my $aa2 = $codon_2->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_3 = Bio::Seq->new(
				-seq      => $i . $j . "T",
				-alphabet => 'dna'
			);
	
			my $aa3 = $codon_3->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_4 = Bio::Seq->new(
				-seq      => $i . $j . "G",
				-alphabet => 'dna'
			);

			my $aa4 = $codon_4->translate(
				-frame         => 0,
				-codontable_id => 1
			);
	
			if (   $aa1->seq eq $aa2->seq
				&& $aa1->seq eq $aa3->seq
				&& $aa1->seq eq $aa4->seq )
			{
				$degenerate_3pos{ $i . $j } = 0;
			}
	
		}
	}

	#foreach my $val ( keys %degenerate_3pos ) {
	#	print "degenerate 3-d position: "
	#	  . $val . "   "
	#	  . $degenerate_3pos{$val} . "\n";
	#}


return %degenerate_3pos;
	
}	

#2) nondegenerate at 1st position:
	
sub nondegenerate_codon_1pos{	
	
	my %nondegenerate_1pos;
	
	foreach my $i (qw(A C T G)) {
		foreach my $j (qw(A C T G)) {
			my $codon_1 = Bio::Seq->new(
				-seq      => "A" . $i . $j,
				-alphabet => 'dna'
			);
	
			my $aa1 = $codon_1->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_2 = Bio::Seq->new(
				-seq      => "C" . $i . $j,
				-alphabet => 'dna'
			);
	
			my $aa2 = $codon_2->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_3 = Bio::Seq->new(
				-seq      => "T" . $i . $j,
				-alphabet => 'dna'
			);
	
			my $aa3 = $codon_3->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_4 = Bio::Seq->new(
				-seq      => "G" . $i . $j,
				-alphabet => 'dna'
			);
	
			my $aa4 = $codon_4->translate(
				-frame         => 0,
				-codontable_id => 1
			);

			if (   $aa1->seq ne $aa2->seq
				&& $aa1->seq ne $aa3->seq
				&& $aa1->seq ne $aa4->seq
				&& $aa2->seq ne $aa3->seq
				&& $aa2->seq ne $aa4->seq
				&& $aa3->seq ne $aa4->seq )
			{
				$nondegenerate_1pos{ $i . $j } = 0;
			}
	
		}
	}

	delete $nondegenerate_1pos{AA};
	delete $nondegenerate_1pos{AG};
		
	#foreach my $val ( keys %nondegenerate_1pos ) {
	#	print "nondegenerate 1-st position: "
	#	  . $val . "   "
	#	  . $nondegenerate_1pos{$val} . "\n";
	#}

	return %nondegenerate_1pos;

}

#3) nondegenerate at 2nd position:

sub nondegenerate_codon_2pos{

	my %nondegenerate_2pos;
	
	foreach my $i (qw(A C T G)) {
		foreach my $j (qw(A C T G)) {
			my $codon_1 = Bio::Seq->new(
				-seq      => $i . "A" . $j,
				-alphabet => 'dna'
			);
	
			my $aa1 = $codon_1->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_2 = Bio::Seq->new(
				-seq      => $i . "C" . $j,
				-alphabet => 'dna'
			);
	
			my $aa2 = $codon_2->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_3 = Bio::Seq->new(
				-seq      => $i . "T" . $j,
				-alphabet => 'dna'
			);

			my $aa3 = $codon_3->translate(
				-frame         => 0,
				-codontable_id => 1
			);
			my $codon_4 = Bio::Seq->new(
				-seq      => $i . "G" . $j,
				-alphabet => 'dna'
			);

			my $aa4 = $codon_4->translate(
				-frame         => 0,
				-codontable_id => 1
			);

			if (   $aa1->seq ne $aa2->seq
				&& $aa1->seq ne $aa3->seq
				&& $aa1->seq ne $aa4->seq
				&& $aa2->seq ne $aa3->seq
				&& $aa2->seq ne $aa4->seq
				&& $aa3->seq ne $aa4->seq )
			{
				$nondegenerate_2pos{ $i . $j } = 0;
			}

		}
	}
	
	
	delete $nondegenerate_2pos{TG}; 

	#foreach my $val ( keys %nondegenerate_2pos ) {
	#	print "nondegenerate 2-nd position: "
	#	  . $val . "   "
	#	  . $nondegenerate_2pos{$val} . "\n";
	#}
	
	return %nondegenerate_2pos;

}



# The gap_filtration_2 function marks regions around gaps as bad (returns an array of 0s and 1s)


sub gap_filtration_2{
	
	my $prot_obj1 = shift(@_);
    my $prot_obj2 = shift(@_);
    my $prot_obj3 = shift(@_);
    my $prot_obj4 = shift(@_);

	my @gaps     = ();
	my @decision = ();

	#записываем гэпы

	foreach my $i ( 1 .. ( $prot_obj1->length ) ) {
		
		# handle every column of alignment, search for gaps
		my @column = (
			$prot_obj1->subseq( $i, $i ),
			$prot_obj2->subseq( $i, $i ),
			$prot_obj3->subseq( $i, $i ),
			$prot_obj4->subseq( $i, $i )
		);

		if ( grep( /[X\-]/, @column ) ) {
			push @gaps, "1";
		}
		else {
			push @gaps, "0";
		}
	}

	#we do not use the beginning of alignment anyway
	for ( my $i = 0 ; $i < 10 ; $i++ ) {
		push @decision, 0;
	}

	#handling of alignment
	foreach my $i ( 10 .. ( $#gaps - 10 ) ) {

		if (   $gaps[$i] == 0
			&& $gaps[ $i - 1 ] == 0
			&& $gaps[ $i - 2 ] == 0
			&& $gaps[ $i - 3 ] == 0
			&& $gaps[ $i - 4 ] == 0
			&& $gaps[ $i - 5 ] == 0
			&& $gaps[ $i - 6 ] == 0
			&& $gaps[ $i - 7 ] == 0
			&& $gaps[ $i - 8 ] == 0
			&& $gaps[ $i - 9 ] == 0
			&& $gaps[ $i - 10 ] == 0
			&& $gaps[ $i + 1 ] == 0
			&& $gaps[ $i + 2 ] == 0
			&& $gaps[ $i + 3 ] == 0
			&& $gaps[ $i + 4 ] == 0
			&& $gaps[ $i + 5 ] == 0
			&& $gaps[ $i + 6 ] == 0
			&& $gaps[ $i + 7 ] == 0
			&& $gaps[ $i + 8 ] == 0
			&& $gaps[ $i + 9 ] == 0
			&& $gaps[ $i + 10 ] == 0)
		{
			push @decision, 1;
		}
		else {
			push @decision, 0;
		}
	}

	#we do not use the end of alignment anyway
	for ( my $i = 0 ; $i <= 10 ; $i++ ) {
		push @decision, 0;
	}
	
	return @decision;
}


# The function conserv_filtration_2 marks sites, which are surrounded by highly polymorphic sites, as bad
# (returns an array of 0s and 1s)
 
sub conserv_filtration_2{

	my $seq_obj1 = shift(@_);
	my $seq_obj2 = shift(@_);
	my $seq_obj3 = shift(@_);
	my $seq_obj4 = shift(@_);

	my @seq_str=($seq_obj1->seq,$seq_obj2->seq,$seq_obj3->seq,$seq_obj4->seq);
	my @seq_AA=();
	my @seq_str_AA=();

	foreach my $seq_num (0 .. 3){
				my $temp_seq = Bio::Seq->new(
					-seq      => $seq_str[$seq_num],
					-alphabet => 'dna'
					);
				$seq_AA[$seq_num]=$temp_seq->translate(
					-frame         => 0,
					-codontable_id => 1
					);
				$seq_str_AA[$seq_num]=$seq_AA[$seq_num]->seq;
	}
		
	my $use_flank=10;

	my @cons_AA=();

	foreach my $i ( 0 .. ( $seq_AA[0]->length-1 ) ) { 
		my $first_AA=substr($seq_str_AA[0],$i,1);
		my $cons=1;
			foreach my $seq_num (1 .. 3){
				if(substr($seq_str_AA[$seq_num],$i,1) ne $first_AA){
					$cons=0;
					last;
				}
			}
		push(@cons_AA,$cons);
	}

	my @conserv_decision = ();

	foreach my $i ( 0 .. ( $seq_AA[0]->length - 1) ) {

		if ( $i < $use_flank || $i > ($seq_AA[0]->length-$use_flank-1) ) {
			push @conserv_decision, "0"; #будет негодный, т.к. лень обрабатывать
		}
		else{
			my $ambigs=0;
			my $mismatches=$use_flank*2-sum(@cons_AA[(($i-$use_flank .. $i-1),($i+1 .. $i+$use_flank))]);

			if($mismatches<3){
				push @conserv_decision, "1";    # good AA position
			}else{
				push @conserv_decision, "0";    # bad
			}
		}
			
			
			
	}

	return @conserv_decision;	
}





# The GC_filtration function marks sites, which are precided by C and followed by G as bad

sub GC_filtration{

	my $seq_obj1 = shift(@_);
    my $seq_obj2 = shift(@_);
    my $seq_obj3 = shift(@_);
    my $seq_obj4 = shift(@_);

    my @seq_str=($seq_obj1->seq,$seq_obj2->seq,$seq_obj3->seq,$seq_obj4->seq);

	my @seq_GC;
	push @seq_GC,0;

	foreach my $i ( 1 .. ( length($seq_str[0])-2 ) ) { 
		my $GC=1;
			foreach my $seq_num (0 .. 3){
				my $prev_NT=substr($seq_str[$seq_num],$i-1,1); 
				my $next_NT=substr($seq_str[$seq_num],$i+1,1); 
				if(uc($prev_NT) eq "C" && uc($next_NT) eq "G" ){
					$GC=0;
					last;
				}
			}
		push(@seq_GC,$GC);
	}
	push @seq_GC,0;

	return @seq_GC;
}


# The polymorph_type function is needed for classification of substitution types, as in P test
# different types of substitutions are counted separately

sub polymorph_type{
	my $nuc1 = shift(@_);
	my $nuc2 = shift(@_);
	
	
	#six types of substitutions: AC AT AG CT CG TC
					if (   ( $nuc1 eq "A" && $nuc2 eq "C" ) || ( $nuc1 eq "C" && $nuc2 eq "A" ) )
					{
						return "AC";
					}
					if (   ( $nuc1 eq "A" && $nuc2 eq "T" ) || ( $nuc1 eq "T" && $nuc2 eq "A" ) )
					{
						return "AT";
					}
					if (   ( $nuc1 eq "A" && $nuc2 eq "G" )	|| ( $nuc1 eq "G" && $nuc2 eq "A" ) )
					{
						return "AG";
					}
					if (   ( $nuc1 eq "C" && $nuc2 eq "T" )	|| ( $nuc1 eq "T" && $nuc2 eq "C" ) )
					{
						return "CT";
					}
					if (   ( $nuc1 eq "C" && $nuc2 eq "G" )	|| ( $nuc1 eq "G" && $nuc2 eq "C" ) )
					{
						return "CG";
					}
					if (   ( $nuc1 eq "T" && $nuc2 eq "G" ) || ( $nuc1 eq "G" && $nuc2 eq "T" ) )
					{
						return "TG";
					}
					if ( $nuc1 eq "A" && $nuc2 eq "A" ) {
						return "AA";
					}
					if ( $nuc1 eq "C" && $nuc2 eq "C" ) {
						return "CC";
					}
					if ( $nuc1 eq "T" && $nuc2 eq "T" ) {
						return "TT";
					}
					if ( $nuc1 eq "G" && $nuc2 eq "G" ) {
						return "GG";
					}
}





################################################################################################
# THE MAIN BODY OF THE SCRIPT: 
#      ..                    ..                    ..                    ..   
#     00\\\\\\\\~~~~~~~     00\\\\\\\\~~~~~~~     00\\\\\\\\~~~~~~~     00\\\\\\\\~~~~~~~
#     "  ||   ||            "  ||   ||            "  ||   ||            "  ||   ||         



opendir( RAWDIR, $rawdata_dir ) or die $!;

#initialisations

#1)For parallel substitutions count

#path 1:

my %parall_path1_1pos;

$parall_path1_1pos{AC} = 0;
$parall_path1_1pos{AT} = 0;
$parall_path1_1pos{AG} = 0;
$parall_path1_1pos{CT} = 0;
$parall_path1_1pos{CG} = 0;
$parall_path1_1pos{TG} = 0;

my %parall_path1_2pos;

$parall_path1_2pos{AC} = 0;
$parall_path1_2pos{AT} = 0;
$parall_path1_2pos{AG} = 0;
$parall_path1_2pos{CT} = 0;
$parall_path1_2pos{CG} = 0;
$parall_path1_2pos{TG} = 0;

my %parall_path1_3pos;

$parall_path1_3pos{AC} = 0;
$parall_path1_3pos{AT} = 0;
$parall_path1_3pos{AG} = 0;
$parall_path1_3pos{CT} = 0;
$parall_path1_3pos{CG} = 0;
$parall_path1_3pos{TG} = 0;



#path 2(parall):

#1st pos

my %parall_path2_1pos;

$parall_path2_1pos{AT_AT} = 0;
$parall_path2_1pos{AT_AA} = 0;
$parall_path2_1pos{AT_TT} = 0;

$parall_path2_1pos{AC_AC} = 0;
$parall_path2_1pos{AC_AA} = 0;
$parall_path2_1pos{AC_CC} = 0;

$parall_path2_1pos{AG_AG} = 0;
$parall_path2_1pos{AG_AA} = 0;
$parall_path2_1pos{AG_GG} = 0;

$parall_path2_1pos{CT_CT} = 0;
$parall_path2_1pos{CT_CC} = 0;
$parall_path2_1pos{CT_TT} = 0;

$parall_path2_1pos{CG_CG} = 0;
$parall_path2_1pos{CG_CC} = 0;
$parall_path2_1pos{CG_GG} = 0;

$parall_path2_1pos{TG_TG} = 0;
$parall_path2_1pos{TG_TT} = 0;
$parall_path2_1pos{TG_GG} = 0;


#2nd pos

my %parall_path2_2pos;

$parall_path2_2pos{AT_AT} = 0;
$parall_path2_2pos{AT_AA} = 0;
$parall_path2_2pos{AT_TT} = 0;

$parall_path2_2pos{AC_AC} = 0;
$parall_path2_2pos{AC_AA} = 0;
$parall_path2_2pos{AC_CC} = 0;

$parall_path2_2pos{AG_AG} = 0;
$parall_path2_2pos{AG_AA} = 0;
$parall_path2_2pos{AG_GG} = 0;

$parall_path2_2pos{CT_CT} = 0;
$parall_path2_2pos{CT_CC} = 0;
$parall_path2_2pos{CT_TT} = 0;

$parall_path2_2pos{CG_CG} = 0;
$parall_path2_2pos{CG_CC} = 0;
$parall_path2_2pos{CG_GG} = 0;

$parall_path2_2pos{TG_TG} = 0;
$parall_path2_2pos{TG_TT} = 0;
$parall_path2_2pos{TG_GG} = 0;


#3rd pos 

my %parall_path2_3pos;

$parall_path2_3pos{AT_AT} = 0;
$parall_path2_3pos{AT_AA} = 0;
$parall_path2_3pos{AT_TT} = 0;

$parall_path2_3pos{AC_AC} = 0;
$parall_path2_3pos{AC_AA} = 0;
$parall_path2_3pos{AC_CC} = 0;

$parall_path2_3pos{AG_AG} = 0;
$parall_path2_3pos{AG_AA} = 0;
$parall_path2_3pos{AG_GG} = 0;

$parall_path2_3pos{CT_CT} = 0;
$parall_path2_3pos{CT_CC} = 0;
$parall_path2_3pos{CT_TT} = 0;

$parall_path2_3pos{CG_CG} = 0;
$parall_path2_3pos{CG_CC} = 0;
$parall_path2_3pos{CG_GG} = 0;

$parall_path2_3pos{TG_TG} = 0;
$parall_path2_3pos{TG_TT} = 0;
$parall_path2_3pos{TG_GG} = 0;




#overall (path 2):

#1pos

my %overall_1pos;

$overall_1pos{AC} = 0;
$overall_1pos{AT} = 0;
$overall_1pos{AG} = 0;
$overall_1pos{CT} = 0;
$overall_1pos{CG} = 0;
$overall_1pos{TG} = 0;

$overall_1pos{AA} = 0;
$overall_1pos{CC} = 0;
$overall_1pos{TT} = 0;
$overall_1pos{GG} = 0;

#2pos

my %overall_2pos;

$overall_2pos{AC} = 0;
$overall_2pos{AT} = 0;
$overall_2pos{AG} = 0;
$overall_2pos{CT} = 0;
$overall_2pos{CG} = 0;
$overall_2pos{TG} = 0;

$overall_2pos{AA} = 0;
$overall_2pos{CC} = 0;
$overall_2pos{TT} = 0;
$overall_2pos{GG} = 0;

#3pos

my %overall_3pos;

$overall_3pos{AC} = 0;
$overall_3pos{AT} = 0;
$overall_3pos{AG} = 0;
$overall_3pos{CT} = 0;
$overall_3pos{CG} = 0;
$overall_3pos{TG} = 0;

$overall_3pos{AA} = 0;
$overall_3pos{CC} = 0;
$overall_3pos{TT} = 0;
$overall_3pos{GG} = 0;



#hash keeps info about degenerate/nondegenerate codons
	
my %nondegenerate_1pos = nondegenerate_codon_1pos();
my %nondegenerate_2pos = nondegenerate_codon_2pos();
my %degenerate_3pos = degenerate_codon_3pos();





#Alignments handling

#Reads all files in the dir (it should only contain aligned genes)

while ( my $file = readdir(RAWDIR) ) {

	if ( $file =~ /^\..*/ ) { next; }    # Skip everything, that begins with dot

	#print $file . "\n";

	#reads the alignment
	my $alignment = Bio::AlignIO->new(
		-file        => '<' . $rawdata_dir . $file,
		-format      => 'fasta',
		-alphabet => "dna"
	);

	my $alnObj = $alignment->next_aln;

	#uppercase
	$alnObj->uppercase;

    #reads sequences of four target species:

	my $seq_obj1 = Bio::Seq->new();
	my $seq_obj2 = Bio::Seq->new();
	my $seq_obj3 = Bio::Seq->new();
	my $seq_obj4 = Bio::Seq->new();

	my $count_seqs = 0; #counter for number of specs
	my @tmp_seq_array = ();
	
	foreach my $seq_obj ( $alnObj->each_seq ) {

		my $seq_name = $seq_obj->display_id;

		#print $seq_name, "_\n";
		if ( $seq_name eq $sp_name_1 ) {
			$seq_obj1 = $seq_obj;    #path 1
			$count_seqs++;
			push @tmp_seq_array, $seq_obj1;
		}
		if ( $seq_name eq $sp_name_2 ) {
			$seq_obj2 = $seq_obj;    #path 1
			$count_seqs++;
			push @tmp_seq_array, $seq_obj2;
		}
		if ( $seq_name eq $sp_name_3 ) {
			$seq_obj3 = $seq_obj;    #path 2
			$count_seqs++;
			push @tmp_seq_array, $seq_obj3;
		}
		if ( $seq_name eq $sp_name_4 ) {
			$seq_obj4 = $seq_obj;    #path 2
			$count_seqs++;
			push @tmp_seq_array, $seq_obj4;
		}
	}
	
	#here we skip the gene, if it contains stops of too high proportion of gaps
	my $count_good_genes=0;
	
	if($count_seqs!=4){
		next;
	}elsif($count_seqs==4){

		foreach my $seq_from_arr (@tmp_seq_array){
	
			my $seq=$seq_from_arr->seq;
			$seq = uc($seq);
		
			my $gene_stops_existence =0;
			
			#stops   	     
			$gene_stops_existence = $seq =~ /^([ATGCN-]{3})*(TGA|TAA|TAG).*[ATGCN-]{3}/;
          
          
			#gaps
			my $gene_gap = $seq =~ tr/-//;
			#gene length
			my $gene_len = length($seq);
			#gap proportion 
			my $gap_proportion=$gene_gap/$gene_len;
			my $too_many_gaps=0;
			if($gap_proportion>0.7){
				$too_many_gaps=1;
			}
			
			#control point
			if ($too_many_gaps==0 &&  $gene_stops_existence ==0){ 
			      $count_good_genes++;
			} 
		}
		
	}
	
	if($count_good_genes!=4){
		next;
	}
	

	#convertion to amino acid sequences:
	my $prot_obj1 = $seq_obj1->translate(
		-frame         => 0,
		-codontable_id => 1
	);
	my $prot_obj2 = $seq_obj2->translate(
		-frame         => 0,
		-codontable_id => 1
	);
	my $prot_obj3 = $seq_obj3->translate(
		-frame         => 0,
		-codontable_id => 1
	);
	my $prot_obj4 = $seq_obj4->translate(
		-frame         => 0,
		-codontable_id => 1
	);

	#arrays for filtration, calling of filtering functions
	
	my @GC_dec=();
	my @conserv_dec=();

	if(GC_FILTR_ON){
		@GC_dec=GC_filtration($seq_obj1, $seq_obj2, $seq_obj3, $seq_obj4);
	}else{
		foreach my $i ( 1 .. ( $seq_obj1->length ) ) {
			push @GC_dec, "1";	
		}
	}
	
	
	if(CONSERV_FILTR_ON){
		@conserv_dec=conserv_filtration_2($seq_obj1, $seq_obj2, $seq_obj3, $seq_obj4);
	}else{
		foreach my $i ( 1 .. ( $seq_obj1->length ) ) {
			push @conserv_dec, "1";	
		}
	}
	
	my @gap_dec=gap_filtration_2($prot_obj1,$prot_obj2,$prot_obj3,$prot_obj4);
	

	#CODON HANDLING (main part of script)
	#we read alignment codon by codon
	#$i - codon number, 1 - based
	foreach my $i ( 1 .. ( ( $seq_obj1->length ) / 3 ) ) {

		#reads each seq
		
		#amino acids
		my $prot_subseq1 = $prot_obj1->subseq( $i, $i );
		my $prot_subseq2 = $prot_obj2->subseq( $i, $i );
		my $prot_subseq3 = $prot_obj3->subseq( $i, $i );
		my $prot_subseq4 = $prot_obj4->subseq( $i, $i );
		#codons
		my $nuc_subseq1 = $seq_obj1->subseq(( $i * 3 ) - 2, $i * 3);
		my $nuc_subseq2 = $seq_obj2->subseq(( $i * 3 ) - 2, $i * 3);
		my $nuc_subseq3 = $seq_obj3->subseq(( $i * 3 ) - 2, $i * 3);
		my $nuc_subseq4 = $seq_obj4->subseq(( $i * 3 ) - 2, $i * 3);
        
        #nucleotides of path I, separated:
        #species 1 and 2
		my $path_1_nucA_1 = substr $nuc_subseq1, 0, 1;
		my $path_1_nucB_1 = substr $nuc_subseq2, 0, 1;
		my $path_1_nucA_2 = substr $nuc_subseq1, 1, 1;
		my $path_1_nucB_2 = substr $nuc_subseq2, 1, 1;
		my $path_1_nucA_3 = substr $nuc_subseq1, 2, 1;
		my $path_1_nucB_3 = substr $nuc_subseq2, 2, 1;
		# nucleotides of path II, separated:
		# species 3 and 4
		my $path_2_nucA_1 = substr $nuc_subseq3, 0, 1;
		my $path_2_nucB_1 = substr $nuc_subseq4, 0, 1;
		my $path_2_nucA_2 = substr $nuc_subseq3, 1, 1;
		my $path_2_nucB_2 = substr $nuc_subseq4, 1, 1;
		my $path_2_nucA_3 = substr $nuc_subseq3, 2, 1;
		my $path_2_nucB_3 = substr $nuc_subseq4, 2, 1;


		#SKIP CODON, IF IT IS CLOSE TO GAP
		if($gap_dec[$i] == 0) {
			next;
		}
		
		
		
		#two of three positions in the codon should be fully conservative, lets check it:
		
		my $pos1_cons=0;
		my $pos2_cons=0;
		my $pos3_cons=0;
		
		
		
		# 1st pos consevatism check
		if(( $path_1_nucA_1 eq $path_1_nucB_1 )
		&& ( $path_1_nucA_1 eq $path_2_nucA_1 ) && ($path_1_nucA_1 eq $path_2_nucB_1 )){
			$pos1_cons=1;
		}
		
		# 2nd pos consevatism check
		if(( $path_1_nucA_2 eq $path_1_nucB_2 )
		&& ( $path_1_nucA_2 eq $path_2_nucA_2 ) && ($path_1_nucA_2 eq $path_2_nucB_2 )){
			$pos2_cons=1;
		}
		
		# 3rd pos consevatism check
		if(($path_1_nucA_3 eq $path_1_nucB_3 )
		&& ( $path_1_nucA_3 eq $path_2_nucA_3 ) && ($path_1_nucA_3 eq $path_2_nucB_3 )){
			$pos3_cons=1;
		}
		    
		#PARALLEL SUBSTITUTIONS COUNT (for P test)

		#POS_1 - count parallel and nonparallel changes in 1st pos of the codon
		
		#check that 2nd and 3rd positions are conservative in four specs:
		if($pos2_cons==1 && $pos3_cons==1){
		    
		    #only nondegenerate sites
			my $dinuc = $path_1_nucA_2 . $path_1_nucA_3;
			if (exists $nondegenerate_1pos{$dinuc}){
				#site should not be in CpG context and it should be surrounded by more or less conservative region:
				if ( $GC_dec[ ( $i * 3 ) - 3 ] == 1
					&& $conserv_dec[$i-1] == 1 ){
		
					#1)Different nucleotides in 1st position of path I
					if ($path_1_nucA_1 ne $path_1_nucB_1){
						my $polymorph_path_1 = polymorph_type($path_1_nucA_1,$path_1_nucB_1); #orders a pair of nucleotides
						$parall_path1_1pos{$polymorph_path_1}++; #учет

						#2)look at path II:
						
					    my $polymorph_path_2 = polymorph_type($path_2_nucA_1,$path_2_nucB_1);  #orders a pair of nucleotides
						my $hash_index=$polymorph_path_1 . "_" . $polymorph_path_2;   #strange index
						
						#a)if in path I we see the same substitution as in path II - ABAB or ABBA:
						if ($polymorph_path_1 eq $polymorph_path_2){  
							$parall_path2_1pos{$hash_index}++;    #count

						}
						#b)If there is no substitution in path II (ABAA or ABBB)
						elsif((($path_1_nucA_1 eq $path_2_nucA_1) && ($path_1_nucA_1 eq $path_2_nucB_1))||
						(($path_1_nucB_1 eq $path_2_nucA_1) && ($path_1_nucB_1 eq $path_2_nucB_1))){
							$parall_path2_1pos{$hash_index}++;    #учет
							
						}
					}
				}
			}
		}

		#POS_2
		
		#check that 1st and 3rd positions are conservative in four specs:
		if($pos1_cons==1 && $pos3_cons==1){
		    
		    #only non-degenerate sites
			my $dinuc = $path_1_nucA_1 . $path_1_nucA_3;
			if (exists $nondegenerate_2pos{$dinuc}){
				#site should not be in CpG context and it should be surrounded by more or less conservative region:
				if (   $GC_dec[ ( $i * 3 ) - 2 ] == 1
					&& $conserv_dec[$i-1] == 1 ){
					
					#1)Different nucleotides at 2nd position of path I
					if ($path_1_nucA_2 ne $path_1_nucB_2){
						my $polymorph_path_1 = polymorph_type($path_1_nucA_2,$path_1_nucB_2);  #orders a pair of nucleotides
						$parall_path1_2pos{$polymorph_path_1}++;    #count
						
						#2)look at path II:
					    my $polymorph_path_2 = polymorph_type($path_2_nucA_2,$path_2_nucB_2);
						my $hash_index=$polymorph_path_1 . "_" . $polymorph_path_2;
						
						#a)if in path I we see the same substitution as in path II - ABAB or ABBA:
						if ($polymorph_path_1 eq $polymorph_path_2){  
							$parall_path2_2pos{$hash_index}++;
						}
						#b)If there is no substitution in path II (ABAA or ABBB)
						elsif((($path_1_nucA_2 eq $path_2_nucA_2) && ($path_1_nucA_2 eq $path_2_nucB_2))||
						(($path_1_nucB_2 eq $path_2_nucA_2) && ($path_1_nucB_2 eq $path_2_nucB_2))){
							$parall_path2_2pos{$hash_index}++;
						}
					}
				}
			}
		}


		#POS_3
		
		#check that 1st and 3rd positions are conservative in four specs:
		if($pos1_cons==1 && $pos2_cons==1){
		    
		    #only 4-fold degenerate sites
			my $dinuc = $path_1_nucA_1 . $path_1_nucA_2;
			if (exists $degenerate_3pos{$dinuc}){
				#site should not be in CpG context and it should be surrounded by more or less conservative region:
				if (   $GC_dec[ ( $i * 3 ) - 1 ] == 1
					&& $conserv_dec[$i-1] == 1 ){
		
					#1)Different nucleotides at 3rd position of path I
					if ($path_1_nucA_3 ne $path_1_nucB_3){
						my $polymorph_path_1 = polymorph_type($path_1_nucA_3,$path_1_nucB_3);
						$parall_path1_3pos{$polymorph_path_1}++;

						#2)look at path II:
						
					    my $polymorph_path_2 = polymorph_type($path_2_nucA_3,$path_2_nucB_3);
						my $hash_index=$polymorph_path_1 . "_" . $polymorph_path_2;
						
						#a)if in path I we see the same substitution as in path II - ABAB or ABBA:
						if ($polymorph_path_1 eq $polymorph_path_2){  
							$parall_path2_3pos{$hash_index}++;
						}
						#b)If there is no substitution in path II (ABAA or ABBB)
						elsif((($path_1_nucA_3 eq $path_2_nucA_3) && ($path_1_nucA_3 eq $path_2_nucB_3))||
						(($path_1_nucB_3 eq $path_2_nucA_3) && ($path_1_nucB_3 eq $path_2_nucB_3))){
							$parall_path2_3pos{$hash_index}++;
						}
					}
				}
			}
		}



		# OVERALL SUBSTITUTIONS COUNT (for trivial dn/ds test analog), based on path II:

		#POS_1

		# 2nd and 3rd positions should be conservative:
		if ($pos2_cons==1 && $pos3_cons==1){
        	# only non-degenerate sites
			my $dinuc = $path_2_nucA_2 . $path_2_nucA_3;
			if (exists $nondegenerate_1pos{$dinuc}) {
				# site should not be in CpG context and it should be surrounded by more or less conservative region:
				if ($GC_dec[ ( $i * 3 ) - 3 ] == 1 && $conserv_dec[$i-1] == 1){
                #if (   $GC_dec[ ( $i * 3 ) - 3 ] == 1){
					my $polymorph = polymorph_type($path_2_nucA_1, $path_2_nucB_1);
					$overall_1pos{$polymorph}++;
				}
			}
		}

		#POS_2

		# 1st and 3rd positions should be conservative:
		if ($pos1_cons==1 && $pos3_cons==1){
			# only non-degenerate sites
			my $dinuc = $path_2_nucA_1 . $path_2_nucA_3;
			if (exists $nondegenerate_2pos{$dinuc}) {
				#site should not be in CpG context and it should be surrounded by more or less conservative region:
				if ($GC_dec[ ( $i * 3 ) - 2 ] == 1 && $conserv_dec[$i-1] == 1 ){
                #if (   $GC_dec[ ( $i * 3 ) - 2 ] == 1){
                	
					my $polymorph = polymorph_type($path_2_nucA_2, $path_2_nucB_2);
					$overall_2pos{$polymorph}++;
				}
			}
		}

		#POS_3
		
		# 1st and 2nd positions should be conservative:
		if ($pos1_cons==1 && $pos2_cons==1){
   			# only 4-fold degenerate sites
			my $dinuc = $path_2_nucA_1 . $path_2_nucA_2;
			if (exists $degenerate_3pos{$dinuc}) {
				#site should not be in CpG context and it should be surrounded by more or less conservative region:
				if (   $GC_dec[ ( $i * 3 ) - 1 ] == 1 && $conserv_dec[$i-1] == 1 ){
					my $polymorph = polymorph_type($path_2_nucA_3, $path_2_nucB_3);
					$overall_3pos{$polymorph}++;
				}
			}
			    
		}
		
	}# the end of the file handling

}#the end of the dir handling



#OUTPUT:

my $filename =
"P_test_raw_counts_" 
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";
  
my $filename2 =
"P_test_with_threshold_"
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";

 my $filename3 =
 "P_test_with_pseudocounts_"
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";


my $file = ">>" . $output_dir . $filename;
open( FILE_INFO, $file ) or die("Could not open file");

my $file2 = ">>" . $output_dir . $filename2;
open( FILE_INFO_2, $file2 ) or die("Could not open file");

my $file3 = ">>" . $output_dir . $filename3;
open( FILE_INFO_3, $file3 ) or die("Could not open file");


#0) output about distances between last common ancestors of species pair I and species pair II

my @arr=avg_dist($sp_name_1,$sp_name_2,$sp_name_3,$sp_name_4);
my $in_dist=$arr[0]; 
my $out_dist=$arr[1];
my $path1_dist=$arr[2];
my $path2_dist=$arr[3];
my $lca_dist=$arr[4];



print FILE_INFO "sp_name_1    sp_name_2    sp_name_3     sp_name_4     LCA_dist     AC_AC_1pos     AC_AA_1pos     AC_CC_1pos     AC_AC_2pos     AC_AA_2pos     AC_CC_2pos     AC_AC_3pos     AC_AA_3pos     AC_CC_3pos     AG_AG_1pos     AG_AA_1pos     AG_GG_1pos     AG_AG_2pos     AG_AA_2pos     AG_GG_2pos     AG_AG_3pos     AG_AA_3pos     AG_GG_3pos     AT_AT_1pos     AT_AA_1pos     AT_TT_1pos     AT_AT_2pos     AT_AA_2pos     AT_TT_2pos     AT_AT_3pos     AT_AA_3pos     AT_TT_3pos     CG_CG_1pos     CG_CC_1pos     CG_GG_1pos     CG_CG_2pos     CG_CC_2pos     CG_GG_2pos     CG_CG_3pos     CG_CC_3pos     CG_GG_3pos     CT_CT_1pos     CT_CC_1pos     CT_TT_1pos     CT_CT_2pos     CT_CC_2pos     CT_TT_2pos     CT_CT_3pos     CT_CC_3pos     CT_TT_3pos     TG_TG_1pos     TG_TT_1pos     TG_GG_1pos     TG_TG_2pos     TG_TT_2pos     TG_GG_2pos     TG_TG_3pos     TG_TT_3pos     TG_GG_3pos\n";
print FILE_INFO "$sp_name_1    $sp_name_2    $sp_name_3     $sp_name_4     $lca_dist     ";


print FILE_INFO_2 "sp_name_1     sp_name_2     sp_name_3     $sp_name_4     LCA_dist     AC     AG     AT     CG     CT     TG\n";
print FILE_INFO_2 "$sp_name_1     $sp_name_2     $sp_name_3     $sp_name_4     $lca_dist     ";


print FILE_INFO_3 "sp_name_1     sp_name_2     sp_name_3     $sp_name_4     LCA_dist     AC     AG     AT     CG     CT     TG\n";
print FILE_INFO_3 "$sp_name_1     $sp_name_2     $sp_name_3     $sp_name_4     $lca_dist     ";

#1) output about parallel substitutions (for P test)



foreach my $val (sort keys %parall_path1_1pos) {

	#RAW COUNTS OF DIFFERENT NUCLEOTIDE PATTERNS:

	(my $first_nuc, my $second_nuc)=(split("",$val))[0,1];
	my $homo_1=$val . "_" . $first_nuc . $first_nuc;
	my $homo_2=$val . "_" . $second_nuc . $second_nuc;
	my $heter=$val . "_" . $val;
	
	print FILE_INFO "$parall_path2_1pos{$heter}     $parall_path2_1pos{$homo_1}     $parall_path2_1pos{$homo_2}    ";
	print FILE_INFO "$parall_path2_2pos{$heter}     $parall_path2_2pos{$homo_1}     $parall_path2_2pos{$homo_2}    ";
	print FILE_INFO "$parall_path2_3pos{$heter}     $parall_path2_3pos{$homo_1}     $parall_path2_3pos{$homo_2}    ";

	


	#RATIOS WITH THRESHOLD:

	my $threshold = 2;
	
	my $p_test;

	#dSp

	my $dsp;
	
	#if denominator is less than $threshold
	if(($parall_path2_3pos{$heter}) < $threshold){
		$dsp = "NA";
	}else{
		$dsp = $parall_path2_3pos{$heter} / ($parall_path2_3pos{$heter} + $parall_path2_3pos{$homo_1} + $parall_path2_3pos{$homo_2});
	}


	#dNp

	my $dnp;
	#my $num=$parall_path2_1pos{$heter}+$parall_path2_2pos{$heter};
	
	if(($parall_path2_1pos{$heter} + $parall_path2_2pos{$heter}) < $threshold){
		$dnp = "NA";
	}else{
		$dnp = ($parall_path2_1pos{$heter} + $parall_path2_2pos{$heter}) / ($parall_path2_1pos{$heter} + $parall_path2_1pos{$homo_1} + $parall_path2_1pos{$homo_2} + $parall_path2_2pos{$heter} + $parall_path2_2pos{$homo_1} + $parall_path2_2pos{$homo_2});
	}
	if($dsp eq "NA" || $dnp eq "NA"){
		$p_test = "NA";
	}
	else{
		$p_test = $dnp / $dsp;
	}
	
	if($p_test eq "NA"){
		print FILE_INFO_2 "NA     "; 
	}else{
		printf FILE_INFO_2 "\t%.4f     ", $p_test;
	}
	
	#RATIOS WITH PSEUDOCOUNTS:

	my $p_test_pseudocounts;

	#dSp

	my $dsp_pseudocounts = ($parall_path2_3pos{$heter} + 1) / ($parall_path2_3pos{$heter} + $parall_path2_3pos{$homo_1} + $parall_path2_3pos{$homo_2} + 1);

	#dNp

	my $dnp_pseudocounts = ($parall_path2_1pos{$heter} + $parall_path2_2pos{$heter} + 1) / ($parall_path2_1pos{$heter} + $parall_path2_1pos{$homo_1} + $parall_path2_1pos{$homo_2} + $parall_path2_2pos{$heter} + $parall_path2_2pos{$homo_1} + $parall_path2_2pos{$homo_2} + 1);
	
	$p_test_pseudocounts = $dnp_pseudocounts / $dsp_pseudocounts;
	
	
	printf FILE_INFO_3 "\t%.4f     ", $p_test_pseudocounts;
	


}


printf FILE_INFO "\n";
printf FILE_INFO_2 "\n";
printf FILE_INFO_3 "\n";



#####################

#2) output about nonparallel substitutions (for estimation of trivial dn/ds analog)




my $filename4 =
"dn_ds_analog_raw_counts_" 
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";
  
my $filename5 =
"dn_ds_analog_with_threshold_"
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";

 my $filename6 =
 "dn_ds_analog_with_pseudocounts_"
  . $sp_name_1 . "_"
  . $sp_name_2 . "_"
  . $sp_name_3 . "_"
  . $sp_name_4 . ".csv";


my $file4 = ">>" . $output_dir . $filename4;
open( FILE_INFO_4, $file4 ) or die("Could not open file");

my $file5 = ">>" . $output_dir . $filename5;
open( FILE_INFO_5, $file5 ) or die("Could not open file");

my $file6 = ">>" . $output_dir . $filename6;
open( FILE_INFO_6, $file6 ) or die("Could not open file");




print FILE_INFO_4 "sp_name_1    sp_name_2    sp_name_3     sp_name_4     LCA_dist     AC_1pos    AA_1pos    CC_1pos    AC_2pos    AA_2pos    CC_2pos    AC_3pos    AA_3pos    CC_3pos     AG_1pos    AA_1pos    GG_1pos    AG_2pos    AA_2pos    GG_2pos    AG_3pos    AA_3pos    GG_3pos     AT_1pos    AA_1pos    TT_1pos    AT_2pos    AA_2pos    TT_2pos    AT_3pos    AA_3pos    TT_3pos     CG_1pos    CC_1pos    GG_1pos    CG_2pos    CC_2pos    GG_2pos    CG_3pos    CC_3pos    GG_3pos     CT_1pos    CC_1pos    TT_1pos    CT_2pos    CC_2pos    TT_2pos    CT_3pos    CC_3pos    TT_3pos     TG_1pos    TT_1pos    GG_1pos    TG_2pos    TT_2pos    GG_2pos    TG_3pos    TT_3pos    GG_3pos\n";
print FILE_INFO_4 "$sp_name_1    $sp_name_2    $sp_name_3     $sp_name_4     $lca_dist     ";


print FILE_INFO_5 "sp_name_1     sp_name_2     sp_name_3     $sp_name_4     LCA_dist     AC     AG     AT     CG     CT     TG\n";
print FILE_INFO_5 "$sp_name_1     $sp_name_2     $sp_name_3     $sp_name_4     $lca_dist     ";


print FILE_INFO_6 "sp_name_1     sp_name_2     sp_name_3     $sp_name_4     LCA_dist     AC     AG     AT     CG     CT     TG\n";
print FILE_INFO_6 "$sp_name_1     $sp_name_2     $sp_name_3     $sp_name_4     $lca_dist     ";





foreach my $val (sort keys %parall_path1_1pos) {
	
	(my $n1, my $n2) = (split("",$val))[0,1];
	my $double_1=$n1 . $n1;
	my $double_2=$n2 . $n2;
	

	print FILE_INFO_4 "$overall_1pos{$val}     $overall_1pos{$n1 . $n1}     $overall_1pos{$n2 . $n2}     ";
	print FILE_INFO_4 "$overall_2pos{$val}     $overall_2pos{$n1 . $n1}     $overall_2pos{$n2 . $n2}     ";
	print FILE_INFO_4 "$overall_3pos{$val}     $overall_3pos{$n1 . $n1}     $overall_3pos{$n2 . $n2}     ";
	


	#RATIOS WITH THRESHOLD:

	my $threshold = 2;
	my $dnds_analog;

	my $ds;
	
	if(($overall_3pos{$val}) < $threshold){
		$ds = "NA";
	}else{
		$ds = $overall_3pos{$val}/($overall_3pos{$val} + $overall_3pos{$n1 . $n1} + $overall_3pos{$n2 . $n2});	
	}
	
	my $dn;

	if(($overall_1pos{$val} + $overall_2pos{$val}) < $threshold){
		$dn = "NA";
	}else{
		$dn = ($overall_1pos{$val}+$overall_2pos{$val})/($overall_1pos{$val} + $overall_1pos{$n1 . $n1} + $overall_1pos{$n2 . $n2} + $overall_2pos{$val} + $overall_2pos{$n1 . $n1} + $overall_2pos{$n2 . $n2});
	}
	if($ds eq "NA" || $dn eq "NA"){
		$dnds_analog = "NA";
	}else{
		$dnds_analog = $dn / $ds;
	}


	if($dnds_analog eq "NA"){
		print FILE_INFO_5 "NA     "; 
	}else{
		printf FILE_INFO_5 "\t%.4f     ", $dnds_analog;
	}

	
	

	#RATIOS WITH PSEUDOCOUNTS:

	
	my $dnds_analog_pseudocounts;

	my $ds_pseudocounts = ($overall_3pos{$val} + 1)/($overall_3pos{$val} + $overall_3pos{$n1 . $n1} + $overall_3pos{$n2 . $n2} + 1);	
	my $dn_pseudocounts = ($overall_1pos{$val} + $overall_2pos{$val}+1)/($overall_1pos{$val} + $overall_1pos{$n1 . $n1} + $overall_1pos{$n2 . $n2} + $overall_2pos{$val} + $overall_2pos{$n1 . $n1} + $overall_2pos{$n2 . $n2} + 1);
	
	$dnds_analog_pseudocounts = $dn_pseudocounts / $ds_pseudocounts;

	printf FILE_INFO_6 "\t%.4f     ", $dnds_analog_pseudocounts;
	
		


}

printf FILE_INFO_4 "\n";
printf FILE_INFO_5 "\n";
printf FILE_INFO_6 "\n";


close FILE_INFO;
close FILE_INFO_2;
close FILE_INFO_3;
close FILE_INFO_4;
close FILE_INFO_5;
close FILE_INFO_6;
closedir(RAWDIR);

