#!/usr/bin/perl -w

use Getopt::Long;
use Pod::Usage;
use FindBin;
use Data::Dump;
use File::Basename;
use Bio::Seq;
use Bio::SeqIO;

=head1 Function

        miRDeep.pl
        This programming for find species miRNA and predict their target gene
        This programming wish to analyze deep sequencing data mapping to target species genome for known and novel miRNA genes.

=head1 Date_Version_Author

        version 1.0

=head1 Usage and Parameter description

        $0 -reads reads.fa -genome genome.fa [options]
        --species    species name abbreviation
        --sp         species name
        --genome     a fasta file with the reference genome [needed]
        --index      build query sequence index
        --mature     a fasta file with mature sequence of this species [needed]
        --precursor  a fasta file with precursor of this species [needed]
        --mature_o   a fasta file with precurosr of other species [optional]
        --reads      a fasta file with the deep sequencing reads
        --length_filter  filter sequence short this length [Default:18]
        --help       help information

=head1 Example
perl step1_find_miRNA.pl -genome cel_cluster.fa -index  cel -reads  reads.fa -adapter  TCGTATGCCGTCTTCTGCTTGT -sp special_sp -precursor precursors_ref_this_species.fa -mature mature_ref_this_species.fa -species cel -mo mature_ref_other_species.fa -step 1234

=cut

GetOptions (
"help|h" => \$help,
"species=s" => \$species,
"sp" => \$special_sp,
"genome|g=s"  => \$genome,
"index=s" => \$index,
"mature|m=s" => \$mature,
"precursor|p=s" =>\$precursor,
"mature_o|mo:s" =>\$mature_others,
"reads|r=s" => \$reads,
"adapter:s" => \$adapter,
"length_filter:i" => \$length_filter,
"step=s" => \$step,
);
print "$step";

die `pod2text $0` if $help;
###  default parameter
$length_filter ||= 18;
$step ||= "1234" ;
my $index ||= "index";
###
#================================================================
#
#        step1 : build index
#
#================================================================
# usage : bowtie-build <ref> <ebwt_out>
# Default <ref> is fasta

if ($step =~ /1/){

`bowtie-build $genome $index`;

}
#===============================================================
#
#       step2 : process reads and map them to the genome
#
#===============================================================
#The -c option designates that the input file is a fasta file (for other input formats, see the README file). The -j options removes entries withnon-canonical letters (letters other than a,c,g,t,u,n,A,C,G,T,U,N). The -k option clips adapters. The -l option discards reads shorter than 18 nts.
#The -m option collapses the reads. The -p option maps the processed reads against the previously indexed genome (cel_cluster). The -s optiondesignates the name of the output file of processed reads and the -t option designates the name of the output file of the genome mappings. Last,-v gives verbose output to the screen.

#mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p cel_cluster -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v
if ($step =~ /2/){

        if (defined $adapter){

`perl mapper.pl $reads -c -j -k  $adapter -l $length_filter -m -p $index -s collapsed_$reads -t collapsed_vs_genome.$species.arf -v`;

}else{

`perl mapper.pl $reads -c -j -l $length_filter -m -p $index -s collapsed_reads -t collapsed_vs_genome.arf -v`;

        }
}
#==============================================================
#
#      step3 : fast quantitation of reads mapping to known miRBase precursors
#
#==============================================================
if ($step =~/3/){
my $y = `date +%m_%d`;
`perl quantifier.pl -p  $precursor -m $mature -r  collapsed_reads`;

}

#=============================================================
#
#      step4: identification of known and novel miRNAs in the deep sequencing data:
#
#=============================================================d
if ($step =~/4/){

`perl miRDeep2.pl  collapsed_reads $genome  collapsed_vs_genome.arf  $mature  $mature_others  $precursor 2> report.log`;

}
#=============================================================
#
#      step5: browse the results
#
#============================================================
#open the results.html using an internet browser. Notice that cel-miR-37 is predicted twice, since both potential precursors excised from this locus
#can fold into hairpins. However, the annotated hairpin scores much higher than the non-annotated one (miRDeep2 score 6.1e+4 vs. -0.2).
