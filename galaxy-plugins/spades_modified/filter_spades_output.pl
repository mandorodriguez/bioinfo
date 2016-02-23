#!/usr/bin/env perl -w

=head1 SYNOPSIS

filter_spades_output.pl - Filters contigs or scaffolds based on contig length and coverage.

=head1 USAGE

filter_spades_output.pl [-c|--coverage-cutoff] [-l|--length-cutoff] [-o|--filtered-out out.fasta] -t|--tab stats.tab seqs.fasta

=head1 INPUT

=head2 [-c|--coverage-cutoff]

Mininum coverage. Contigs with lower coverage will be discarded. Default 10.

=head2 [-l|--length-cutoff]

Mininum coverage. Smaller ontigs will be discarded. Default 500.

=head2 -t|--tab stats.tab

A tabular file, with three columns: contig name, length, and coverage:

NODE_1	31438	24.5116
NODE_2	31354	2316.96
NODE_3	26948	82.3294

Such a file is produced by spades.xml. Contigs should be in the same order as in the fasta file.

=head2 [-o|--filtered-out out.fasta]

If specified, filtered out sequences will be written to this file.

=head2 seqs.fasta

Sequences in fasta format. Start of IDs must match ids in the tabular file. 

=head1 OUTPUT

A fasta file on stdout.

=head1 AUTHOR

Lionel Guy (lionel.guy@icm.uu.se)

=head1 DATE

Thu Aug 29 13:51:13 CEST 2013

=cut

# libraries
use strict;
use Getopt::Long;
use Bio::SeqIO;

my $coverage_co = 10;
my $length_co = 500;
my $out_filtered;
my $tab_file;

GetOptions(
    'c|coverage-cutoff=s' => \$coverage_co,
    'l|length-cutoff=s' => \$length_co,
    'o|filtered-out=s' => \$out_filtered,
    't|tab=s' => \$tab_file,
);
my $fasta_file = shift(@ARGV);
die ("No tab file specified") unless ($tab_file);
die ("No fasta file specified") unless ($fasta_file);

## Read tab file, discard rows with comments
open TAB, '<', $tab_file or die "$?";
my @stats; 
while (<TAB>){
    chomp;
    push @stats, $_ unless (/^#/);
}

## Read fasta
my $seq_in = Bio::SeqIO->new(-file => $fasta_file, 
			     -format => 'fasta');
my $seq_out = Bio::SeqIO->new(-fh => \*STDOUT, 
			      -format => 'fasta');
my $seq_out_filt = Bio::SeqIO->new(-file => ">$out_filtered", 
				   -format => 'fasta') if ($out_filtered);
while (my $seq = $seq_in->next_seq){
    my $stat = shift @stats;
    die "Less rows in tab than sequences in seq file" unless $stat;
    my ($id_tab, $length, $coverage) = split(/\t+/, $stat);
    die "id, length or coverate not defined at $stat\n" 
	unless ($id_tab && $length && $coverage);
    my $id_seq = $seq->id;
    die "Unmatched ids $id_seq and $id_tab\n" unless ($id_seq =~ /^$id_tab/i);
    if ($length >= $length_co && $coverage >= $coverage_co){
	$seq_out->write_seq($seq);
    } elsif ($out_filtered){
	$seq_out_filt->write_seq($seq);
    } else {
	# do nothing
    }
}
die "More rows in tab than sequences in seq file" if (scalar(@stats) > 0);
exit 0;

