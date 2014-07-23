#!/usr/bin/env perl
##
##	/shr/home/dherman/bin/merge-bed.pl			DSH, 04/04/09
##
##------------------------------------------------------

use Getopt::Long;
use File::Basename;
use strict;
use warnings;
use List::Util qw(max min);
do '/home/genetics/working/hermands/bin/seqAnalTools.pl';

my $this_script = basename($0);
my $helptext  = <<ENDHELP;
USAGE:  $this_script --ref <reference genome fasta> -shift <#> <.bed file>

OPTIONS
--method <merge|split|merge_same> # use merge_same to merge adjacent entries with the same name and ignore adjacent others
HISTORY
10-18-2010 	-skipped from no version to v3 (in order to skip existing v2)
			-added a no 'chr' --nochr flag
			-minimum chr position of 1 and maximum determined by .dict file that corresponds to the specified reference file

TODO
	- not properly handling same name entries that are non-adjacent in the 'merge_same' option

ENDHELP


my ($debug,$help,$shift_,$ref) = 0;
my $chr_=1;
my $SHIFT = 0;
my $method = 'merge';
&GetOptions(
    "help"  => 	\$help,
    "debug" =>	\$debug,
    "shift=i" => \$SHIFT,   #expand bed file
    "method=s" => \$method,   #default is 'merge'
    "chr!" => \$chr_,
    "ref=s" => \$ref,
);

die "Must specify reference!!!\n" unless ($ref);
print STDERR "Shifting .bed file by $SHIFT bases\n";
print STDERR "Reference = $ref\n";
my $DICT = $ref; $DICT =~ s/\.fasta/.dict~/;
#------------------------------------------------------
if ($help) {
	print STDERR $helptext;
	exit 1;
}

#------------------------------------------------------
#load AMPs
sub loadBed() {
    my ($bed_to_load,$frag_href,$dict) = @_;

    print STDERR "Loading in amplimers from $$bed_to_load\n";
    my $BEDin=();
    my $n=0;
    open $BEDin, $$bed_to_load or die "Cannot load in bed file $$bed_to_load\n";
    while (<$BEDin>) {
		next if($chr_ && ! m/^chr/);
		chomp;
		my ($chr,$start,$end,$name,@rest) = split;
		if ($chr eq "chr100") { $chr = "chrX" }
			elsif ($chr eq "chr1000") {	$chr = "chrY" }

		push @{$frag_href->{$chr}{NAMES}}, $name;
		push @{$frag_href->{$chr}{STARTS}}, max(1,$start - $SHIFT);
		push @{$frag_href->{$chr}{ENDS}}, min($end + $SHIFT,$dict->{$chr});
		push @{$frag_href->{$chr}{REST}}, join(";",@rest);
		$n++;
    }
    print STDERR "$n amplimers loaded\n";
    close $BEDin;
}
#--------------------
sub printBed() {
	my ($frags,$fout) = @_;

	open FOUT, ">$fout" or die "Cannot open $fout to print out.....\n";
	foreach my $chr (sort {$a cmp $b} keys %$frags) {
		for(my $i=0; $i <= $#{$frags->{$chr}{STARTS}};$i++) {
			print FOUT join("\t",($chr, ${$frags->{$chr}{STARTS}}[$i], ${$frags->{$chr}{ENDS}}[$i], ${$frags->{$chr}{NAMES}}[$i], split(/;/,${$frags->{$chr}{REST}}[$i]))), "\n";
		}
	}
	close FOUT;
}
#------------
my %frags =();
my $suffix = $method eq 'split' ? "${SHIFT}s" :
											$method eq 'merge' ? "${SHIFT}m" : "${SHIFT}ms";
my $bedIN = $ARGV[0];
my $bedOUT = $bedIN; $bedOUT =~ s/\.bed$/_${suffix}.bed/;

my %dictionary = &loadDict($DICT);
&loadBed(\$bedIN,\%frags,\%dictionary);

foreach my $chr (keys %frags) {
	my @order =(sort {${$frags{$chr}{STARTS}}[$a] <=> ${$frags{$chr}{STARTS}}[$b]} 0..$#{$frags{$chr}{STARTS}});
	my @starts = @{$frags{$chr}{STARTS}}[@order];
	my @ends = @{$frags{$chr}{ENDS}}[@order];
	my @names = @{$frags{$chr}{NAMES}}[@order];
	my @rest = @{$frags{$chr}{REST}}[@order];

	for(my $i=0;$i < $#starts;$i++) {
		print STDOUT "Processing $chr $names[$i] $starts[$i] $ends[$i] <=> $starts[$i+1] $ends[$i+1]\n" if ($debug);
		if ($ends[$i] >= $starts[$i+1]) {
			if ($method eq 'merge') { #merge
				print STDOUT "(1) STARTS: i=$i ", join(",",@starts[($i-1)..($i+2)]),"\n(1) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(1) Names: ",join(",", @names[($i-1)..($i+2)]),"\n" if ($debug);
				splice(@starts,$i+1,1);
				splice(@ends,$i,1);
				splice(@rest,$i+1,1);

				if ($names[$i+1] eq $names[$i]) {
					splice(@names,$i,1);
				} else {
					$names[$i+1] =~ s/^[^_]*_//;
	#Get rid of extra characters
					my ($n0,$n1,$n2,$n3) = @names[$i-1..$i+2];
					if ($i > 0) {
						$n0 =~ s/^.*[_-]//; $n0 =~ s/(?<=\d)\D.*$//;
					 } else { $n0 = (); }
					$n2 =~ s/^.*[_-]//; $n2 =~ s/(?<=\d)\D.*$//;
					$n1 =~ s/^.*_//; $n1 =~ s/-.*$//; $n1 =~ s/(?<=\d)\D.*$//;
					if ($i < $#starts - 1) {
						$n3 =~ s/^.*_//; $n3 =~ s/-.*$//; $n3 =~ s/(?<=\d)\D.*$//;
					} else { $n3 = ();}
					print STDOUT "Comparing $n0 cmp $n1  -- $n2 cmp $n3\n" if ($debug);
					if (($i > 0 && $n0 != $n1) || $i == 0) {
						if ($n1 =~ m/-/) {
							$names[$i] =~ s/(?<=\d)[^\d_][^_]*(?=-)//;
						} else {
							$names[$i] =~ s/(?<=\d)[^\d_][^_]*$//;
						}
					}
					if (($i < $#starts -1 && $n2 != $n3) || $i == $#starts - 1) {
						$names[$i+1] =~ s/(?<=\d)[^\d_-][^_-]*$//;
					}

					my $fname = join("-",@names[$i..($i+1)]);
					$fname =~ s	/-.*-/-/g;
					$fname =~ s/(?<=_)(.*)-\1$/$1/;
					splice(@names,$i,2,$fname);
				}
				print STDOUT "(2) STARTS: i=$i ", join(",",@starts[($i-1)..($i+2)]),"\n(2) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(2) Names: ",join(",", @names[($i-1)..($i+2)]),"\n\n" if ($debug);
				$i--;
			} elsif ($method eq 'merge_same') {
					print STDOUT "(1) STARTS: ", join(",",@starts[($i-1)..($i+2)]),"\n(1) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(1) Names: ",join(",", @names[($i-1)..($i+2)]),"\n" if ($debug);
					if ($names[$i+1] eq $names[$i]) { #merge if same name
						splice(@names,$i,1);
						splice(@starts,$i+1,1);
						splice(@ends,$i,1);
						splice(@rest,$i+1,1);
						$i--;
					}
					print STDOUT "(2) STARTS: ", join(",",@starts[($i-1)..($i+2)]),"\n(2) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(2) Names: ",join(",", @names[($i-1)..($i+2)]),"\n\n" if ($debug);
			} elsif ($method eq 'split') {
					print STDOUT "(1) STARTS: ", join(",",@starts[($i-1)..($i+2)]),"\n(1) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(1) Names: ",join(",", @names[($i-1)..($i+2)]),"\n" if ($debug);
					my $mid = sprintf("%d",($starts[$i+1] + $ends[$i])/2);
					$starts[$i+1] = $mid + 1;
					$ends[$i] = $mid;
					print STDOUT "(2) STARTS: ", join(",",@starts[($i-1)..($i+2)]),"\n(2) Ends: ",join(",", @ends[($i-1)..($i+2)]),"\n(2) Names: ",join(",", @names[($i-1)..($i+2)]),"\n\n" if ($debug);
			}
		}
	}
	@{$frags{$chr}{STARTS}} = @starts;
	@{$frags{$chr}{ENDS}} = @ends;
	@{$frags{$chr}{NAMES}} = @names;
	@{$frags{$chr}{REST}} = @rest;
	}
#---------------------------

&printBed(\%frags,$bedOUT);


__END__
