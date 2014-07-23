#!/usr/bin/env perl
use warnings;
do "/home/dsh7/bin/print_tools.pl";
#--------------
sub frag2reg() {
	my ($frag_href,$chr,$i) = @_;
	return $chr . ":" . ${$frag_href->{$chr}{STARTS}}[$i] . "-" . ${$frag_href->{$chr}{ENDS}}[$i];
}
#----------------------------
sub timeStamp() { 
	my ($sec,$min,$hour,$mday,$mon,$year)=localtime(time);  
	return(sprintf"%s/%s/%s %02d:%02d:%02d",$mon,$mday,1900+$year,$hour,$min,$sec) 
}
#----------------------
sub median() {
	my ($a) = @_;
	my $l = scalar @{$a};
	if ($l % 2) {
		my $m = sprintf("%d", $l / 2);
		return ${$a}[$m];
	} else {
		return sprintf("%.1f",(${$a}[($l/2-1)] + ${$a}[($l/2)])/2);
	}
}
#--------------
sub mean() {
	my ($a) =@_;
	my $total =0;
 	foreach my $el (@{$a}) { 	$total += $el }
	return ($total / (scalar @{$a}));
}
#---------------
sub fillinvector() {
#define all undefined values (0..$n_vals) as $val
	my ($vec,$val,$n_vals) = @_;
	$n_vals = $#{$vec} + 1 if (! defined($n_vals));
	for(my $i=0;$i < $n_vals; $i++) {
		$vec->[$i]=$val if (! defined($vec->[$i]));
	}
}
#---------------------------
sub fillinmatrix() {
#define all undefined values in 0..($n_cols - 1) columns and 0..$#{$mat} rows in matrix as 0
	my ($mat, $val,$n_cols) = @_;
	for(my $i =0;$i <=$#{$mat};$i++) {
		&fillinvector(\@{$mat->[$i]},$val,$n_cols);
	}
}
#-------------------------
sub findMedian() {

	my ($vec,$n_values) = @_;
	$n_values = sum(@$vec) if (! defined($n_values));
	return 0 if (! $n_values);
	
	my $is_odd = ($n_values % 2);
	my $midpoint = sprintf("%d",$n_values / 2) + 1;
	
	my $cur_count = 0;
	my $i=0;
	while($cur_count < $midpoint) {
		$cur_count += $vec->[$i];
		$i++;
	}
	if (! $is_odd && $vec->[$i-1] == 1) {
		return $i - 3/2;
	}
	return $i -1;
}
#---------------------------------
sub numBasesBelow() {
#sum of base counts below limit
	my ($vec,$limit) = @_;
	
	my $cur_count = 0;
	my $max = min(scalar @$vec,$limit);
	for(my $i=0;$i < $max;$i++) {
		$cur_count += $vec->[$i];
	}
	return $cur_count;
}
#---------------------------------
sub numBasesAbove() {
#sum of base counts above limit
	my ($vec,$limit) = @_;
	
	my $cur_count = 0;
	my $min = max(0,$limit);
	for(my $i=$#{$vec};$i > $min;$i--) {
		$cur_count += $vec->[$i];
	}
	return $cur_count;
}
#------------------------------
sub numBasesOutside() {
#sum of base counts between a & b (inclusive)
# REQUIRE: a < b
	my ($vec,$a,$b) = @_;
	return &numBasesBelow($vec,$a) + &numBasesAbove($vec,$b);
}
#------------------------
sub rPos() {
	my ($pos,$start,$end) = @_;
	#1-based position
	my $pos_l = $pos - $start + 1;
	my $pos_r = $end - $pos + 1;
	return $pos_l <= $pos_r ? $pos_l : -1 * $pos_r;
}
#---------------------------------
sub uniquifyArray() {
	my %seen = ();
	grep { ! $seen{$_} ++ } @_;
}
#===================
sub chr_sort {
#Only handle #, X, Y
	my ($A, $B) = ($a,$b);
	$A =~ s/^chr//;	$B =~ s/^chr//;
	my $a_numeric = ($A =~ m/^\d+$/);
	my $b_numeric = ($B =~ m/^\d+$/);
	
	if ($a_numeric && $b_numeric) {
		return $A <=> $B;
	} elsif (! $a_numeric && ! $b_numeric) {
#		return 0 if ($a eq $b);
		return -1 if ($A eq "X" && $B eq "Y");
		return 1 if ($A eq "Y" && $B eq "X");
	} elsif ($a_numeric) {
		return -1;
	} else {
		return 1;
	}
}
#=======================================================-=-=-=-=-=
sub refQuery() {
	#handle whether 'chr' in reference seek and shift from 0-based non-inclusive to 1-based inclusive coordinates
	my ($fai,$chr,$start,$end) = @_;

	$chr =~ s/^chr//	if (! $chr_keep);
	return($fai->fetch($chr . ':' . ($start + 1) . '-' . $end));
}
#---------------------
sub connectDB() {

	my ($DB,$host,$user) = @_;

	#==============
	#Get password
	open PF, "/home/$user/.my.cnf" or die "Cannot open .my.cnf\n";
	my @l = <PF>;
	close PF;
	my $l = join("",@l);
	$l =~ s/ //g;
	$l =~ m/(?<=password=)(\S*)(?=\n)/;
	my $password = $1;
	if (! $password) {
		&printSTDERR("Enter MySQL user password for USER ${user} DB $DB @ $host: \n"); 
		$password = <STDIN>;chomp $password;
	}
	#====================
	my $genome = DBI->connect("DBI:mysql:$DB:$host",$user,$password);
	if (! $genome) {
		&printSTDERR("Enter MySQL user password for USER ${user} DB $DB @ $host: \n"); 
		$password = <STDIN>;chomp $password;
		$genome = DBI->connect("DBI:mysql:$DB:$host",$user,$password) or die "Cannot connect to $DB @ $host\n";
	}
	print STDERR "Connected to $DB at $host with user $user\n";
	return $genome;
}
#-------------------
sub ntComp() {
	my $seq = $_[0];
    $seq =~ tr/ACGTacgt/TGCATGCA/;
    return $seq;
}
#-------------------
sub RC() {
    my $seq = $_[0];
	$seq =~ tr/ACGTacgt/TGCATGCA/;
    $seq = reverse $seq if ($seq =~ m/^[ACGT]+$/);
    return $seq;
}
#=================
sub loadDict() {

	my ($dict_file) = @_;
	my %dict=();
	open DF, "$dict_file" or die "Cannot open dictionary $dict_file\n";
	while (<DF>) {
		my ($name,$length) = (m/SN:(\S*)\s+LN:(\S*)/);
		$dict{$name} = $length;	
	}
	close DF;
	return %dict;
}

1;