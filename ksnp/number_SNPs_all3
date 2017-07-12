#!/usr/bin/perl
#v3.03

no warnings 'deprecated';
my %fors=();

# number_SNPs_all.pl SNPs_all 

my $all_snps_file=$ARGV[0];
my $out=$all_snps_file."_labelLoci";



# Print all SNPs 

open(ALL_SNPS, ">$out") || die "Can't open [$out] for writing $!\n";

my $count=-1;
open IN, "$all_snps_file" or die "Cannot open all $all_snps_file: $!\n";
my $primer="X";
my $previous_primer="";
while (my $line = <IN>){
    
    if ($line =~ /(\S*)\t(\S*)\t(\d+(?:\sF|\sR)?|x)\t(.*)/ ) { 
	
	#print "$1\t$2\t$3\t$4\t$5\n";
	$primer = $1;
	my $id = $4;
	chomp($id);
	my $snp_char = $2;
	my $position = $3;

	if ($previous_primer ne $primer && $previous_primer ne "") {
	    print ALL_SNPS "\n";
	    $count++;
	    foreach my $id0 (sort keys %fors) {
		foreach my $position0 (sort {$a cmp $b} keys %{$fors{$id0}}) {
		    printf ALL_SNPS "$count\t$previous_primer\t$fors{$id0}{$position0}\t$position0\t$id0\n";
		}
	    }
	    %fors=();
	    $previous_primer=$primer;
	} elsif ($previous_primer eq "") {
	    $previous_primer=$primer;
	}
	$fors{$id}{$position} = $snp_char;
    }
}

# print last one        
print ALL_SNPS "\n";
$count++;
foreach my $id (sort keys %fors) {
    foreach my $position (sort {$a cmp $b} keys %{$fors{$id}}) {
	printf ALL_SNPS "$count\t$previous_primer\t$fors{$id}{$position}\t$position\t$id\n";
    }
}

close ALL_SNPS or die "Cannot close $out: $!.\n";
$count++; 
`echo "Number_SNPs: $count" > COUNT_SNPs`;
print "Number_SNPs: $count\n";
close IN;
