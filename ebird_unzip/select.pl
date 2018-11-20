use strict;
use 5.10.0;

my ($S, $F) = @ARGV;

open S, $S;
open F, $F;

chomp(my $sel = <S>);
my @sel = split/,\s*/, $sel;

chomp(my $head = <F>);
my @head = split/\t/, $head;

my @cols;

foreach my $x (0..$#sel){
	my $s = $sel[$x];
	foreach my $y (0..$#head){
		my $h = $head[$y];
		if ($s eq $h){
			die "More than one match for $s in header" if defined $cols[$x];
			$cols[$x] = $y
		}
	}
	die "No match for $s in header" unless defined $cols[$x];
}

say join "\t", @head[@cols];

while (<F>){
	my @ln = split /\t/;
	say join "\t", @ln[@cols];
}
