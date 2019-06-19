#!/usr/bin/perl -w

# cgmaptools - combineCGStrands.pl
#
# Copyright (C) Benjamin Berman
# Contact: Benjamin Berman <benbfly@gmail.com>
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;

=pod

=head1 DESCRIPTION

    Takes a single CGmap input file, and combines the two strands for all
    palindromic pairs (i.e. will combine the "C" and the adjacent "G" for
    a CG context). Outputs a CGmap file with the two lines combined into
    one, with counts summed.

    In principle, should work for CHG context as well, but this is not yet
    implemented. Only lines with context CG are merged

=head1 USAGE

    cgmaptools combinestrands [options] -i <input> -o <output>

    Options:
    -i   Input file. Name end with .CGmap or .CGmap.gz.
    -o   Output results to file [default: stdout]
    -h   Help message.

=head1 AUTHOR

    Contact:     Berman, Benjamin; benbfly@gmail.com
    Last update: 2019-05-30

=cut

## Parsing arguments from command line
my ( $infn, $outfn, $help);

GetOptions(
    'i=s' => \$infn,
    'o:s' => \$outfn,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 2, -output => \*STDERR } ) if ( $help );
pod2usage() unless ($infn and $outfn);


die "Can't write to output file ${outfn}\n" unless (open(OUTFH, ">$outfn"));
die "Can't read input file ${infn}\n" unless (open(INFH, "$infn"));


# State vars
my @prev_flds = map {0} (1..8);  # Use all 0s as an empty line

LINE: while (my $line = <INFH>)
{
    chomp $line;
    my @f = split("\t",$line);
    die "Following line does not have 8 fields: $line\n" unless (scalar(@f)==8);

    my ($chr, $nuc, $pos, $cont, $dinuc, $meth, $mc, $nc) = @f;
    my ($prev_chr, $prev_nuc, $prev_pos, $prev_cont, $prev_dinuc, $prev_meth, $prev_mc, $prev_nc) = @prev_flds;

    # Combine logic
    my $combine = 0;

    # Check CpG
    if (($cont eq "CG") &&
	($prev_cont eq "CG") &&
	($nuc eq "G") && 
	($prev_nuc eq "C") &&
	($prev_chr eq $chr) &&
	(($pos - $prev_pos) == 1))
    {
	$combine = 1;
    }
    
    # print STDERR join("\t",$nuc,$prev_nuc,$prev_chr, $chr, $pos, $prev_pos, $combine)."\n"; # DEBUG

    # If combine, we output the combined line.  If not, we output the previous line
    # (this could be different if CHG is implemented).
    my @outf = @prev_flds;
    if ($combine)
    {
	$outf[6] = $prev_mc + $mc;
	$outf[7] = $prev_nc + $nc;
	$outf[5] = sprintf("%0.2f",$outf[6]/$outf[7]);

	# This line does not become a previous line, because we've already taken care of it.
	# So zero it out;
	@f = map {0} (1..8);
    }
    
    print OUTFH join("\t",@outf)."\n" if ($outf[0]); # unless it's an empty line
    
    # Update state vars
    @prev_flds = @f;
}

# We must output the final one if it was not combined.
print OUTFH join("\t",@prev_flds)."\n" if (@prev_flds);

close(INFH);
close(OUTFH);
