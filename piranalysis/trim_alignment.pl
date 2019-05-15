#!/usr/local/bin/perl -w
#
# trim_alignment.pl
#
# Identify blocks of relatively ungapped sequence in an alignment
#
# AUTHOR: Adam Reid
# Copyright (C) 2019 Genome Research Ltd.
# This program is distributed under the terms of the GNU General Public License

use strict;

my($file) = @ARGV;

my $min_ungapped = 0.7;
my $min_within_block_bad = 3;
my $min_block_size = 10;
my $link_max = 1;

my $gap_char = '-';

open(F, "<$file") or die "$!";
my @seqs = <F>;
my $seqs = parse_fasta_id(\@seqs);
undef @seqs;
close F;

my @ungapped = ();

my $total_seqs = scalar keys %$seqs;

my %starts = ();
my %ends = ();

foreach my $id (sort keys %$seqs)
{
	my @seq = split //, $seqs->{$id};

	my $c = 0;

	my $current_end = '';

	foreach my $pos (@seq)
	{
		$ungapped[$c]++ unless $pos eq $gap_char;

		if($pos ne $gap_char &&!exists $starts{$id})
		{
			$starts{$id} = $c;
		}
		
		if($pos ne $gap_char)
		{
			$current_end = $c - 1;
		}

		$c++;
	}

	$ends{$id} = $current_end;
}

my @best_superblock = ();

my $window_size = $min_block_size;

my $superblock_start = 0;
my $link_count = 0;
my $in_block = 0;
my $best_end = 0;

for(my $i=0;$i<scalar (@ungapped);$i+=$window_size)
{
	my $good = 0;

	for(my $j=0;$j<$window_size;$j++)
	{
		my $pos = $i+$j;

		if ($pos >= scalar @ungapped)
		{
			# Ignore end of alignment if less than window size
			$good = 0;
			last;
		}

		$good++ if ($ungapped[$pos]/$total_seqs) >= $min_ungapped;
	}
	
	if(($good + $min_within_block_bad) >= $window_size)
	{
		print STDERR "$i\t$good\t$window_size\tGood\n";
		
		if(!$in_block)
		{
			$superblock_start = $i;
			$in_block = 1;
		}

		$best_end = $i+$window_size;
	}
	else
	{
		if($in_block)
		{
			$link_count++;

			if($link_count > $link_max)
			{
				# Break superblock
				print STDERR "SUPERBLOCK $superblock_start $best_end\n";

				if (!defined $best_superblock[0] || (($best_end + 1 - $superblock_start) > ($best_superblock[1] +1 - $best_superblock[0])))
				{
					@best_superblock = ($superblock_start, $best_end);
				}
		
				$in_block = 0;
				$superblock_start = 0;
				$link_count = 0;
				$best_end = 0;
			}
		}
	}
}

if($in_block)
{
	print STDERR "SUPERBLOCK $superblock_start $best_end\n";

	if(!defined $best_superblock[0] || (($best_end + 1 - $superblock_start) > ($best_superblock[1] +1 - $best_superblock[0])))
	{
		@best_superblock = ($superblock_start, $best_end);
	}
}

if(!defined $best_superblock[0])
{
	print STDERR "No blocks found in $file!\n";
	exit;
}

print STDERR "Best SUPERBLOCK: $best_superblock[0] to $best_superblock[1]\n";

foreach my $id (sort keys %$seqs)
{
	if ($starts{$id} <= $best_superblock[0] && $ends{$id} >= $best_superblock[1])
	{
		my $new_seq = substr($seqs->{$id}, $best_superblock[0], ($best_superblock[1] + 1 - $best_superblock[0]));

		print ">$id\n$new_seq\n";
	}
	else
	{
		print STDERR "Removing: $id\t$starts{$id}\t$ends{$id}\n";
	}
}

sub parse_fasta_id
{
        $/ = "\n";

        my $input = shift;
        my %seqs;
        my $header;
        my $seq;

        foreach my $line (@$input)
        {
                chomp $line;

                if ($line =~ /^>(\S+)/ && defined $seq)
                {
                        $seqs{$header} = $seq;
                        $header = $1;
                	$seq = '';
                }
       	        elsif ($line =~ /^>(\S+)/)
                {
	                $header = $1;
                }
                else
                {
                        $seq .= $line;
                }
        }

        $seqs{$header} = $seq;

        return \%seqs;
}

