#!/usr/bin/perl

use strict;
use Spreadsheet::ParseExcel;
my $oBook = Spreadsheet::ParseExcel::Workbook->Parse(shift);
my($iR, $iC, $oWkS, $oWkC);

$oWkS = $oBook->{Worksheet}[0];
for(my $iR = $oWkS->{MinRow} ;
    defined $oWkS->{MaxRow} && $iR <= $oWkS->{MaxRow} ; $iR++) {
    for(my $iC = $oWkS->{MinCol} ; defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ; $iC++) {
	$oWkC = $oWkS->{Cells}[$iR][$iC];
	print "," unless $iC == $oWkS->{MinCol};
	print $oWkC->Value if($oWkC);
    }
    print "\n";
}
