#!/usr/bin/perl

use File::Find;

if(@ARGV[0])
{
    while(@ARGV[0])
    {
	$f= pop @ARGV;
	process($f);
    }
}
else
{
    find(\&wanted, '.');
}

sub wanted
{
    $f = $_; # current filename

    if(!(-d $f))  # not a directory
    {
	process($f);

    }
}

sub process($f)
{
        $g= $f . ".pc";
	
	print "processing: ",$f,"\n";
	`/home2/tcox/Tools/SnapSwapEndian/Swap/Swap $f $g`;
	`/home2/tcox/Tools/SnapSwapEndian/SwapDblSnap/Swap $g`;
	`mv $g $f`;
}














