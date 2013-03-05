#!/usr/bin/perl

open(FILE, "cpu.txt");
open(OUTFILE, ">cpu_idl.txt");

$count = 0;

while($line=<FILE>)
{
    $count++;

    if($count == 25)
    {
	chop $line;
	($first, $second)= split(":", $line);
	($time, $second)= split(",", $second);
	print OUTFILE "$time ";
    }

    for($i=0; $i<33; $i++)
    {
        $line=<FILE>;

	if($count == 25)
	{
	    @fields = split ' ' , $line;
	    print OUTFILE "$fields[1] ";
	}
    }

    $line=<FILE>;   # empty line

    if($count == 25)
    {
	print OUTFILE "\n";
	$count = 0;
    }
	
}


