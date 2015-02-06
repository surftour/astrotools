#!/usr/local/bin/perl

open(FILE, "cpu.txt");
open(OUTFILE, ">cpu_bench.txt");


for($step = 0; $step < 3; $step++)
{
    $line=<FILE>;   # Step number
    for($i=0; $i<27; $i++)
    {
        $line=<FILE>;
        @fields = split ' ' , $line;

        if($step == 0)
        {
            $names[$i] = $fields[0];
            $cpu0[$i] = $fields[1];
        }
        
        if($step == 2)
        {
            $cpu2[$i] = $fields[1];
        }
    }

    $line=<FILE>;   # empty line
}

for($i=0; $i<27; $i++)
{
    printf OUTFILE "%15s     %8.3f  %7.1f\%\n", $names[$i], 0.5 * ($cpu2[$i] - $cpu0[$i]), 
                             100.0*($cpu2[$i] - $cpu0[$i])/($cpu2[0] - $cpu0[0]);
}



