

proc vqWriter {selection charges fStart fStop oFName} {
	set n [$selection num]
	set fo [open $oFName w]
	for {set i $fStart} {$i <= $fStop} {incr i 1} {
		$selection frame $i
		$selection update 
		set pos [$selection get {x y z}]
		puts -nonewline $fo "$i\t\t"			
		foreach xyz $pos charge $charges {
 			set xq [vecscale $charge $xyz]
	 		puts -nonewline $fo "${xq}\t\t"			
  		}	
  		puts -nonewline $fo "\n"
	}
	close $fo
}
