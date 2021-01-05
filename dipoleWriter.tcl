#Name:
#dipoleWriter
#Synopsis:
#	A Tcl scripts to write total dipole moments of selection (dipoleWriterSel) or position multiplied by charge (vqWriterSel) to ASCI file.
#Parameters:
#	selection - VMD selection for which the dipole will be outputed
#	sFileName - name of file to write down the dipole moments
#	charges - list of charges
#	start - frame index to start
#	stop - frame index to stop
#	step - step 
#Example:
#   set sel [atomselect top protein]
#	set charge [$sel get charge]
# 	dipoleWriterSel $sel protein_part.dip
#   vqWriterSel  $sel protein_part.vq $charge 0 999 2
##################################################################################################################
source getVolume.tcl 
puts "Loaded dipoleWriter.tcl"
puts "avail:"
puts "dipoleWriterSel { selection sFileName {start 0} {stop -1} {step 1}}"
puts "vqWriterSel {selection sFileName charges {start 0} {stop -1} {step 1} }"

#Write the dipole moment of given $selection from frame $start to $stop with $step to file $sFileName 
proc dipoleWriterSel { selection sFileName {start 0} {stop -1} {step 1} } {
set molID [$selection molid]
set fp [open $sFileName w]

puts $fp "\# BoxVolume: [getVol $molID 0]        A^3"
puts $fp "frame    dip_x     dip_y      dip_z    |dip|"
if {$stop==-1} {
	set nf [molinfo $molID get numframes]
} else {
	set nf $stop
}
for {set i $start} {$i <= $nf} {incr i $step} {
	$selection frame $i
    $selection update
    if {! [catch {measure dipole $selection -debye -geocenter} vector]} {
	    puts $fp "$i\t\t[lindex $vector 0] [lindex $vector 1] [lindex $vector 2] [veclength $vector]"
	}
}
close $fp
puts "Dipole moment of [$selection text] in molecule $molID from frame $start to frame $nf written to $sFileName."
}

proc vqWriterSel {selection sFileName charges {start 0} {stop -1} {step 1} } {
	set n [$selection num]
	set fo [open $sFileName w]
	set molID [$selection molid]
	puts $fp "\# BoxVolume: [getVol $molID 0]        A^3"
	puts -nonewline $fp "frame\t\t"
	foreach atom_idx [$selection get index] {
	   puts -nonewline "$atom_idx_xq     $atom_idx_yq      $atom_idx_z\t\t"
	}
	puts -nonewline $fo "\n"
	if {$stop==-1} {
		set nf [molinfo $molID get numframes]
	} else {
		set nf $stop
	}
	
	for {set i $start} {$i <= $nf} {incr i $step} {
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
