#Name:
#alignTraj
#Synopsis:
#	A Tcl scripts for reduction of rotational degrees of freedom by align the molecules. 
#	Iterates thru all frames available for slected ID molecule. And transform the position of all atoms in the frame by the transformation matrix resulting from alignement procedure of selection "centerMolSelText" between actual frame and "patternFrameNum" frame.
#Parameters:
#	mID - VMD molecule ID
#	centerMolSelText - VMD selection text
#	patternFrameNum - number of tamplate frame 
#	weights - possible weights for aligning
#Example:
#   centerByMol 0 "backbone" 1
##################################################################################################################
puts "avail:"
puts "centerByMol { mID centerMolSelText patternFrameNum {weights -1} }"

proc centerByMol { mID centerMolSelText patternFrameNum {weights -1} } {
set molID $mID
set sel1 [atomselect $molID $centerMolSelText]
set sel2 [atomselect $molID $centerMolSelText]
set sel3 [atomselect $molID "all"]
$sel1 frame $patternFrameNum
set nf [molinfo $molID get numframes]
for {set i 0} {$i < $nf} {incr i 1} {
	$sel2 frame $i
	if {$weights == -1} {
		set M [measure fit $sel2 $sel1]
	} else {
	set M [measure fit $sel2 $sel1 weight $weights]
	}
    $sel3 frame $i
	$sel3 move $M
}
}
