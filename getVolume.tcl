#Name:
#getVolume
#Synopsis:
#	A wrapper Tcl script to get system volume utilizing PBCTools Plugin of VMD.
#Note:
#	The output value is volume in A^3.
#Parameters:
#	molID - molecule id
#	frame_idx - index of frame 
#Output:
#	Volume of box of system with id "molID" at frame "frame_idx".
#Example:
# 	set volume [getVol 0 now]
##################################################################################################################
puts "Loaded getVolume.tcl"
puts "avail:"
puts "getVol {molID}"

proc getVol {molID frame_idx} {
   set dimensions [lindex [pbc get -molid $molID -first $frame_idx -last $frame_idx] 0]
   set alpha [cosine [lindex $dimensions 3]]
   set beta [cosine [lindex $dimensions 4]]
   set gama [cosine [lindex $dimensions 5]]
   return [expr [lindex $dimensions 0]*[lindex $dimensions 1]*[lindex $dimensions 2]*[unsqrt $alpha $beta $gama]]
   }
   
proc deg2rad {deg} {
  return [expr $deg*4.0*atan(1.0)/180]
}

proc rad2deg {rad} {
  return [expr 180*$rad/(4.0*atan(1.0))]
}

proc cosine {degree} {
	expr {cos([deg2rad $degree])}
}

proc unsqrt {alpha beta gama} {
	expr {sqrt(1+2*$alpha*$beta*$gama-$alpha*$alpha-$beta*$beta-$gama*$gama)}
}
