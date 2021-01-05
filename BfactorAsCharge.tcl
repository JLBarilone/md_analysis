# tcl skript for assingning charge in VMD from B-factor value in pdb
# USAGE: just load the pdb in VMD as top and source this script
puts "Loaded BfactorAsCharge.tcl"
puts "avail: assignChargeFromBeta {}"
proc assignChargeFromBeta {molID} {
	set ev [atomselect $molID all]
    set indices [$ev get index]
   
    foreach ind $indices {
	    set at [atomselect $molID "index $ind"]
        $at set charge [$at get beta]
		$at delete
    } 
$ev delete
}
