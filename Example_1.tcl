# This example script will read three trajetories to VMD and write down dipole moment time evolution of protein to ASCII file.

source BfactorAsCharge.tcl
source dipoleWriter.tcl
	
for {set n 1} {$n <= 3} {incr n} {
	# load trajectory
	set molID [mol new ../trajectories/topo.pdb type {pdb} first 0 last -1 step 1 waitfor -1]
	mol addfile ../trajectories/md_1ns_${n}.xtc type {xtc} first 0 last -1 step 1 waitfor -1
	# remove first frame that corresponds to topo.pdb
	animate delete  beg 0 end 0 skip 0 $molID
	# assign charges
	assignChargeFromBeta $molID
	# make selection
	set sel [atomselect $molID "not water"]
	# write dipole moment evolution of protein to file
	dipoleWriterSel $sel md_1ns_${n}.dip
	# delete molecule
	mol delete $molID
}
