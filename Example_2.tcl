# This example script will read three trajetories to VMD, align the molecules along trajectory according to molecule in topo.pdb and save the trajecetories to new file.

source alignTraj.tcl
	
for {set n 1} {$n <= 3} {incr n} {
	# load trajectory
	set molID [mol new ../trajectories/topo.pdb type {pdb} first 0 last -1 step 1 waitfor -1]
	mol addfile ../trajectories/md_1ns_${n}.xtc type {xtc} first 0 last -1 step 1 waitfor -1
	# align the molecules by backbone
	centerByMol $molID "backbone" 0
	# save the trajectory (except the 0th frame)
	animate write dcd md_1ns_aligned_${n}.dcd beg 1 end -1 skip 1 top
	# delete molecule
	mol delete $molID
}
