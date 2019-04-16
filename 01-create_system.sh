#!/bin/bash

for conc in 0.5 1.0 2.0 3.0 4.0 5.0
do
	Folder=${conc}M
	if [ ! -d "${Folder}" ] ; then
		echo "Folder does not exists"
		mkdir ${Folder}
	else
		echo "Folder does exists"
		rm -r ${Folder}
		mkdir ${Folder}
	fi
	cd ${Folder}/
	
	# Create System
	vmd -dispdev text <<EOF
package require solvate

solvate -minmax {{0 0 0}  {48 48 24}}  -o bottom
solvate -minmax {{0 0 24} {48 48 72}}  -o middle
solvate -minmax {{0 0 72} {48 48 108}} -o top

mol delete all
resetpsf

package require autoionize

autoionize -psf middle.psf -pdb middle.pdb -sc ${conc} -between 3.0 -cation SOD -anion CLA -o middle-ionized

resetpsf
mol new middle-ionized.pdb

set W [atomselect top "water"]
set I [atomselect top "ions"]
\$W writepdb water.pdb
\$I writepdb ions.pdb

mol delete all

package require psfgen

topology ../common/toppar_water_ions.str

segment WTT {
	pdb top.pdb
}
coordpdb top.pdb WTT

segment WTM {
	pdb water.pdb
}
coordpdb water.pdb WTM

segment WTB {
	pdb bottom.pdb
}
coordpdb bottom.pdb WTB

segment ION {
	pdb ions.pdb
}
coordpdb ions.pdb ION

writepsf NaCl-${conc}M.psf
writepdb NaCl-${conc}M.pdb

mol delete all

mol new     NaCl-${conc}M.psf
mol addfile NaCl-${conc}M.pdb

set all [atomselect top "all"]
\$all set beta 0.0
\$all set occupancy 0.0

set Ion [atomselect top "ions"]
\$Ion set occupancy 1.0

\$all writepdb tagged_atoms.pdb
\$all set beta 0.0
\$all set occupancy 0.0

mol delete all
EOF
	# Clean up
	rm middle* ions* bottom* top* water*
	
	cd ../
done

