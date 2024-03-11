#!/bin/bash

#Author: Thibault Tubiana, PhD, 2020
#Version: 1.3.0

#activate acpype environment
#conda activate acpype

#---------  FILE SETUP  -----------
#PDB FILE, will be replace in the code later by $PDB
FILE="" #<<<<<<<<<<<<<<<<<<<<<<<<< PUT THE PDB NAME HERE (without the extension)
# LIGAND NAME, if you have a ligand, it will be parametrize with acpype and the ligand name will be replace by "LIG".
LIGNAME="" #<<<<<<<<<<<<<<<<<<<<<<  #PUT LIGAND NAME HERE, leave it blank if no ligand.

#---------  SIMU SETUP  -----------
BOXSIZE=1.5 #cubic simulation boxsiwe
BOXTYPE=cubic #Box type
NT=20 #Number of core.
WATER=tip3p #Water type
NUMBEROFREPLICAS=1 #Number of replica
FF=amber99sb-ildn #Force field
SIMULATIONTIME=200 #Simulation time in nanosec. Will be converted in fs and modified in the mdp file.


	#Setup simulationtime
if [ ! -z "$SIMULATIONTIME" ]
then
python_command=$(python <<EOF
import re
restep = re.compile("nsteps *= *(\d*)")
redt = re.compile("dt *= *(\d*.\d*)")
dt = 0
simulationtime = float($SIMULATIONTIME) *1000 #Time in ps
outputLines = []
with open("mdp/md_prod.mdp",'r') as f:
    mdp = f.readlines()
    #find first the timestep
    for line in mdp: 
        dtmatch = redt.match(line)
        if dtmatch:
            dt = float(dtmatch.group(1))
            break
    for line in mdp:
        stepmatch = restep.match(line)
        if stepmatch and float(dt) > 0:
            nsteps = int(simulationtime)/dt
            line = "nsteps            = {}        ; {} * {} = {} ps or {} ns\n".format(int(nsteps),dt,nsteps, dt*nsteps, simulationtime/1000)
        outputLines.append(line)
    with open("mdp/md_prod.mdp",'w') as f:
        for line in outputLines:
            f.write(line)
EOF
)
fi


if [ ! -z "$LIGNAME" ]
then
      #Create parameters directories
      mkdir param
      cp $FILE".pdb" param/
      cd param

      grep 'ATOM  ' $FILE".pdb" --color=none > receptor.pdb

      #Extract ligand and connect
      python_command=$(python <<EOF
ligand_atom = []
keepLine = []
with open("$FILE.pdb","r") as file:
    lines = file.readlines()
    for line in lines:
        if '$LIGNAME' in line[17:20]:
            line = line[:17]+"LIG"+line[20:]
            keepLine.append(line)
            ligand_atom.append(int(line[6:11]))
        elif "CONECT" in line[0:6]:
            idx = [int(x) for x in line.split()[1:]]
            if any(id in idx for id in ligand_atom):
                keepLine.append(line)
with open("ligand.pdb","w") as file:
    for line in keepLine:
        file.write(line)
EOF
)


    #Convert in mol2 while adding hydrogens.
    obabel -ipdb ligand.pdb -omol2 -h > ligand.mol2

    #use ACPYPE to create the ligand topology.
    #DISCLAMER! This is a "quick and dirty method", it has to be optimised with ACPYPE parameters of course and adapted to ligands
    #if you see strange MD behaviours.
    # You may also consider Automated Topology Builder (ATB) (webserver) Or LibParGen (webserver & standalone tools)
    acpype -i ligand.mol2
    mkdir ligand
    mv ligand* ligand/
    mkdir receptor
    mv receptor.pdb receptor/
    cd receptor

    #Preparing topology
    $GMX pdb2gmx -f receptor.pdb -o receptor_GMX.pdb -water $WATER -ignh -ff $FF

    #Copy files from receptors/ligand folders.
    cd ../../
    cp param/receptor/*.itp param/receptor/topol.top .
    cp param/ligand/ligand.acpype/ligand_GMX.itp ligand.itp
    grep -h ATOM param/receptor/receptor_GMX.pdb param/ligand/ligand.acpype/ligand_NEW.pdb > complex.pdb

    #Add the ligand topology inside "topol.top"
    cp topol.top topol.bak
    cat topol.top | sed '/forcefield\.itp\"/a\
  #include "ligand.itp"
  ' > topol2.top
    mv topol2.top topol.top
    echo "ligand   1" >> topol.top

#Generate idx group for ligand without hydrogens (for restraints)
ndx=$($GMX make_ndx -f param/ligand/ligand.acpype/ligand_NEW.pdb -o lig_noh.ndx <<EOF
r LIG & !a H*
name 3 LIG-H
q
EOF
)	
    #echo "LIG-H" | $GMX genrestr -f param/ligand/ligand.acpype/ligand_NEW.pdb -o posre_ligand.itp -n lig_noh.ndx -fc 1000 1000 1000
	#copying position restrained
	cp param/ligand/ligand.acpype/posre_ligand.itp .
	
	#Include posre_ligand.itp AT THE END!!!!!!! of  ligand.itp
	echo '
	
	 ; Include Position restraint file
#ifdef POSRES
#include "posre_ligand.itp"
#endif'  >> ligand.itp

    $GMX editconf -f  complex.pdb -o complex_newbox.gro -d $BOXSIZE -bt $BOXTYPE
	PDB=complex
else

	#set "PDB" name (all the simulation filenames are based on this variable).
	PDB=$FILE
	########################
	##   TOPOLOGIE Creation
	#######################
	$GMX pdb2gmx -f $PDB".pdb" -o $PDB"_processed.gro" -water $WATER -ignh -ff $FF
	########################
	##   Solvatation
	#######################
	#Default editconf, changer if you want...
	$GMX editconf -f  $PDB"_processed.gro" -o $PDB"_newbox.gro" -d $BOXSIZE -bt $BOXTYPE
fi



#SOLVATATION
$GMX solvate -cp $PDB"_newbox.gro" -cs spc216.gro -o $PDB"_solv.gro" -p topol.top

#######################
## ADDING IONS
#######################
$GMX grompp -f mdp/ions.mdp -c $PDB"_solv.gro" -p topol.top -o ions.tpr --maxwarn 1

echo "SOL" | $GMX genion -s ions.tpr -o $PDB"_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral

#######################
# MODIFYING FILES IF LIGAND IS DETECTED
#######################
if [ ! -z "$LIGNAME" ]
then
#1. Add constraints to the ligand
    ndx=$($GMX make_ndx -f complex_solv_ions.gro -o index.ndx <<EOF
1 | r LIG
r SOL | r CL | r NA
q
EOF
)

#2. Rename protein+ligand group and Water+ions group
$(python <<EOF
import re
with open('index.ndx', 'r') as file:
    content = file.read()
matches = re.findall(r'\[ \w+ \]', content)
if matches:
    content = content.replace(matches[-1], '[ Water_Ions ]')
    content = content.replace(matches[-2], '[ Protein_Ligand ]')
    with open('index.ndx', 'w') as file:
        file.write(content)
EOF
)

    #Now replace the groups in Equilibration and production mdp files
    sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/nvt_300.mdp
    sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/npt.mdp
    sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/md_prod.mdp


    INDEX="-n index.ndx"
else
INDEX=""

fi

mkdir minimization
cd minimization 
	
cp -R ../mdp .
cp ../$PDB"_solv_ions.gro" .
cp ../topol.top .
cp ../*.itp .
cp ../index.ndx . 2> /dev/null


#######################
## MINIMISATION
#######################
$GMX grompp -f mdp/em.mdp -c $PDB"_solv_ions.gro" -p topol.top -o em.tpr $INDEX
$MDRUNmini -v -deffnm em 

echo "Protein" | $GMX trjconv -f em.gro -s em.tpr -o "${PDB}_minimized.pdb"