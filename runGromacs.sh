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
BOXSIZE=1.2 #cubic simulation boxsiwe
BOXTYPE=cubic #Box type
NT=8 #Number of core.
WATER=tip3p #Water type
NUMBEROFREPLICAS=3 #Number of replica
FF=amber99sb-ildn #Force field
SIMULATIONTIME=100 #Simulation time in nanosec. Will be converted in fs and modified in the mdp file.

#---------  HPC SETUP  -----------
MPI="" #If you have to submit jobs with MPI softwares like "mpirun -np 10". Add the command here
GMX=gmx #GMX command (can be "$GMX_mpi" sometimes. Just change it here

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
    echo "LIG-H" | $GMX genrestr -f param/ligand/ligand.acpype/ligand_NEW.pdb -o posre_ligand.itp -n lig_noh.ndx -fc 1000 1000 1000
	
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


 
for ((i=0; i<$NUMBEROFREPLICAS; i++))
	do 
	echo ">>>>> replica_"$((i+1))
	mkdir "replica_"$((i+1))
	cd "replica_"$((i+1))
	cp -R ../mdp .
	cp ../$PDB"_solv_ions.gro" .
	cp ../topol.top .
	cp ../*.itp .


	#######################
	## MINIMISATION
	#######################
	$GMX grompp -f mdp/em.mdp -c $PDB"_solv_ions.gro" -p topol.top -o em.tpr
	$MPI $GMX mdrun -v -deffnm em -nt $NT


	#Cleaning
	mkdir -p results/mini
	mkdir -p gro
	mkdir graph
	mv em* results/mini/
	mv mdout.mdp results/mini/
	mv *.gro gro/
	#potential energy graph
	echo "11 0" | $GMX energy -f results/mini/em.edr -o graph/mini_"$PDB"_pot.xvg





	#######################
	## temperature 300
	#######################
	$GMX grompp -f mdp/nvt_300.mdp -c results/mini/em.gro -r results/mini/em.gro  -p topol.top -o nvt_300.tpr -maxwarn 2
	$MPI $GMX mdrun -deffnm nvt_300 -v -nt $NT
	#temperature_graph

	#cleaning
	mkdir -p results/nvt
	mv nvt* results/nvt/ 2> /dev/null
	mv mdout.mdp results/nvt/

	#Temparture graph
	echo "16 0" | $GMX energy -f results/nvt/nvt_300.edr -o graph/temperature_nvt_300.xvg

	#######################
	## Pression
	#######################
	$GMX grompp -f mdp/npt.mdp -c results/nvt/nvt_300.gro -r results/nvt/nvt_300.gro -t results/nvt/nvt_300.cpt -p topol.top -o npt_ab.tpr -maxwarn 2
	$MPI $GMX mdrun -deffnm npt_ab -v -nt $NT

	#cleaning
	mkdir -p results/npt
	mv npt* results/npt/ 2> /dev/null
	mv mdout.mdp results/npt_ab/
	#Pression and density graph
	echo "17 0" | $GMX energy -f results/npt/npt_ab.edr -o graph/npt_"$PDB"_pressure.xvg
	echo "22 0" | $GMX energy -f results/npt/npt_ab.edr -o graph/npt_"$PDB"_volume.xvg



	#######################
	## Production 
	#######################
	
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
	
	$GMX grompp -f mdp/md_prod.mdp -c results/npt/npt_ab.gro -t results/npt/npt_ab.cpt -p topol.top -o "md_"$PDB"_prod.tpr" -maxwarn 2
	$MPI $GMX mdrun -deffnm "md_"$PDB"_prod"  -v -nt $NT

	mkdir -p results/prod
	mv md_* results/prod 2> /dev/null
	mv mdout.mdp results/prod/

	echo "backbone backbone" | $GMX rms -s "results/prod/md_"$PDB"_prod.tpr" -f "results/prod/md_"$PDB"_prod.trr" -o graph/prod_"$PDB"_rmsd.xvg -tu ns

	cd results/prod

	echo "Protein System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_prod.trr" -o "md_"$PDB"_clean_temp.xtc" -pbc nojump -ur compact -center

	echo "Protein System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_full.xtc" -fit rot+trans

	echo "Protein non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.xtc" -fit rot+trans

	echo "Protein non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.pdb" -pbc nojump -ur compact -center -b 0 -e 0

	rm "md_"$PDB"_clean_temp.pdb"
	
	echo "non-Water" | $GMX convert-tpr -s "md_"$PDB"_prod.tpr" -o tpr_nowat.tpr

	cd ../../
	#Final graph
	echo "backbone backbone" | $GMX rms -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsd_ca.xvg" -tu ns
	echo "protein protein" | $GMX rms -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsd_all.xvg" -tu ns
	echo "backbone" | $GMX gyrate -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_gyrate.xvg"

	echo "backbone" | $GMX rmsf -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsf_ref.xvg" -res
ls
	#####
	# TESTING DSSP Installation
	#####
	export DSSP=`which dssp`
	if [ -z "$DSSP" ]
	then
	      export DSSP=`which mkdssp`
	else
		echo "DSSP found. Good!"
	fi

	if [ -z "$DSSP" ]
	then
	      echo "DSSP SOFTWARE NOT FOUND. If you have anaconda or miniconda install, please install it with this command line:"
		  echo "conda install -c conda-forge -c salilab dssp"
		  echo "note that the name will be mkdssp"
		  echo "Trying to install it right now..."
		  conda install -y -c conda-forge -c salilab dssp
		  export DSSP=`which mkdssp`
		if [ -z "$DSSP" ]
			then
			      echo "Installation faillure..."
			else
				echo "1" |  $GMX do_dssp -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_ss.xpm" -tu ns -dt 0.05 -ver 3
				$GMX xpm2ps -f "graph/prod_"$PDB"_ss.xpm" -o "graph/prod_"$PDB"_ss.ps" -by 10 -bx 3
			fi
	else
		echo "Protein" | $GMX do_dssp -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_ss.xpm" -ver 3 -tu ns -dt 0.05
		$GMX xpm2ps -f "graph/prod_"$PDB"_ss.xpm" -o "graph/prod_"$PDB"_ss.ps" -by 10 -bx 3
	fi

	cd ../
done
