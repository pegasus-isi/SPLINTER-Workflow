#!/usr/bin/python

from Pegasus.DAX3 import *

import sys
import os
from stat import *

# the top level workflow generator is passing a bunch of stuff to us
base_dir =  sys.argv[1]
id = sys.argv[2]
priority = sys.argv[3]
receptors_fname = sys.argv[4]
receptor_dir = sys.argv[5]
ligand_dir = sys.argv[6]
subdax_fname = sys.argv[7]

# globals
rec_files = []
lig_pdbqt_files = []
lig_mol2_files = {}

def find_ligands(dir):
    global lig_pdbqt_files
    global lig_mol2_files
    for entry in os.listdir(dir):
        path = os.path.join(dir, entry)
        mode = os.stat(path).st_mode
        if S_ISDIR(mode):
            find_ligands(path)
        else:
            f_base, f_ext = os.path.splitext(entry)
            if f_ext == ".pdbqt":
                lig_pdbqt_file = File(entry)
                lig_pdbqt_file.addPFN(PFN("file://" + path, "local"))
                lig_pdbqt_files.append(lig_pdbqt_file)
            elif f_ext == ".mol2":
                mol2_file = File(entry)
                mol2_file.addPFN(PFN("file://" + path, "local"))
                lig_mol2_files[f_base] = mol2_file
    
# make lists of all the input files - doing this once and reusing
# the lists speeds up the dax generation process
recList = [ x.strip() for x in open(receptors_fname) ]
for rec in recList:
    rec_file = File(rec + '.pdbqt')
    rec_file.addPFN(PFN("file://" + receptor_dir + "/" + rec + '/' + rec + '.pdbqt', "local"))
    rec_files.append(rec_file)

find_ligands(ligand_dir)

# Create a abstract dag
subdax = ADAG("splinter-" + id)

# Add executables to the DAX-level replica catalog
wrapper = Executable(name="vina_wrapper.sh", arch="x86_64", installed=False)
wrapper.addPFN(PFN("file://" + base_dir + "/vina_wrapper.sh", "local"))
wrapper.addProfile(Profile(Namespace.CONDOR, "priority", priority))
wrapper.addProfile(Profile(Namespace.PEGASUS, "clusters.size", 5))
subdax.addExecutable(wrapper)

# sub-executables (added as input files)
vina_file = File("vina")
vina_file.addPFN(PFN("file://" + base_dir + "/vina", "local"))
subdax.addFile(vina_file)
endscriptFile = File("pdbqt2mol2.py")
endscriptFile.addPFN(PFN("file://" + base_dir + "/pdbqt2mol2.py", "local"))
subdax.addFile(endscriptFile)

innerCount = 0
for rec_file in rec_files:
    subdax.addFile(rec_file)
    innerCount += 1
    for lig_pdbqt_file in lig_pdbqt_files:

        mol2_file = lig_mol2_files[lig_pdbqt_file.name[:-6]]

        # only add the files to the dax one
        if innerCount == 1:
            subdax.addFile(lig_pdbqt_file)
            subdax.addFile(mol2_file)

        # load arguments
        center = open(receptor_dir + "/" + rec_file.name[:-6] + '/' + rec_file.name[:-6] + '.center').readline().strip().split()
        
        # Output file
        out_file = File(rec_file.name[:-6] + '-' + lig_pdbqt_file.name[:-6] + '.mol2')

        # Add job to dax
        job = Job(name="vina_wrapper.sh")
        job.addArguments(rec_file.name[:-6], lig_pdbqt_file.name[:-6], center[0], center[1], center[2])
        job.uses(vina_file, link=Link.INPUT)
        job.uses(endscriptFile, link=Link.INPUT)
        job.uses(rec_file, link=Link.INPUT)
        job.uses(lig_pdbqt_file, link=Link.INPUT)
        job.uses(mol2_file, link=Link.INPUT)
        job.uses(out_file, link=Link.OUTPUT)
        subdax.addJob(job)

# Write the DAX
f = open(base_dir + "/" + subdax_fname, "w")
subdax.writeXML(f)
f.close()

