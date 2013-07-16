#!/usr/bin/python

from Pegasus.DAX3 import *

import sys
import os

# the top level workflow generator is passing a bunch of stuff to us

base_dir =  sys.argv[1]
id = sys.argv[2]
priority = sys.argv[3]
receptors_fname = sys.argv[4]
receptorDir = sys.argv[5]
ligandpdbqtDir = sys.argv[6]
ligandmol2Dir = sys.argv[7]
subdax_fname = sys.argv[8]

# make lists of all the input files - doing this once and reusing
# the lists speeds up the dax generation process
recFiles = []
ligFiles = []
mol2Files = {}

recList = [ x.strip() for x in open(receptors_fname) ]
for rec in recList:
    recFile = File(rec + '.pdbqt')
    recFile.addPFN(PFN("file://" + receptorDir + "/" + rec + '/' + rec + '.pdbqt', "local"))
    recFiles.append(recFile)

for lig in os.listdir(ligandmol2Dir):
    ligFile = File(lig.replace('.mol2','.pdbqt'))
    ligFile.addPFN(PFN("file://" + ligandpdbqtDir + "/" + lig.replace('.mol2','.pdbqt'), "local"))
    mol2File = File(lig)
    mol2File.addPFN(PFN("file://" + ligandmol2Dir + "/" + lig, "local"))
    ligFiles.append(ligFile)
    mol2Files[ligFile] = mol2File

# Create a abstract dag
subdax = ADAG("vina-" + id)

# Add executables to the DAX-level replica catalog
wrapper = Executable(name="vina_wrapper.sh", arch="x86_64", installed=False)
wrapper.addPFN(PFN("file://" + base_dir + "/vina_wrapper.sh", "local"))
wrapper.addProfile(Profile(Namespace.CONDOR, "priority", priority))
wrapper.addProfile(Profile(Namespace.PEGASUS, "clusters.size", 5))
subdax.addExecutable(wrapper)

# sub-executables (added as input files)
vinaFile = File("vina")
vinaFile.addPFN(PFN("file://" + base_dir + "/vina", "local"))
subdax.addFile(vinaFile)
endscriptFile = File("pdbqt2mol2.py")
endscriptFile.addPFN(PFN("file://" + base_dir + "/pdbqt2mol2.py", "local"))
subdax.addFile(endscriptFile)

innerCount = 0
for recFile in recFiles:
    subdax.addFile(recFile)
    innerCount += 1
    for ligFile in ligFiles:

        mol2File = mol2Files[ligFile]

        # only add the files to the dax one
        if innerCount == 1:
            subdax.addFile(ligFile)
            subdax.addFile(mol2File)

        # load arguments
        #center = open(receptorDir + "/" + recFile.name[:-6] + '/center.txt').readline().strip().split()
        center = open(receptorDir + "/" + recFile.name[:-6] + '/' + recFile.name[:-6] + '.center').readline().strip().split()
        
        # Output file
        outFile = File(recFile.name[:-6] + '-' + ligFile.name[:-6] + '.mol2')

        # Add job to dax
        job = Job(name="vina_wrapper.sh")
        job.addArguments(recFile.name[:-6], ligFile.name[:-6], center[0], center[1], center[2])
        job.uses(vinaFile, link=Link.INPUT)
        job.uses(endscriptFile, link=Link.INPUT)
        job.uses(recFile, link=Link.INPUT)
        job.uses(ligFile, link=Link.INPUT)
        job.uses(mol2File, link=Link.INPUT)
        job.uses(outFile, link=Link.OUTPUT)
        subdax.addJob(job)

# Write the DAX
f = open(base_dir + "/" + subdax_fname, "w")
subdax.writeXML(f)
f.close()

