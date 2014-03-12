#!/usr/bin/python

from Pegasus.DAX3 import *

import sys
import os
from stat import *

# the top level workflow generator is passing a bunch of stuff to us
base_dir =  sys.argv[1]
id = int(sys.argv[2])
priority = sys.argv[3]
work_fname = sys.argv[4]
subdax_fname = sys.argv[5]

receptors_pdbqt_saved = {}
receptors_mol2_saved = {}
ligands_pdbqt_saved = {}
ligands_mol2_saved = {}

# Create a abstract dag
subdax = ADAG("splinter-%06d" % id)

# Managed stageout jobs
pt = Executable(namespace="pegasus", name="transfer", arch="x86_64", installed=True)
pt.addPFN(PFN("/usr/bin/pegasus-transfer", "local"))
pt.addProfile(Profile(Namespace.PEGASUS, "style", "condor"))
pt.addProfile(Profile(Namespace.CONDOR, "universe", "vanilla"))
pt.addProfile(Profile(Namespace.CONDOR, "requirements", "TARGET.FileSystemDomain =?= \"osg-xsede.grid.iu.edu\""))
pt.addProfile(Profile(Namespace.CONDOR, "+RunOnSubmitNode", "True"))
subdax.addExecutable(pt)

# Add executables to the DAX-level replica catalog
wrapper = Executable(name="vina_wrapper.sh", arch="x86_64", installed=False)
wrapper.addPFN(PFN("file://" + base_dir + "/vina_wrapper.sh", "local"))
wrapper.addProfile(Profile(Namespace.CONDOR, "priority", priority))
wrapper.addProfile(Profile(Namespace.PEGASUS, "clusters.size", 40))
subdax.addExecutable(wrapper)

# sub-executables (added as input files)
vina_file = File("vina")
vina_file.addPFN(PFN("file://" + base_dir + "/vina", "local"))
subdax.addFile(vina_file)
endscriptFile = File("pdbqt2mol2.py")
endscriptFile.addPFN(PFN("file://" + base_dir + "/pdbqt2mol2.py", "local"))
subdax.addFile(endscriptFile)
endscript2File = File("genfet.py")
endscript2File.addPFN(PFN("file://" + base_dir + "/genfet.py","local"))
subdax.addFile(endscript2File)
endscript2datFile = File("KB_stat_76.dat")
endscript2datFile.addPFN(PFN("file://" + base_dir + "/KB_stat_76.dat","local"))
subdax.addFile(endscript2datFile)

work_file = open(work_fname)
for line in work_file:
    line = line.strip()
    rec_name, rec_location, lig_name, lig_pdbqt_location, lig_mol2_location = line.split("\t")
        
    # only add the files to the dax one
    if not rec_name in receptors_pdbqt_saved:
        f = File(rec_name + '.pdbqt')
        f.addPFN(PFN("file://" + rec_location + '/' + rec_name + '.pdbqt', "local"))
        subdax.addFile(f)
        receptors_pdbqt_saved[rec_name] = f
        #Add mol2 rec file
        f = File(rec_name + '.mol2')
        f.addPFN(PFN("file://" + rec_location + '/' + rec_name + '.mol2', "local"))
        subdax.addFile(f)
        receptors_mol2_saved[rec_name] = f

    if not lig_name in ligands_pdbqt_saved:
        f = File(lig_name + '.pdbqt')
        f.addPFN(PFN("file://" + lig_pdbqt_location, "local"))
        subdax.addFile(f)
        ligands_pdbqt_saved[lig_name] = f
        # and also the mol2 file
        f = File(lig_name + '.mol2')
        f.addPFN(PFN("file://" + lig_mol2_location, "local"))
        subdax.addFile(f)
        ligands_mol2_saved[lig_name] = f

    # load arguments (center.txt or [rec].center)
    center = None
    if os.path.isfile(rec_location + '/center.txt'):
        center = open(rec_location + '/center.txt').readline().strip().split()
    else:
        center = open(rec_location + '/' + rec_name + '.center').readline().strip().split()

    # sometimes the first column is the receptor name
    if rec_name == center[0]:
        center[0] = center[1]
        center[1] = center[2]
        center[2] = center[3]
    
    # Output file
    out_file = File("%s-%s.mol2" % (rec_name, lig_name))

    # Output feature file
    out_file2 = File("%s-%s.dat" % (rec_name, lig_name))

    # Add job to dax
    job = Job(name="vina_wrapper.sh")
    job.addArguments("%06d" % id, rec_name, lig_name, center[0], center[1], center[2])
    job.uses(vina_file, link=Link.INPUT)
    job.uses(endscriptFile, link=Link.INPUT)
    job.uses(endscript2File, link=Link.INPUT)
    job.uses(endscript2datFile, link=Link.INPUT)
    job.uses(receptors_pdbqt_saved[rec_name], link=Link.INPUT)
    job.uses(receptors_mol2_saved[rec_name], link=Link.INPUT)
    job.uses(ligands_pdbqt_saved[lig_name], link=Link.INPUT)
    job.uses(ligands_mol2_saved[lig_name], link=Link.INPUT)
    job.uses(out_file, link=Link.OUTPUT)
    job.uses(out_file2, link=Link.OUTPUT)
    subdax.addJob(job)

# Write the DAX
f = open(base_dir + "/" + subdax_fname, "w")
subdax.writeXML(f)
f.close()

