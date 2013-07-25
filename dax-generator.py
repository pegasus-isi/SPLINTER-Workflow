#!/usr/bin/python

from Pegasus.DAX3 import *

import sys
import os

base_dir = os.getcwd()

##############################

## Number of receptors to group together in a sub workflow
num_recs_per_sub_wf = 1

## Location of receptor PDBQTs
receptor_dir = base_dir + '/inputs/rec'

## Location of ligand PDBQTs and MOL2s
ligand_dir = base_dir + '/inputs/lig'

##############################


def splitlist(lst, slicelen):
    for i in range(0,len(lst),slicelen):
        yield lst[i:i + slicelen]


def add_subwf(dax, id, id_formatted, receptors):

    priority = 1000 - id
    receptors_fname = "receptors-%s.txt" % id_formatted
    subdax_fname = "dax-%s.xml" % id_formatted
    
    # create the receptors file
    f = open(receptors_fname, "w")
    for rec in receptors:
        f.write(rec)
        f.write("\n")
    f.close()

    receptors_file = File(receptors_fname)
    receptors_file.addPFN(PFN("file://%s/%s" % (base_dir, receptors_fname), "local"))
    dax.addFile(receptors_file)
    
    subdax_file = File(subdax_fname)
    subdax_file.addPFN(PFN("file://%s/%s" % (base_dir, subdax_fname), "local"))
    dax.addFile(subdax_file)

    # job to generate the subdax
    subdax_gen = Job(name="subdax-generator.py")
    subdax_gen.addArguments(base_dir,
                            "%d" % id,
                            "%d" % priority,
                            receptors_fname,
                            receptor_dir,
                            ligand_dir,
                            subdax_fname)
    subdax_gen.uses(receptors_file, link=Link.INPUT)
    subdax_gen.addProfile(Profile("dagman", "PRIORITY", "%d" % (priority)))
    subdax_gen.addProfile(Profile("hints", "executionPool", "local"))
    subdax_gen.addProfile(Profile("env", "PATH", os.environ['PATH']))
    subdax_gen.addProfile(Profile("env", "PYTHONPATH", os.environ['PYTHONPATH']))
    subdax_gen.invoke("on_error", "/usr/share/pegasus/notification/email")
    dax.addJob(subdax_gen)

    # job to run subwf
    subwf = DAX(subdax_fname, id="sub-%s" % (id_formatted))
    subwf.addArguments("-Dpegasus.catalog.site.file=%s/sites.xml" % (base_dir),
                       "--cluster", "horizontal",
                       "--sites", "condorpool",
                       "--basename", id_formatted,
                       "--force",
                       "--output-site", "local")
    subwf.uses(subdax_file, link=Link.INPUT, register=False)
    subwf.addProfile(Profile("dagman", "PRIORITY", "%d" % (priority)))
    subwf.addProfile(Profile("dagman", "CATEGORY", "subworkflow"))
    subwf.invoke("on_error", "/usr/share/pegasus/notification/email")
    dax.addDAX(subwf)
    dax.depends(parent=subdax_gen, child=subwf)



# build a list of the receptors we want to process
recList = []
for rec in os.listdir(receptor_dir):
    recList.append(rec)

# top level workflow
dax = ADAG("splinter")
    
# notifcations on state changes for the dax
dax.invoke("all", "/usr/share/pegasus/notification/email")
    
# Add executables to the DAX-level replica catalog
subdax_generator = Executable(name="subdax-generator.py", arch="x86_64", installed=False)
subdax_generator.addPFN(PFN("file://" + base_dir + "/subdax-generator.py", "local"))
dax.addExecutable(subdax_generator)

subwfID = 1
for recSubList in splitlist(recList, num_recs_per_sub_wf):
    print "  generarating subworkflow %06d" % (subwfID)
    add_subwf(dax, subwfID, "%06d" % subwfID, recSubList)
    subwfID += 1

# Write the DAX
f = open("dax-top.xml", "w")
dax.writeXML(f)
f.close()


