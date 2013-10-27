#!/usr/bin/python

from Pegasus.DAX3 import *

import sys
import os
from stat import *

##############################
## Settings

# Maximum tasks in a sub workflow
max_tasks_per_sub_wf = 10000

##############################

base_dir = os.getcwd()

receptor_dir = sys.argv[1]
ligand_dir = sys.argv[2]

# globals
rec_names = []  
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
                lig_pdbqt_files.append(path)
            elif f_ext == ".mol2":
                lig_mol2_files[f_base] = path


def splitlist(lst, slicelen):
    for i in range(0,len(lst),slicelen):
        yield lst[i:i + slicelen]


def add_subwf(dax, id):

    priority = 1000 - id
    subdax_fname = "dax-%06d.xml" % id

    work_file = File("work-%06d.txt" % id)
    work_file.addPFN(PFN("file://%s/work-%06d.txt" % (base_dir, id), "local"))
    dax.addFile(work_file)
    
    subdax_file = File(subdax_fname)
    subdax_file.addPFN(PFN("file://%s/%s" % (base_dir, subdax_fname), "local"))
    dax.addFile(subdax_file)

    # job to generate the subdax
    subdax_gen = Job(name="subdax-generator.py")
    subdax_gen.addArguments(base_dir,
                            "%d" % id,
                            "%d" % priority,
                            "work-%06d.txt" % (id),
                            subdax_fname)
    subdax_gen.uses(work_file, link=Link.INPUT)
    subdax_gen.addProfile(Profile("dagman", "PRIORITY", "%d" % (priority)))
    subdax_gen.addProfile(Profile("hints", "executionPool", "local"))
    subdax_gen.addProfile(Profile("env", "PATH", os.environ['PATH']))
    subdax_gen.addProfile(Profile("env", "PYTHONPATH", os.environ['PYTHONPATH']))
    subdax_gen.invoke("on_error", "/usr/share/pegasus/notification/email")
    dax.addJob(subdax_gen)

    # job to run subwf
    subwf = DAX(subdax_fname, id="sub-%06d" % (id))
    subwf.addArguments("-Dpegasus.catalog.site.file=%s/sites.xml" % (base_dir),
                       "--cluster", "horizontal",
                       "--sites", "condorpool",
                       "--basename", "%06d" % id,
                       "--force",
                       "--force-replan",
                       "--output-site", "local")
    subwf.uses(subdax_file, link=Link.INPUT, register=False)
    subwf.addProfile(Profile("dagman", "PRIORITY", "%d" % (priority)))
    subwf.addProfile(Profile("dagman", "CATEGORY", "subworkflow"))
    subwf.invoke("on_error", "/usr/share/pegasus/notification/email")
    dax.addDAX(subwf)
    dax.depends(parent=subdax_gen, child=subwf)



# build a list of the receptors and ligands we want to process
for rec in os.listdir(receptor_dir):
    rec_names.append(rec)
find_ligands(ligand_dir)

# top level workflow
dax = ADAG("splinter")
    
# notifcations on state changes for the dax
dax.invoke("all", "/usr/share/pegasus/notification/email")
    
# Add executables to the DAX-level replica catalog
subdax_generator = Executable(name="subdax-generator.py", arch="x86_64", installed=False)
subdax_generator.addPFN(PFN("file://" + base_dir + "/subdax-generator.py", "local"))
dax.addExecutable(subdax_generator)

subwf_id = 1
subwf_task_count = 0
work_file = open("work-%06d.txt" %(1), 'w')
for rec_name in rec_names:
 
    rec_location = receptor_dir + "/" + rec_name

    for lig_pdbqt_file in lig_pdbqt_files:

        lig_name = os.path.splitext(os.path.basename(lig_pdbqt_file))[0]
    
        if subwf_task_count > max_tasks_per_sub_wf:
            work_file.close()
            print "  generarating subworkflow %06d" % (subwf_id)
            add_subwf(dax, subwf_id)
            
            subwf_id += 1
            work_file = open("work-%06d.txt" % subwf_id, 'w')
            subwf_task_count = 0

        subwf_task_count += 1
        work_file.write("%s\t%s\t%s\t%s\t%s\n" %(rec_name, rec_location, lig_name, lig_pdbqt_file, lig_mol2_files[lig_name]))

if subwf_task_count > 0:
    work_file.close()
    print "  generarating subworkflow %06d" % (subwf_id)
    add_subwf(dax, subwf_id)

# Write the DAX
f = open("dax-top.xml", "w")
dax.writeXML(f)
f.close()


