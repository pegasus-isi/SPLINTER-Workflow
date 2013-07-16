#!/usr/bin/python
"""
Usage: ./pdbqt.py <old.pdbqt> <old.mol2> <new.mol2>
"""
def randString():
   from random import choice
   import string
   chars = string.letters + string.digits
   newstring = ''
   for i in range(8):
       newstring += choice(chars)
   return newstring

class pdbqt2mol2:
   """
   Loads values from PDBQT file into a mol2 file and gives new mol2 file.
   """
   def __init__(self,oldPDBQTfilename,oldmol2filename,newmol2filename):

      import os,tempfile,sys,shutil

      tempmol2filename = '/tmp/pdbqt2mol2_'+randString()

      PDBQTlines = open(oldPDBQTfilename).readlines()
      atomlocations = {}                # new atom coordinates
      remapatoms = {}                   # new atom number, for bond remap
      adscore = ''
      for line in PDBQTlines:
         if line.startswith('REMARK VINA RESULT:'):
            adscore,rmsd1,rmsd2 = line.split(':')[1].split()[:3]
         if line.find('ATOM') == 0:
            atomname = line[11:16].split()[0]
            xcoord = line[30:38].split()[0]
            ycoord = line[38:46].split()[0]
            zcoord = line[46:54].split()[0]

            atomlocations[atomname] = [xcoord, ycoord, zcoord]

      oldmol2 = open(oldmol2filename)
      newmol2 = open(tempmol2filename, 'w')

      atomnum = 0
      bondnum = 0

      counter = 0
      start_mol = False
      while 1:
         if counter == 0 and adscore != '':
            print >> newmol2,'### USER    Estimated Free Energy of Binding    =   %s kcal/mol'%(adscore,)
         counter += 1
         line = oldmol2.readline()

         if line.startswith('@<TRIPOS>MOLECULE'): start_mol = True
         if start_mol:
            newmol2.write(line)
         else:
            continue
         if line.startswith('@<TRIPOS>ATOM'): break

      oldFilePos = 'New'
      newFilePos = oldmol2.tell()
      while oldFilePos != newFilePos:
         line = oldmol2.readline()
         oldFilePos = newFilePos
         newFilePos = oldmol2.tell()
         if oldFilePos == newFilePos: break

         if line.find('<TRIPOS>BOND') != -1:
            newmol2.write(line)
            break
         contents = line.split()
 
        # common @<TRIPOS>ATOM line: 
        #      1 S1          6.5100   -2.3824    0.4691 S.o2      1 LIG         0.8327 
        #    []0 1           2        3          4      5         6 7           8
 
         try:
            newcoords = atomlocations[contents[1]]
            atomnum += 1
            remapatoms[contents[0]] = str(atomnum)
            modline = '  ' + ' '.join([str(atomnum), contents[1],
                                       newcoords[0], newcoords[1],
                                       newcoords[2], contents[5],
                                       contents[6], contents[7],
                                       contents[8]]) + '\n'
            llist = modline.split()
            nxx = float(llist[2])
            nyy = float(llist[3])
            nzz = float(llist[4])
            achg = float(llist[8])
            newmol2.write("%7s %-8s%10.4f%10.4f%10.4f %-8s%3s %-8s%10.4f\n"%(llist[0],llist[1],nxx,nyy,nzz,llist[5],llist[6],llist[7],achg))
         except KeyError:
            continue

      while oldFilePos != newFilePos:
         line = oldmol2.readline()
         oldFilePos = newFilePos
         newFilePos = oldmol2.tell()
         if oldFilePos == newFilePos: break

         if line.find('<TRIPOS>SUBSTRUCTURE') != -1:
            newmol2.write(line)
            break
         try:
            bondfirst = line.split()[1]
            bondsecond = line.split()[2]
         except IndexError:
            newmol2.write(line)
            continue

        # common @<TRIPOS>BOND line:
        #     5    4    6 ar
        #   []0    1    2 3

         try:
            newfirst = remapatoms[bondfirst]
            newsecond = remapatoms[bondsecond]
         except KeyError:
            continue
         bondnum += 1
         modline = '  ' + ' '.join([str(bondnum), newfirst, newsecond,
                                   line.split()[3]]) + '\n'
         llist = modline.split()
         newmol2.write(" %5s%5s%5s %-5s\n"%(llist[0],llist[1],llist[2],llist[3]))

      while oldFilePos != newFilePos:
         line = oldmol2.readline()
         oldFilePos = newFilePos
         newFilePos = oldmol2.tell()
         if oldFilePos == newFilePos: break
         newmol2.write(line)

      oldmol2.close()
      newmol2.flush()
      newmol2.close()

      # Open it right back up -- fix the header
      oldmol2 = open(tempmol2filename)
      oldFilePos = 'New'
      newFilePos = oldmol2.tell()
      newmol2 = open(newmol2filename, 'w')

      while oldFilePos != newFilePos:
         line = oldmol2.readline()
         oldFilePos = newFilePos
         newFilePos = oldmol2.tell()

        # common header lines:
        # @<TRIPOS>MOLECULE
        # 10113978a
        # 55   58    1
        # SMALL
        # USER_CHARGES

         if line.find('<TRIPOS>MOLECULE') != -1:
            newmol2.write(line)
            newmol2.write(oldmol2.readline())
            line = oldmol2.readline()
            try:
               tp1,tp2,tp3 = line.split()[2:5]
               tp1 = int(tp1)
               tp2 = int(tp2)
               tp3 = int(tp3)
            except:
               tp1,tp2,tp3 = 1,0,0
            newmol2.write('%5d%6d%6d%6d%6d\n'%(atomnum,bondnum,tp1,tp2,tp3,))
            break
         newmol2.write(line)

      while oldFilePos != newFilePos:
         line = oldmol2.readline()
         oldFilePos = newFilePos
         newFilePos = oldmol2.tell()
         newmol2.write(line)

      oldmol2.close()
      newmol2.flush()
      newmol2.close()
      os.unlink(tempmol2filename)

def splitpdbqt(pdbqt):
   f = open(pdbqt,'r')
   fout = open(pdbqt + '-new','w')
   f.next()
   for line in f:
      if line.startswith('ENDMDL'):
         break
      fout.write(line)
   f.close()
   fout.close()
   return pdbqt + '-new'

if __name__ == '__main__':
   import os,sys
   if len(sys.argv) != 4:
      print>>sys.stderr, __doc__
      raise SystemExit
   oldPDBQT = sys.argv[1]
   oldMOL2 = sys.argv[2]
   newMOL2 = sys.argv[3]
   oldPDBQT = splitpdbqt(oldPDBQT)
   a = pdbqt2mol2(oldPDBQT,oldMOL2,newMOL2)
