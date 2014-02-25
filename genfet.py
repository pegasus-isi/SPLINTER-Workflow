#!/usr/bin/env python
"""
    Copyright (c) 2012, Indiana University School of Medicine
    All rights reserved.

    Author: Liwei Li
    Updated 2014/02/20: Combine into 1 file and fixes for OSG-XSEDE
    Last Updated Date: 02/20/2014
"""
import sys
import os
import math
from collections import deque

"""Class ChemKit"""
def atomic_distance2(aAtom,bAtom):
    try:
       return (aAtom.x-bAtom.x)*(aAtom.x-bAtom.x)+(aAtom.y-bAtom.y)*(aAtom.y-bAtom.y)+(aAtom.z-bAtom.z)*(aAtom.z-bAtom.z)
    except:
       return None

def atomic_distance(aAtom,bAtom):
    dist2 = atomic_distance2(aAtom,bAtom)
    if dist2 != None:
       return math.sqrt(dist2)
    else:
       return None

class Atom(object):
    """Atom class"""
    def __init__(self,ax=0.0,ay=0.0,az=0.0,atp=None,id=None,sb=None,res=None,custom=None):
        self.x = float(ax)
        self.y = float(ay)
        self.z = float(az)
        self.type = atp
        self.id = None
        if id: self.id = int(id)
        self.symbol = sb
        self.resname = res
        self.connlist = []
        self.ctm = custom

    def setXYZ(self,ax,ay,az):
        self.x = float(ax)
        self.y = float(ay)
        self.z = float(az)

    def getXYZ(self):return [self.x,self.y,self.z]

    def setType(self,atp):
        self.type = atp

    def getType(self):return self.type

    def setSymbol(self,sb):
        self.symbol = sb

    def setResName(self,res):
        self.resname = res

    def setID(self,id):
        self.id = int(id)

    def getID(self,id): return self.id
      
    def getMass(self):
        return Atom.__massDict[self.symbol]

    def getUnitedRadius(self):
        return Atom.__radiusUnitedDict[self.type.lower()]

    def updateConnlist(self,newlist):
        self.connlist = newlist[:]

    def addToConnlist(self,newid):
        self.connlist.append(newid)

    __massDict = {
    'H'  :   1.00794,
    'Li' :   6.941,
    'LI' :   6.941,
    'B'  :  10.811,
    'C'  :  12.0107,
    'N'  :  14.0067,
    'O'  :  15.9994,
    'F'  :  18.9984032,
    'Na' :  22.989770,
    'NA' :  22.989770,
    'Mg' :  24.3050,
    'MG' :  24.3050,
    'Al' :  26.981538,
    'AL' :  26.981538,
    'Si' :  28.0855,
    'SI' :  28.0855,
    'P'  :  30.973761,
    'S'  :  32.065,
    'Cl' :  35.453,
    'CL' :  35.453,
    'K'  :  39.0983,
    'Ca' :  40.078,
    'CA' :  40.078,
    'Cr' :  51.9961,
    'CR' :  51.9961,
    'Mn' :  54.938049,
    'MN' :  54.938049,
    'Fe' :  55.845,
    'FE' :  55.845,
    'Co' :  58.933200,
    'CO' :  58.933200,
    'Ni' :  58.6934,
    'NI' :  58.6934,
    'Cu' :  63.546,
    'CU' :  63.546,
    'Zn' :  65.39,
    'ZN' :  65.39,
    'Br' :  79.904,
    'BR' :  79.904,
    'Mo' :  95.94,
    'MO' :  95.94,
    'Ag' : 107.8682,
    'AG' : 107.8682,
    'I'  : 126.90447,
    'Hg' : 200.59,
    'HG' : 200.59,
    }

    __radiusUnitedDict = {
    'h'    :  0.000,
    'c.3'  :  1.908,
    'c.2'  :  1.908,
    'c.ar' :  1.908,
    'c.cat':  1.908,
    'n.4'  :  1.824,
    'n.am' :  1.824,
    'n.pl3':  1.824,
    'n.2'  :  1.824,
    'n.ar' :  1.824,
    'o.3'  :  1.7210,
    'o.2'  :  1.6612,
    'o.co2':  1.6612,
    'p.3'  :  2.10,
    's.3'  :  2.00,
    'f'    :  1.75,
    'cl'   :  1.948,
    'br'   :  2.02,
    'i'    :  2.150,
    'met'  :  2.00,
    'du'   :  1.90
    }
    
class Molecule(object):
    """class Molecule"""
    def __init__(self):
        self.molID = None
        self.atoms = []
        self.atoms = deque(self.atoms)
        pass

    def __del__(self):
        self.atoms.clear()

    def addAtom(self,oneAtom):
        if isinstance(oneAtom,Atom.Atom):
            self.atoms.append(oneAtom)
        else:
            print >> sys.stderr,'Only atoms can be added into a molecule'

    def delAtom(self,oneAtom):
        if isinstance(oneAtom,Atom):
            self.atoms.remove(oneAtom)
        else:
            print >> sys.stderr,'Only atoms can be added into a molecule'

    def addBond(self,atom1ID,atom2ID):
        index = 0
        for atom in self.atoms:
            if atom.id == atom1ID:
               self.atoms[index].addToConnlist(atom2ID)
            elif atom.id == atom2ID:
               self.atoms[index].addToConnlist(atom1ID)
            index += 1

    def hasAtom(self,aAtom):
        for bAtom in self.atoms:
            if (math.fabs(aAtom.x - bAtom.x)<0.0001 and math.fabs(aAtom.y - bAtom.y)<0.0001 and math.fabs(aAtom.z - bAtom.z)<0.0001 
                and (aAtom.id == bAtom.id) and (aAtom.id != None)):
                return True
        return False

class SMmolecule(Molecule):
    def addAtom(self,oneAtom):
        atp = oneAtom.type.lower()
        if atp[0] == 'h': return
        if atp == 'lp': return
        sb = atp.split('.')[0].upper()
        if len(sb) == 2: sb = sb[0]+sb[1].lower()
        if atp in ['c.1']: atp = 'c.ar'
        if atp in ['n.3']: atp = 'n.4'
        if atp in ['s.o','s.o2']: atp = 'p.3'
        if atp in ['s.2']: atp = 'o.2'
        if atp in ['n.1']: atp = 'n.2'
        if atp in ['br','i']: atp = 'cl'
        oneAtom.type = atp
        oneAtom.symbol = sb
        self.atoms.append(oneAtom)
        return

class PTmolecule(Molecule):
    def addAtom(self,oneAtom):
        atp = oneAtom.type.upper()
        if atp[0] == 'H': return
        if atp == 'LP': return
        sb = atp.split('.')[0].upper()
        if len(sb) == 2: sb = sb[0]+sb[1].lower()
        if atp in ['S.O','S.O2',]: atp = 'P.3'
        if atp in ['C.1']: atp = 'C.AR'
        if atp in ['S.2']: atp = 'O.2'
        if atp in ['N.3']: atp = 'N.4'
        if atp in ['N.1','N.AR']: atp = 'N.2'
        if atp in ['CA','C0','FE','ZN','MN','MG','NI','HG','CO','CO.OH']: atp = 'MET'
        oneAtom.type = atp
        oneAtom.symbol = sb
        self.atoms.append(oneAtom)
        return

class Molecules(object):
    """ class Molecules """
    def __init__(self):
        self.molecules = []
        self.molecules = deque(self.molecules)

    def __del__(self):
        self.molecules.clear()

    def addMolecule(self,molecule):
        self.molecules.append(molecule)

    def readCompdFromMol2(self,mol2Name):
        mollines = open(mol2Name)
        inmol = False
        inatom = False
        inbond = False
        midflag = False
        count = 1
        for line in mollines:
            if inmol and inbond:
                if line.startswith('#') or line.startswith('@') or line.strip() == '':
                    inmol = False
                    self.addMolecule(newmol)
                    del newmol
            if line.startswith('@<TRIPOS>MOLECULE'):
                inmol = True
                inatom = False
                inbond = False
                newmol = SMmolecule()
                midflag = True
            if line.startswith('@<TRIPOS>BOND'):
                inbond = True
                inatom = False
            if line.startswith('@<TRIPOS>ATOM'):
                inatom = True

            if line.startswith('@<TRIPOS>'): continue

            if midflag:
                mid = line.split()[0]
                if mid.find('*') != -1: mid = 'TMP'+str(count)
                newmol.molID = mid
                midflag = False
                continue

            if inmol and inatom:
                atomid,tmp,ax,ay,az,atp,resi,res = line.split()[:8]
                sb = atp.split('.')[0].upper()
                newatom = Atom(ax,ay,az,atp,atomid,sb,res)
                newmol.addAtom(newatom)
            elif inmol and inbond:
                bdid,iatom1,iatom2,bdtype = line.split()[:4]
                iatom1 = int(iatom1)
                iatom2 = int(iatom2)
                newmol.addBond(iatom1,iatom2)
            else:
                continue
        if inmol:
            self.addMolecule(newmol)
            del newmol
            inmol = False
            count += 1

    def readProtFromMol2(self,mol2Name):
        mollines = open(mol2Name)
        inmol = False
        inatom = False
        inbond = False
        midflag = False
        for line in mollines:
            if inmol and inbond:
                if line.startswith('#') or line.startswith('@') or line.strip() == '':
                    inmol = False
                    self.addMolecule(newmol)
                    del newmol
            if line.startswith('@<TRIPOS>MOLECULE'):
                inmol = True
                inatom = False
                inbond = False
                newmol = PTmolecule()
                midflag = True
            if line.startswith('@<TRIPOS>BOND'):
                inbond = True
                inatom = False
            if line.startswith('@<TRIPOS>ATOM'):
                inatom = True

            if line.startswith('@<TRIPOS>'): continue

            if midflag:
                mid = line.split()[0]
                if mid.find('*') != -1: mid = mol2Name.split(os.sep)[-1].split('.mol2')[0]
                newmol.molID = mid
                midflag = False
                continue

            if inmol and inatom:
                atomid,tmp,ax,ay,az,atp,resi,res = line.split()[:8]
                sb = atp.split('.')[0].upper()
                newatom = Atom(ax,ay,az,atp,atomid,sb,res)
                newmol.addAtom(newatom)
            elif inmol and inbond:
                bdid,iatom1,iatom2,bdtype = line.split()[:4]
                iatom1 = int(iatom1)
                iatom2 = int(iatom2)
                newmol.addBond(iatom1,iatom2)
            else:
                continue
        if inmol:
            self.addMolecule(newmol)
            del newmol
            inmol = False

""" Compute solvent accessible surface area"""
class SASA(object):
    """class of SASA"""
    PI = 3.1415926536
    def __init__(self,prot,compd,probeRadius=1.4):
        self.__init_direction_matrix()
        self.__probeRad = probeRadius
        self.__prot = prot
        self.__compd = compd
        pass

    def __init_direction_matrix(self):
        PI = SASA.PI
        #self.__directions = zeros((600,3))
        self.__directions = [[0.0]*3 for x in range(600)]
        interval = 0.05
        mpoints = int(1.0/interval+0.5)
        self.__ndirection = 0
        for i in range(mpoints+1):
            radialArc = 2.0 * math.sin((PI*float(i))/float(mpoints))
            longMax = int(radialArc/interval + 0.5)
            if longMax == 0: longMax = 1
            for j in range(longMax):
                fi = (2.0*PI*float(j))/longMax
                self.__directions[self.__ndirection][0] = math.sin((PI*float(i))/float(mpoints)) * math.cos(fi)
                self.__directions[self.__ndirection][1] = math.sin((PI*float(i))/float(mpoints)) * math.sin(fi)
                self.__directions[self.__ndirection][2] = math.cos((PI*float(i))/float(mpoints))
                self.__ndirection += 1

    def get_area(self):
        """ compute the area"""
        prot = self.__prot
        compd = self.__compd
        PI = SASA.PI
        mark_dict = {}
        for ligAtom in compd.atoms:
            ligAtom_radius = ligAtom.getUnitedRadius()
            neighbors = []
            for aAtom in compd.atoms:
                cutoff2 = (ligAtom_radius+aAtom.getUnitedRadius()+2.0*self.__probeRad)**2
                dist2 = atomic_distance2(ligAtom,aAtom)
                if (dist2 < cutoff2) and (dist2 > 0.001):
                   neighbors.append(aAtom)
            for aAtom in prot.atoms:
                cutoff2 = (ligAtom_radius+aAtom.getUnitedRadius()+2.0*self.__probeRad)**2
                dist2 = atomic_distance2(ligAtom,aAtom)
                if dist2 < cutoff2:
                   neighbors.append(aAtom)
           
            lx,ly,lz = ligAtom.getXYZ() 
            marks = 0
            for ii in range(self.__ndirection):
                xx = (ligAtom_radius+self.__probeRad)*self.__directions[ii][0]
                yy = (ligAtom_radius+self.__probeRad)*self.__directions[ii][1]
                zz = (ligAtom_radius+self.__probeRad)*self.__directions[ii][2]
                freeMark = True
                for aAtom in neighbors:
                    nx,ny,nz = aAtom.getXYZ()
                    nnx = nx - lx
                    nny = ny - ly
                    nnz = nz - lz
                    dist2 = (nnx-xx)*(nnx-xx) + (nny-yy)*(nny-yy) + (nnz-zz)*(nnz-zz)
                    cmp2 = (aAtom.getUnitedRadius()+self.__probeRad)**2
                    if dist2 < cmp2:
                        freeMark = False
                        break
                if freeMark: marks += 1

            atomic_area = float(marks)/float(self.__ndirection)*(4.0*PI*(ligAtom_radius+self.__probeRad)**2)
       
            atype = ligAtom.type
            if mark_dict.has_key(atype):
                mark_dict[atype] += atomic_area
            else:
                mark_dict[atype]  = atomic_area
        return mark_dict
        
class PLdescriptors(object):
    """Compute knowledge-based features using protein-ligand complexes for SVM model development"""
    # KB_stat_76.dat'
    # KB_stat_100.dat'
    # KB_stat_120.dat'
    # KB_stat_224.dat'

    __KB_data_name = 'KB_stat_76.dat'
    interface_cutoff = 9.5

    def __init__(self,protMol2,ligMol2,ligLabel):

        #read the Knowledge-based potentials
        self.__read_potentials()

        #read the protein mol2 file
        proteins = Molecules()
        proteins.readProtFromMol2(protMol2)
        prot = proteins.molecules[0]
        

        #read the ligand mol2 file
        self.__ligs = Molecules()
        self.__ligs.readCompdFromMol2(ligMol2)

        #assign the compound label
        self.__label = ligLabel

        self.__prot = self.__get_pocket(prot)

    def run(self,outputName='Protein-Compound.KB.features.out',protocol=1):
        """run the designed protocol"""

        #format the label
        labl = float(self.__label)
        if math.fabs(labl-1.000) < 0.0001:
           labl = '1'
        elif math.fabs(labl+1.000) < 0.0001:
           labl = '-1'
        elif math.fabs(labl) < 0.0001:
           labl = '0'
        else:
           labl = str('%.3f'%(labl))

        #open up the output file
        try:
            outputFile = open(outputName,'w')
        except IOError:
            print >> sys.stderr,'Not able to open ',outputName,'for its output'
            raise SystemExit

        prot = self.__prot
        #do the calculations
        if protocol == 1:
            for compd in self.__ligs.molecules: 
                results = self.__compute_pairwise_features(prot,compd,labl)
                

                #output the results for one molecule
                outputFile.write('%s '%(labl,))
                index = 0
                for feat in results:
                    if (feat > 0.0001) or (feat < -0.0001):
                       outputFile.write('%i:%.3f '%(index+1,feat,))
                    index += 1
                outputFile.write(' #%s:%s\n'%(compd.molID,prot.molID))
        elif protocol == 2:
            for compd in self.__ligs.molecules: 
                results = self.__compute_pairwise_features(prot,compd,labl)
                props = self.__compute_compd_features(compd)
                sasa = SASA.SASA(prot,compd)
                sasa_dict =  sasa.get_area()

                #output the results for one molecule
                outputFile.write('%s '%(labl,))
                index = 0
                for feat in results:
                    if (feat > 0.0001) or (feat < -0.0001):
                       outputFile.write('%i:%.3f '%(index+1,feat,))
                    index += 1
                for feat in props:
                    if (feat > 0.0001) or (feat < -0.0001):
                       outputFile.write('%i:%.3f '%(index+1,feat,))
                    index += 1

                for ii in range(20):
                    foundType = False
                    for atmType in self.__ligAtmTypes.keys():
                        if self.__ligAtmTypes[atmType] == ii:
                            foundType = True
                            break
                    if foundType:
                        try:
                            area = sasa_dict[atmType]
                        except KeyError:
                            area = 0.0
                        if area > 0.001:
                           outputFile.write('%i:%.3f '%(index+1,area,))
                        index += 1    
                outputFile.write(' #%s:%s\n'%(compd.molID,prot.molID))
        else:
            pass
        outputFile.close()

    def __compute_pairwise_features(self,prot,compd,labl):
        """compute the protein-ligand potentials"""
        cut_off2 = PLdescriptors.interface_cutoff**2.0
        #results = zeros(self.__num_interact)
        results = [0.0]*self.__num_interact
        for atomProt in prot.atoms:
            pAtp = atomProt.type
            try: 
                pIndex = self.__protAtmTypes[pAtp]
            except KeyError:
                continue
            px,py,pz = atomProt.getXYZ()
            for atomCompd in compd.atoms:
                lAtp = atomCompd.type
                try: 
                    lIndex = self.__ligAtmTypes[lAtp]
                except KeyError:
                    continue
                lx,ly,lz = atomCompd.getXYZ()
                dist2 = (px-lx)*(px-lx) + (py-ly)*(py-ly) + (pz-lz)*(pz-lz)
                if dist2 > cut_off2: continue
                rr = math.sqrt(dist2)
                if rr < 2.0:
                    m = 1
                elif rr < 8.0:
                    m = int((rr - 2.0+0.5)*2.0)+1
                else:
                    m = int(rr-8.0)+14
                poten = self.__pot[pIndex][lIndex][m-1]
                results[int(self.__featMap[pIndex][lIndex])-1] += poten
        return results

    def __compute_compd_features(self,compd):
        prop = []
        sum = 0.0
        for atom in compd.atoms:
            sum += atom.getMass()
        prop.append(sum)
        return prop

    def __read_potentials(self):
        """read the knowledge-based potentials input file"""
        if os.path.exists('.'+os.sep+PLdescriptors.__KB_data_name):
            KBdataName = '.'+os.sep+PLdescriptors.__KB_data_name
        else:
            path = os.getenv('SVMSP_PARAM')
            if not path:
                print >> sys.stderr,'Can not find the Knowledge potential data file'
                raise SystemExit
            else:
                KBdataName = path+os.sep+PLdescriptors.__KB_data_name
        try:
            potFile = open(KBdataName)
        except IOError:
            print >> sys.stderr,'Can not read the Knowledge potential data file'
            raise SystemExit

        num1Max,num2Max,num3Max,num4Max = 0,0,0,0
        self.__protAtmTypes = {}
        self.__ligAtmTypes = {}
        for line in potFile:
            if line.startswith('#'): continue
            pTyp,lTyp,pot,dist,pIndex,lIndex,intIndex = line.split()
            pIndex = int(pIndex)
            lIndex = int(lIndex)
            dist = int(dist)
            self.__protAtmTypes[pTyp] = pIndex
            self.__ligAtmTypes[lTyp] = lIndex
            intIndex = int(intIndex)
            if pIndex > num1Max: num1Max = pIndex
            if lIndex > num2Max: num2Max = lIndex
            if dist   > num3Max: num3Max = dist
            if intIndex >  num4Max: num4Max = intIndex
        self.__num_interact = intIndex
        #self.__pot = zeros([num1Max+1,num2Max+1,num3Max])
        self.__pot = [[[0.0]*int(num3Max) for x in range(int(num2Max+1))] for x in range(int(num1Max+1))]
        #self.__featMap = zeros([num1Max+1,num2Max+1])
        self.__featMap = [[0.0]*int(num2Max+1) for x in range(int(num1Max+1))]

        potFile.seek(0)
        for line in potFile:
            if line.startswith('#'): continue
            pTyp,lTyp,pot,dist,pIndex,lIndex,intIndex = line.split()
            pIndex = int(pIndex)
            lIndex = int(lIndex)
            dist = int(dist)
            pot = float(pot)
            intIndex = int(intIndex)
            self.__pot[pIndex][lIndex][dist-1] = pot
            self.__featMap[pIndex][lIndex] = intIndex
        return
        
    def __get_pocket(self,prot):
        pocket = PTmolecule()
        pocket.molID = prot.molID
        cutoff = PLdescriptors.interface_cutoff + 6.5
        cutoff2 = cutoff * cutoff
        atomIDs = []
        for ii in range(5):
           try:
               compd = self.__ligs.molecules[ii]
           except IndexError:
               break
           for ligAtom in compd.atoms:
               if ligAtom.type == 'h': continue
               lx,ly,lz = ligAtom.getXYZ()
               for protAtom in prot.atoms:
                   if protAtom.type == 'H': continue
                   if protAtom.type == 'LP': continue
                   aID = protAtom.id
                   if aID in atomIDs:
                       continue
                   else:
                       px,py,pz = protAtom.getXYZ()
                       dist2 = atomic_distance2(ligAtom,protAtom)
                       if dist2 < cutoff2:
                          atomIDs.append(aID)
        for protAtom in prot.atoms:
            aID = protAtom.id
            if aID in atomIDs:
                pocket.addAtom(protAtom)
        return pocket

def main():
    if len(sys.argv) != 4 and len(sys.argv) != 5:
       print >> sys.stderr,'Usage: PLdescriptors.py <protein.mol2> <compounds.mol2> <label> [output Name]'
       raise SystemExit
    elif len(sys.argv) == 4:
       features = PLdescriptors(*sys.argv[1:])
       features.run(protocol=1)
    elif len(sys.argv) == 5:
       features = PLdescriptors(*sys.argv[1:4])
       features.run(sys.argv[4],protocol=1)

if __name__ == '__main__':
    main()
