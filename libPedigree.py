#! usr/bin/env python
# Copyright 2009 Harald Grove

import sys

class Pedigree(object):
    
    def __init__(self,infile=None):
        self.families = {} # Based on animals
        self.animals = {}
        self.sires = {}
        self.dams = {}
        self.familyOrder = [] # List to record input order of families
        self.founders = [] # All animals without sire&dam (both 0)
        try:
            self.fped = infile
            self.readPedigreeFile()
        except:
            pass
            
    def addAnimal(self,animal,dam,sire,families,sex='3'):
        """ A '0' for dam or sire indicates missing information """
        if not isinstance(families,list):
            families = [families]
        for family in families:
            if not animal in self.animals:
                self.animals[animal] = {'dam':dam,'sire':sire,'family':[family],'sex':sex}
            else:
                if family not in self.animals[animal]['family']: self.animals[animal]['family'].append(family)
                if self.animals[animal]['dam'] == '0': self.animals[animal]['dam'] = dam
                elif self.animals[animal]['dam'] == dam or dam == '0': pass
                else: 
                    sys.stderr.write('%s is recorded with 2 different dams\n' % (animal))
                    sys.exit(0)
                if self.animals[animal]['sire'] == '0': self.animals[animal]['sire'] = sire
                elif self.animals[animal]['sire'] == sire or sire == '0': pass
                else: 
                    sys.stderr.write('%s is recorded with 2 different sires\n' % (animal))
                    sys.exit(0)
            if family != '0' and not family in self.families:
                self.families[family] = [animal]
                self.familyOrder.append(family)
            elif family in self.families and animal not in self.families[family]:
                self.families[family].append(animal)
            if sire != '0' and not sire in self.sires:
                self.sires[sire] = [animal]
            elif sire in self.sires:
                self.sires[sire].append(animal)
            if dam != '0' and not dam in self.dams:
                self.dams[dam] = [animal]
            elif dam in self.dams:
                self.dams[dam].append(animal)
        
    def readPedigreeFile(self):
        """ Reads pedigree as: animalID sireID damID familyID [sexID] """
        for line in self.fped:
            l = line.strip().split()
            animal = l[0]
            #if len(l) > 1: dam = l[1]
            #else: dam = '0'
            #if len(l) > 2: sire = l[2]
            #else: sire = '0'
            if len(l) > 1: sire = l[1]
            else: sire = '0'
            if len(l) > 2: dam = l[2]
            else: dam = '0'
            if len(l) > 3: family = l[3]
            else: family = 'F0'
            if len(l) > 4: sex = l[4]
            else: sex = '3'
            self.addAnimal(animal, dam, sire,family,sex)
        for anim in self.animals:
            if self.animals[anim]['sex'] == '3':
                if anim in self.sires: self.animals[anim]['sex'] = '1'
                elif anim in self.dams: self.animals[anim]['sex'] = '0'
    
    def getAnimals(self):
        """ Returns a dict containing all animals
        """
        anims = []
        for animal in self.animals:
            anims.append(animal)
        return anims
            
    def getSire(self,animal):
        try:
            return self.animals[animal]['sire']
        except KeyError:
            return '0'
        
    def getDam(self,animal):
        try:
            return self.animals[animal]['dam']
        except KeyError:
            return '0'
        
    def getFamily(self,animal):
        """ Returns a list of families """
        return self.animals[animal]['family']
        
    def getSex(self,animal):
        try:
            return self.animals[animal]['sex']
        except KeyError:
            return '3'
        
    def getNumFamilies(self):
        return len(self.families)
        
    def getFamilies(self):
        return [f for f in self.familyOrder]
        
    def getFamilyMembers(self,family):
        return [m for m in self.families[family]]
        
    def getNumAnimals(self):
        return len(self.animals)

    def getOffspring(self,animal):
        if animal in self.sires:
            return [o for o in self.sires[animal]]
        elif animal in self.dams:
            return [o for o in self.dams[animal]]
        else:
            return []
            
    def getSires(self):
        return [s for s in self.sires]
        
    def getDams(self):
        return [d for d in self.dams]
    
    def getFoundingAnimals(self,family):
        """ Returns founding sires and dams as two lists, defined as having no parents in the pedigree
        """
        animals = self.getFamilyMembers(family)
        sires,dams = [],[]
        for animal in animals:
            if self.getSire(animal) == '0' and self.getSex(animal) == '1': sires.append(animal)
            elif self.getSire(animal) == '0' and self.getSex(animal) == '0': dams.append(animal)
        return [sires,dams]
    
    def findAnimal(self,animal):
        """ Returns True if animal is present in data
        """
        if animal in self.animals: return True
        return False
    
    def setSire(self,animal,sire):    
        self.animals[animal]['sire'] = sire
        
    def setDam(self,animal,dam):
        self.animals[animal]['dam'] = dam

    def setSex(self,animal,sex):
        self.animals[animal]['sex'] = sex

    def printPedigree(self):
        sep = '\t'
        out = ''
        for family in self.getFamilies():
            for anim in self.getFamilyMembers(family):
                if anim in self.animals:
                    out += anim+sep+self.getSire(anim)+sep+self.getDam(anim)+sep+family+sep+self.getSex(anim)+'\n'
        return out

    def gather(self,animal,family):
        """ Collects all animal that are connected in some way
        """
        if animal in self.printedAnimals or animal == '0': return # Already been here, or reached the end
        self.printedAnimals[animal] = 1
        sire,dam = self.getSire(animal),self.getDam(animal)
        children = self.getOffspring(animal)
        self.gather(sire,family)
        self.gather(dam,family)
        for child in children:
            self.gather(child,family)
        print animal+'\t'+sire++'\t'+dam+'\t'+family+'\t'+self.getSex(animal)
        

    def printFullped(self):
        """ Prints a pedigree gathering all connected animals into the same family
        """
        self.printedAnimals = {}
        families = self.getFamilies()
        for family in families:
            members = self.getFamilyMembers(family)
            out = self.gather(members[0],family)
        
    def __str__(self):
        sep = '\t'
        out = ''
        for family in self.families:
            members = self.getFamilyMembers(family)
            for animal in members:
                out += animal+sep+self.getSire(animal)+sep+self.getDam(animal)+sep+family+sep+self.getSex(animal)+'\n'
        return out
        
def main():
    infile = open(sys.argv[1],'r')
    ped = Pedigree(infile)
    ped.printFullped()

if __name__ == '__main__':
    main()
