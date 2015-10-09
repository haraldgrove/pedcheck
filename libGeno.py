#!/usr/bin/env python

# libGeno v1.0
# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys
import libPed
import libMark
import gzip

class Translate(object):
    
    def __init__(self):
        self.t = {'A':'1','C':'2','G':'3','T':'4','DEL':'5','D':'5','1':'1','2':'2','3':'3','4':'4','9':'9','0':'0'}
    
    def __getitem__(self,key):
        try: return self.t[key]
        except: return '0'

class Geno(object):
    def __init__(self, infile,markfile='0',pedfile='0',mode=''):
        """ if markfile and pedfile are not valid files, markobj and pedobj are generated locally and can be exported
        """
        zipfile = False
        if infile[-3:] == '.gz': zipfile = True
        self.genotypes = {} # Genotypes are stored by animalID
        if zipfile: self.fgeno = gzip.open(infile,'r')
        else: self.fgeno = open(infile,'r')
        self.mode = mode
        self.mafinfo = {}
        try:
            fmark = open(markfile,'r')
            self.mark = libMark.Mark(markfile)
            line = self.fgeno.next().strip().split()
            self.localMark = False
            if line[0] != '#':
                if zipfile: self.fgeno = gzip.open(infile,'r')
                else: self.fgeno = open(infile,'r') # No markerline
        except IOError:
            line = self.fgeno.next().strip().split()
            if line[0] != '#':
                if zipfile: self.fgeno = gzip.open(infile,'r')
                else: self.fgeno = open(infile,'r')
                line = ['#']+['M'+str(i) for i in xrange( (len(line)-3)/2 )]
            self.mark = libMark.Mark(line[1:])
            self.localMark = True
        try:
            fped = open(pedfile,'r')
            self.ped = libPed.Ped(pedfile)
            self.localPed = False
        except IOError:
            self.ped = libPed.Ped()
            self.localPed = True
        #self.fgeno = open(infile,'r') # another way to reset the file iterator?
        self.importFile(mode)
        self.fgeno.close()

    def __getitem__(self, key):
        animal,marker = key[0],self.mark[str(key[1])][1] # Input is the markerID as found in the markerobject
        return self.genotypes[animal][marker*2:marker*2+2]

    def __setitem__(self,key,value):
        #if value[0] not in ['0','1','2','3','4','5','A','C','G','T'] or value[1] not in ['0','1','2','3','4','5','A','C','G','T']:
        #    print "Setting illegal genotype",key,value
        #    sys.exit(1)
        try:
            animal,marker = key[0],self.mark[key[1]][1] # Input is the markerID as found in the markerobject
            self.genotypes[animal][marker*2:marker*2+2] = value
            #if value[0] != '0' and value[1] == '0' and self.mode != 'hap':
            #    print key,value
            #    sys.exit(1)
        except (KeyError,ValueError):
            print "Not recognized key",key
            sys.exit(1)
            #id1,id2 = key
            #if id1 == 'genotype' and id2 in self.genotypes: self.genotypes[id2] = value

    def importFile(self,mode):
        """ Assumes unique pedigree in input, or that multiple entries are identical """
        #if not self.localMark:
        #    self.markIndex = [[self.mark[str(mark)][1]*2,self.mark[str(mark)][1]*2+1] for mark in self.mark['markers']]
        count = 0
        for line in self.fgeno:
            if line.startswith('#'):
                count += 1
                continue
            try: animal,sire,dam,geno = line.strip().split(None,3)
            except ValueError:
                sys.stderr.write('Error: Wrong number of info-columns or missing genotypes in line %d: found %d but expected 4\n'\
                                  % (count,len(line.strip().split())))
                sys.exit(1)
            if not self.localPed and animal not in self.ped.getAnimals():
                count += 1
                continue # with external pedigree, skip unlisted animals
            if self.localPed: self.ped.addAnimal(animal,dam,sire,'F0','3')
            if mode != '0': self.addGenotype(animal,geno.split())
            else: self.addGenotype(animal,['0','0']*len(self.mark)) # Used when generating a blank genotype file
            if mode == 'hap':
                for i in xrange(0,len(self.genotypes[animal]),2):
                    a1,a2 = self.genotypes[animal][i:i+2]
                    if a1 != a2: self.genotypes[animal][i:i+2] = ['0','0']
            self.ped.updateSex()
            count += 1
        

    def addGenotype(self,animal,genos):
        """ Used for adding new animals to the data, genos is a list of genotypes
        """
        corr1 = corr2 = corr3 = 0
        if len(genos) == len(self.mark):
            temp = []
            for g in genos:
                temp.append(g[0])
                temp.append(g[1])
            genos = temp
        if animal in self.genotypes: # used when two animals with the same name are present in the data
            for i in xrange(len(self.mark)):
                oa1,oa2 = self.genotypes[animal][2*i:2*i+2]
                na1,na2 = genos[2*i:2*i+2]
                if [na1,na2] == [oa1,oa2] or [na2,na1] == [oa1,oa2]: continue
                if oa1 == oa2 == '0':
                    self.genotypes[animal][2*i:2*i+2] = [na1,na2]
                    corr1 += 1
                elif na1 == na2 == '0': corr2 += 1
                else:
                    self.genotypes[animal][2*i:2*i+2] = ['0','0'] # If it's not possible to fill in later, then don't guess here
                    corr3 += 1
        else:
            self.genotypes[animal] = genos
            #if self.localMark: self.genotypes[animal] = genos
            #else:
            #    self.genotypes[animal] = [genos[mark] for markblock in self.markIndex for mark in markblock]
        if (corr1+corr2+corr3) > 0: print "%d %d %d in %s" % (corr1,corr2,corr3, animal)

    def updateMarkers(self):
        for mark in self.mark.getMarkers():
            a1 = a2 = '0'
            for animal in self.genotypes:
                if a1 != '0' and a2 != '0': break
                try: allele = self.__getitem__([animal,mark])
                except (KeyError,IndexError): continue
                if '0' in allele: continue
                try:
                    if a1 == '0': a1 = allele[0]
                except IndexError:
                    print animal,mark,len(self.genotypes[animal])
                    print len(self.mark.getMarkers())
                    sys.exit(1)
                if a2 == '0' and a1 != '0': 
                    if a1 != allele[0]: a2 = allele[0]
                    elif a1 != allele[1]: a2 = allele[1]
            if a2 == '0': a2 = a1 # Set marker to homozygous if only one allele was found
            if a1 not in ['0','1','2','3','4','5','9'] or a2 not in ['0','1','2','3','4','5','9']:
                pass
                #print "<updateMarker>",a1,a2,mark,animal
            self.mark['alleles',mark] = sorted([a1,a2])
    
    def getMarkers(self):        
        """if key == 'markers':
                self.updateMarkers()
                return self.mark """
        self.updateMarkers()
        return self.mark
                
    def getPedigree(self):
        """
            elif key == 'pedigree': return self.ped
        """
        return self.ped
        
    def getPed(self):
        return self.ped
    
    def getAnimals(self):
        """
        """
        return [anim for anim in self.genotypes]
    
    def getGenotype(self,animal):
        """
        """
        return self.genotypes[animal]
        
    def getHaplotype(self,animal):
        """
        """
        return self.genotypes[animal][0::2][:],self.genotypes[animal][1::2][:]
        
    def getMAlleles(self, marker):
        """
        """
        try: return self.mark.getMAlleles(marker)
        except:
            self.updateMarkers()
            return self.mark.getMAlleles(marker)
    
    def getMAF(self,marker):
        try: return self.mafinfo[marker]
        except: return self.getMAF_(marker)
    
    def getMAF_(self,marker):
        """ Returns [minor allele,major allele, maf]
        """
        m1,m2,c1,tot = '0','0',0,0
        for animal in self.ped.getAnimals():
            a1,a2 = self.__getitem__([animal,marker])
            if '0' in a1+a2: continue
            tot += 2
            if m1 == '0': m1 = a1
            if m2 == '0' and m1 != a1: m2 = a1
            elif m2 == '0' and m1 != a2: m2 = a2
            if a1 == m1: c1 += 1
            if a2 == m1: c1 += 1
        if tot == 0: return '0','0',0
        if c1 > tot/2.0:
            self.mafinfo[marker] = [m2,m1,1-c1/float(tot)]
            return m2,m1,1-c1/float(tot)
        else:
            self.mafinfo[marker] = [m1,m2,c1/float(tot)]
            return m1,m2,c1/float(tot)
            
    def setGenotype(self,animal,geno):
        if animal in self.genotypes: self.genotypes[animal] = geno
            
    def __str__(self):
        sep = '\t'
        out = '#'+sep+(sep+sep).join(self.mark.getMarkers())+'\n'
        for animal in self.ped.getAnimals():
            out += animal+sep+self.ped.getSire(animal)+sep+self.ped.getDam(animal)+sep+sep.join(self.genotypes[animal])+'\n'
        return out
        
    def fastPrint(self,fout):
        sep = '\t'
        fout.write('#'+sep+(sep+sep).join(self.mark.getMarkers())+'\n')
        for animal in self.ped.getAnimals():
            fout.write(animal+sep+self.ped.getSire(animal)+sep+self.ped.getDam(animal)+sep+sep.join(self.genotypes[animal])+'\n')
            
        
    def printData(self, fout, ped = [],markers = []):
        """ Use when modifications have been done to the markerlist
        """
        sep = '\t'
        if ped == []: ped = self.ped
        if markers == []: markers = self.mark
        fout.write('#'+sep+sep.join(markers.getMarkers())+'\n')
        for animal in ped.getAnimals():
            fout.write(animal+sep+ped.getSire(animal)+sep+ped.getDam(animal))
            fout.write(sep+sep.join([sep.join(self.__getitem__([animal,mark])) for mark in markers.getMarkers()])+'\n')
        
    #******************************************************************
    # Old functions to keep compatibility with other functions
    #******************************************************************
        
    def getOldPed(self):
        """ Creates a copy of the pedigree object
        """
        import libPedigree
        outped = libPedigree.Pedigree()
        for animal in self.ped.getAnimals():
            family,sire,dam,sex = self.ped.getFamily(animal),self.ped.getSire(animal),self.ped.getDam(animal),self.ped.getSex(animal)
            for fam in family:
                outped.addAnimal(animal,dam,sire,fam,sex)
        return outped
        
    def getOldMarkers(self,input = None):
        self.updateMarkers()
        ut = []
        if input:
            chr,nr,pos = '99',1,0
            for marker in input:
                try: m1,m2 = self.getMAlleles(marker)
                except KeyError: m1,m2 = '0','0'
                ut.append([marker,str(pos),m1,m2,chr])
                pos += 1
            return ut
        count = 1
        for marker in self.mark.getMarkers():
            chr,nr,id,pos = self.mark[marker]
            m1,m2 = self.getMAlleles(marker)
            ut.append([marker,str(pos),m1,m2,chr])
        return ut
        
    def getGenotype(self,animal):
        return self.genotypes[animal]
        
    def getAlleles(self,animal,marker,opt = 'lett'):
        #m1,m2 = self.__getitem__('alleles',marker)
        a1,a2 = self.__getitem__([animal,marker])
        if opt == 'lett':
            if a1 == '1': a1 = 'A'
            elif a1 == '2': a1 = 'C'
            elif a1 == '3': a1 = 'G'
            elif a1 == '4': a1 = 'T'
            if a2 == '1': a2 = 'A'
            elif a2 == '2': a2 = 'C'
            elif a2 == '3': a2 = 'G'
            elif a2 == '4': a2 = 'T'
        return a1,a2
        
    def getMarkerGenotype(self,marker):
        m1,m2 = self.getMAlleles(marker)
        return m1,m2
        
    def getAssays(self):
        return self.mark.getMarkers()
        
    def getSNPNumbers(self):
        return len(self.mark)
        
    #****************************************************************************************
    # Functions reporting from the object
    #****************************************************************************************
    
    def findAlleleParent(self,animal,marker):
        """ Returns list of alleles that must be inherited from parents, empty if none
        """
        pA = []
        try: s1,s2 = self.__getitem__([self.getSire(animal),marker])
        except KeyError: s1 = s2 = '0'
        try: d1,d2 =self.__getitem__([self.getDam(animal),marker])
        except KeyError: d1 = d2 = '0'
        if s1 == s2 and s1 != '0': pA.append(s1)
        if d1 == d2 and d1 != '0': dA.append(d1)
        return pA
        
    def findAlleleChildren(self,parent,marker,limit = 1000):
        """ Returns list of alleles that must be given to offspring 
            Assumes haplotypes
        """
        offspring = self.ped.getOffspring(parent)
        if len(offspring) == 0: return []
        m1,m2 = self.getMAlleles(marker)
        count = []
        for off in offspring:
            try: o1,o2 = self.__getitem__([off,marker])
            except KeyError: continue
            if o1 == '0': continue
            if o1 not in count: count.append(o1)
            if len(count) > 1: return count
        return []
        
    def wrongMendel(self, animal,sire,dam,marker):
        """ Returns  1,2,3,4 if there's a discord agains sire,dam,either,both
            returns -1 if there's no information to base a decision on
            returns 0 otherwise
        """
        try: a1,a2 = self.__getitem__([animal,marker])
        except KeyError: return -1
        try: s1,s2 = self.__getitem__([sire,marker])
        except KeyError: s1 = s2 = '0'
        try: d1,d2 = self.__getitem__([dam,marker])
        except KeyError: d1 = d2 = '0'
        if s1+d1 == '00': return -2 # Both parents are ungenotyped or missing
        if a1 == a2 == '0': return -1
        if d1 == '0': # only sire has genotype
            if a1 == a2 and s1 == s2 and a1 != s1: return 1
            return 0
        if s1 == '0': # only dam has genotype
            if a1 == a2 and d1 == d2 and a1 != d2: return 2
            return 0
        # Both sire and dam have genotypes
        if a1 != a2:
            if d1 == d2 == s1 == s2: return 3
            return 0
        # animal is hom: a1 = a2
        count = 0
        if a1 in s1+s2: count += 1
        if a1 in d1+d2: count += 2
        if count == 3: return 0
        if count == 2: return 1
        if count == 1: return 2
        if count == 0: return 4
        
    def wrongMendelS(self, animal,sire,marker):
        """ Returns  1 if there's a discord against sire
            returns -1 if there's no information to base a decision on
            returns 0 otherwise
        """
        try: a1,a2 = self.__getitem__([animal,marker])
        except KeyError: return -1
        try: s1,s2 = self.__getitem__([sire,marker])
        except KeyError: return -2
        if s1 == '0': return -2 # Both parents are ungenotyped or missing
        if a1 == '0':  return -1
        if a1 == a2 and s1 == s2 and a1 != s1: return 1
        return 0

    def kinship(self, animal,sire,marker):
        """ Return kinship coefficient between animal and sire at marker
            Returns -2 if there is not enough information to calculate
        """
        def transkin(s1,s2,a,b):
            if s1 != s2: return 0.5
            if s1 == a: return 0
            if s1 == b: return 1
            return 'NA'
        
        try: a1,a2 = self.__getitem__([animal,marker])
        except KeyError: return -2
        try: s1,s2 = self.__getitem__([sire,marker])
        except KeyError: return -2
        m1,m2,maf = self.getMAF(marker)
        if maf == 0: return -2
        ga = transkin(a1,a2,m1,m2)
        gs = transkin(s1,s2,m1,m2)
        if ga == 'NA' or gs == 'NA': return -2
        pk = 1-maf
        value = (ga-pk)*(gs-pk) / (pk*(1-pk))
        return value

    def findLegalOffspring(self, offspring, marker, limit=5):
        """ Returns '0','0' if it's not a unique answer """
        if len(offspring) == 0: return '0','0'
        m1,m2 = self.getMAlleles(marker)
        count = [0,0,0]
        for off in offspring:
            try: o1,o2 = self.__getitem__([off,marker])
            except KeyError: continue
            if o1 == '0': continue
            if o1 != o2: count[1] += 1
            elif o1 == m1: count[0] += 1
            else: count[2] += 1
            if count[0] > 0 and count[2] > 0: return m1,m2
        if count[1] > 0 or count[0] + count[2] == 0: return '0','0'
        if len(offspring) > limit and count[0] == len(offspring): return m1,m1
        if len(offspring) > limit and count[2] == len(offspring): return m2,m2
        return '0','0'
                                                                                                                                                                                                                                                                                                                                                                                    
    def findLegalAlleles(self,offspring,sire,dam,marker):
        """ Only returns 0 or 2 alleles """
        try: s1,s2 = self.__getitem__([sire,marker])
        except KeyError: s1 = s2 = '0'
        try: d1,d2 = self.__getitem__([dam,marker])
        except KeyError: d1 = d2 = '0'
        o1,o2 = self.findLegalOffspring(offspring,marker)
        if o1 == '0':
            if s1 == s2 and s1 != '0' and d1 == d2 and d1 != '0': return s1,d1
            return '0','0'
        if o1 != o2:
            if s1 == s2 == d1 == d2 and s1 != '0': pass # Disagreement between parents and children, children takes precedens
            return o1,o2
        count = 0
        if o1 not in s1+s2: count += 1
        if o1 not in d1+d2: count += 2
        if count == 0: return o1,o2
        return '0','0'
        
    def checkMendel(self,ferr,animals=[],markers=[]):
        # Running this function on a haplotype file is not recommended!!
        if animals == []: animals = self.ped
        if markers == []: markers = self.mark
        report = {}
        animrep = {}
        ferr.write('#Changed_genotypes\n')
        for mark in markers.getMarkers():
            report[mark] = 0
            m1,m2 = self.getMAlleles(mark)
            errorMendel = {}
            fixMendel = {}
            findAllele = []
            # Record all parent-offspring discords
            for anim in animals.getAnimals():
                try: a = self.__getitem__([anim,mark])
                except KeyError: continue
                if a == ['0','0']:
                    findAllele.append(anim)
                    continue
                sire,dam = animals.getSire(anim),animals.getDam(anim)
                if sire == '0' and dam == '0': continue
                if dam == '0': result = self.wrongMendelS(anim,sire,mark)
                else: result = self.wrongMendel(anim,sire,dam,mark)
                if result == 0: continue
                if result in [1,3,4]:
                    errorMendel[sire] = errorMendel.get(sire,0)+1
                    fixMendel[sire] = fixMendel.get(sire,[]) + [anim]
                if result in [2,3,4]:
                    errorMendel[dam] = errorMendel.get(dam,0)+1
                    fixMendel[dam] = fixMendel.get(dam,[]) + [anim]
            # Delete genotypes of parent if more than 2 offsprings disagree, otherwise delete genotypes for offsprings
            for parent in errorMendel:
                if errorMendel[parent] > 2:
                    self.__setitem__([parent,mark],['0','0'])
                    findAllele.append(parent)
                    report[mark] += 1
                    animrep[parent] = animrep.get(parent,0)+1
                    try: ferr.write('P'+'\t'+parent+'\t'+mark+'\t'+'0'+'\n')
                    except: pass
                else:
                    for fixAnim in fixMendel[parent]:
                        self.__setitem__([fixAnim,mark],['0','0'])
                        findAllele.append(fixAnim)
                        report[mark] += 1
                        animrep[fixAnim] = animrep.get(fixAnim,0) + 1
                        try: ferr.write('O'+'\t'+fixAnim+'\t'+mark+'\t'+'0'+'\n')
                        except: pass
            report[mark] = report[mark],report[mark]/float(len(animals))
            # Check all the deleted or missing genotypes if they can be inferred from parents/children
            for animal in findAllele:
                sire,dam = animals.getSire(animal),animals.getDam(animal)
                offspring = animals.getOffspring(animal)
                xa = self.findLegalAlleles(offspring,sire,dam,mark)
                if xa[0] != '0' and xa[1] == '0':
                    print animal,mark,xa
                    sys.exit(1)
                if xa == ['0','0']: continue
                self.__setitem__([animal,mark],xa)
                try: ferr.write('N'+'\t'+animal+'\t'+mark+'\t'+xa+'\n')
                except: pass # Wrong exception
            # Print all animals with a discord > 0
        ferr.write('#Animal_discords\n')
        for anim in animals.getAnimals():
            if anim in animrep and animrep[anim] > 0: ferr.write(anim+'\t'+str(animrep[anim])+'\n')
        return report
        
    def fixImpute(self,ferr,animals=[],markers=[]):
        """ Corrects positions where an offspring has gotten the opposite homozygote alleles to the parent
        """
        if animals == []: animals = self.ped
        if markers == []: markers = self.mark
        for mark in markers.getMarkers():
            m1,m2 = self.getMAlleles(mark)
            for anim in animals.getAnimals():
                try: a = self.__getitem__([anim,mark])
                except KeyError: continue
                sire,dam = animals.getSire(anim),animals.getDam(anim)
                if sire == '0' and dam == '0': continue
                if dam == '0': result = self.wrongMendelS(anim,sire,mark)
                else: result = self.wrongMendel(anim,sire,dam,mark)
                if result == 0: continue
                self.__setitem__([anim,mark],[m1,m2])
                ferr.write(anim+'\t'+mark+'\t'+sire+'\t'+dam+'\n')
        
    def checkPed(self,animal,otherAnimal=''):
        """ Checks animal against either the given otherAnimal or against recorded sire and dam
        """
        if otherAnimal == '':
            try:
                sire = self.ped.getSire(animal)
                dam = self.ped.getDam(animal)
            except KeyError:
                sire = '0'
                dam = '0'
            if sire == '0' and dam == '0': return -1,-1
        else:
            sire = otherAnimal
            dam = '0'
        sc,dc,tot = 0,0,0
        for marker in self.mark.getMarkers():
            result = self.wrongMendel(animal,sire,dam,marker)
            if result == 0: tot += 1
            elif result > 0:
                tot += 1
                if result in [1,3,4]: sc += 1
                if result in [2,3,4]: dc += 1
        if tot == 0: return -1,-1
        return sc*100/float(tot),100*dc/float(tot)
        
    def checkPedS(self,animal,sire):
        """ Checks animal against the given 'parent'
        """
        sc,tot = 0,0
        if sire in ['','0']: return -1
        for marker in self.mark.getMarkers():
            result = self.wrongMendelS(animal,sire,marker)
            if result == 0: tot += 1
            elif result > 0:
                tot += 1
                sc += 1
        if tot == 0: return -1
        return sc*100/float(tot)
        
    def checkIBS(self,animal,sire):
        """ Calculates ibs between animal and sire
        """
        kinc = []
        if sire in ['','0']: return 'NA'
        for marker in self.mark.getMarkers():
            kin = self.kinship(animal,sire,marker)
            if kin > -2: kinc.append(kin)
        if len(kinc) == 0: return 'NA'
        return sum(kinc)/len(kinc)
        
    def anyDams(self):
        for animal in self.ped.getAnimals():
            if self.ped.getDam(animal) != '0': return True
        return False
        
    def switchAlleles(self):
        """ Changes markers genotyped as CT/GT to GA/CA
        """
        for marker in self.mark.getMarkers():
            for animal in self.ped.getAnimals():
                a1,a2 = self.__getitem__([animal,marker])
                count = 0
                if a1 in ['T','4'] and count in [0,2]: count += 1
                elif a2 in ['T','4'] and count in [0,2]: count += 1
                if a1 in ['C','G','2','3'] and count in [0,1]: count += 2
                elif a2 in ['C','G','2','3'] and count in [0,1]: count += 2
                if count == 3:
                    self.switchStrand(marker)
                    break
                        
    def switchStrand(self,mark):
        """ changes from one strand to the other
        """
        for member in self.ped.getAnimals():
            try: [a1,a2] = self.__getitem__([member,mark])
            except KeyError: [a1,a2] = ['0','0']
            if a1 == '4': a1 = '1'
            elif a1 == '2': a1 = '3'
            elif a1 == '3': a1 = '2'
            elif a1 == '1': a1 = '4'
            if a2 == '4': a2 = '1'
            elif a2 == '2': a2 = '3'
            elif a2 == '3': a2 = '2'
            elif a2 == '1': a2 = '4'
            self.__setitem__([member,mark],[a1,a2])
        
    def reportStats(self,animals=[],markers=[]):
        if animals == []: animals = self.ped
        if markers == []: markers = self.mark
        report = {}
        for mark in markers.getMarkers():
            m1,m2 = markers.getMAlleles(mark)
            MAF,propgeno = 0,0
            for anim in animals.getAnimals():
                try: a1,a2 = self.__getitem__([anim,mark])
                except KeyError: continue # Missing allele
                if a1 == '0': continue
                if a1 == m2: MAF += 1
                if a2 == m2: MAF += 1
                propgeno += 1
            try: MAF = min(MAF,propgeno*2-MAF)/(propgeno*2.0)
            except ZeroDivisionError: MAF = 0
            report[mark] = [MAF,float(propgeno)/len(animals)]
        return report
        
    
        
    def removeGenotypes(self,fout,p):
        """ Randomly selects p% of the known genotypes and
            writes out a list with details
        """
        import random
        sep = '\t'
        p = p/100.0
        for member in self.ped.getAnimals():
            for marker in self.mark.getMarkers():
                try: a1,a2 = self.__getitem__([member,marker])
                except KeyError: a1,a2 = '0','0'
                if a1 == '0' and a2 == '0': continue
                r = random.random() # [0,1]
                if r < p: fout.write(marker+sep+member+sep+'0'+sep+'0'+sep+a1+sep+a2+'\n')
    
    def findPhasing(self,a,p1,p2,symbol='io'):
        """ returns 'i','o','-' if a is equal to p1 or p2 or not informative, respectively
        """
        if p1 not in ['0','9'] and p2 not in ['0','9'] and a not in [p1,p2,'0','9']: return 'x' # a is not a legal offspring of p1+p2
        if p1 in ['0','9'] or p2 in ['0','9'] or p1 == p2 or a in ['0','9']: return '-'
        if a == p1: return symbol[0]
        if a == p2: return symbol[1]
    
    def getPhase(self,animal,sire,dam):
        """ Returns i/o phases towards sire and dam
        """
        phases = [[],[]]
        for mark in self.mark.getMarkers():
            [aI,aO] = self.__getitem__([animal,mark])
            try: [sI,sO] = self.__getitem__([sire,mark])
            except: sI,sO = '0','0'
            try: [dI,dO] = self.__getitem__([dam,mark])
            except: dI,dO = '0','0'
            if [sI,sO] == ['0','0']: phases[0].append('-')
            else:
                try: gs1,gs2 = self.__getitem__([self.ped.getSire(sire),mark])
                except: gs1,gs2 = '0','0'
                try: gd1,gd2 = self.__getitem__([self.ped.getDam(sire),mark])
                except: gd1,gd2 = '0','0'
                if (gs1 == gs2 and gs1 != '0') or (gd1 == gd2 and gd1 != '0'): phases[0].append(self.findPhasing(aI,sI,sO,'IO'))
                else: phases[0].append(self.findPhasing(aI,sI,sO))
            if [dI,dO] == ['0','0']: phases[1].append('-')
            else:
                try: gs1,gs2 = self.__getitem__([self.ped.getSire(dam),mark])
                except: gs1,gs2 = '0','0'
                try: gd1,gd2 = self.__getitem__([self.ped.getDam(dam),mark])
                except: gd1,gd2 = '0','0'
                if (gs1 == gs2 and gs1 != '0') or (gd1 == gd2 and gd1 != '0'): phases[1].append(self.findPhasing(aO,dI,dO,'IO'))
                else: phases[1].append(self.findPhasing(aO,dI,dO))
        return phases
        
    def getSafePhase(self,animal,sire,gsire):
        """ Returns I/O phases towards grandsire and granddam
        """
        phases = []
        for mark in self.mark.getMarkers():
            [aI,aO] = self.__getitem__([animal,mark])
            [sI,sO] = self.__getitem__([sire,mark])
            [gsI,gsO] = self.__getitem__([gsire,mark])
            if sI == sO or aI != aO or gsI != gsO or '0' in [aI,aO,sI,sO,gsI,gsO]: phases.append('-')
            elif aI == gsI: phases.append('i')
            else: phases.append('o')
        return phases
    
    def getSafePhase2(self,animal,sire,gsire):
        """ Returns I/O phases towards grandsire and granddam
        """
        phases = []
        for mark in self.mark.getMarkers():
            [aI,aO] = self.__getitem__([animal,mark])
            [sI,sO] = self.__getitem__([sire,mark])
            [gsI,gsO] = self.__getitem__([gsire,mark])
            if aI != aO or gsI != gsO or sI == sO or '0' in [aI,aO,sI,sO,gsI,gsO]: phases.append('-')
            elif aI == gsI: phases.append('I')
            else: phases.append('O')
        return phases
    
    def recomb(self,w):
        """ Returns True if w follows: [a,b,...,b,a] for 1..N b's
        also returns True if w is [a,b]
        """
        if len(w) == 2:
            if w[0] != w[1] and '-' not in w: return True
            else: return False
        if w[0] != w[-1] or '-' in w: return False
        if w[1] == w[0]: return False
        for i in xrange(2,len(w)-1):
            if w[i] != w[1]: return False
        return True
    
    def checkDoubles(self,animal,phases,size):
        """ Look for double recombinations involving 'size' markers
        Special case when size=0: looks for all recombinations, returns position before and after
        """
        def switch(t):
            if t == 'I': return 'i'
            if t == 'O': return 'o'
            return t
        sresult, dresult = [],[]
        swindow,dwindow = ['-']*(size+2),['-']*(size+2)
        spos,dpos = [0]*(size+2),[0]*(size+2)
        for pos in xrange(0,len(phases[0])):
            ph = switch(phases[0][pos])
            if ph not in ['-','x']:
                swindow.append(ph)
                swindow.pop(0)
                spos.append(pos)
                spos.pop(0)
                if self.recomb(swindow):
                    if size > 1:
                        sresult.append(spos[1])
                        sresult.append(spos[size])
                    elif size == 1: sresult.append(spos[1])
                    elif size == 0:
                        sresult.append(spos[0])
                        sresult.append(spos[1])
            # maternal phase
            ph = switch(phases[1][pos])
            if ph not in ['-','x']:
                dwindow.append(ph)
                dwindow.pop(0)
                dpos.append(pos)
                dpos.pop(0)
                if self.recomb(dwindow):
                    if size > 1:
                        dresult.append(dpos[1])
                        dresult.append(dpos[size+1])
                    elif size == 1: dresult.append(dpos[1])
                    elif size == 0:
                        dresult.append(dpos[0])
                        dresult.append(dpos[1])
        return sresult,dresult
        
    def checkDoublesSire(self,animal,phases,size):
        """ Look for double recombinations involving 'size' markers
        Special case when size=0: looks for all recombinations, returns position before and after
        """
        def switch(t):
            if t == 'I': return 'i'
            if t == 'O': return 'o'
            return t
        sresult = []
        swindow = ['-']*(size+2)
        spos = [0]*(size+2)
        for pos in xrange(0,len(phases)):
            ph = switch(phases[pos])
            if ph not in ['-','x']:
                swindow.append(ph)
                swindow.pop(0)
                spos.append(pos)
                spos.pop(0)
                if self.recomb(swindow):
                    if size > 1:
                        sresult.append(spos[1])
                        sresult.append(spos[size])
                    elif size == 1: sresult.append(spos[1])
                    elif size == 0:
                        sresult.append(spos[0])
                        sresult.append(spos[1])
        return sresult
    
    def writePhase(self,fout,testAnim='0'):
        """ Writes i/o-phases for the haplotypes, testAnim is used to limit output to only children of one sire/dam
            mode = 1: only mark homozygouse offspring
        """
        sep = '\t'
        testSex = '9'
        if testAnim != '0': testSex = self.ped.getSex(testAnim)
        for member in self.ped.getAnimals():
            sire,dam,sex = self.ped.getSire(member),self.ped.getDam(member),self.ped.getSex(member)
            if (sire == '0' and dam == '0') or (testAnim != '0' and testAnim not in [sire,dam]) or (sire not in self.ped.getAnimals()): continue
            [ph1,ph2] = self.getPhase(member,sire,dam)
            if ph1.count('x') > 0 or ph2.count('x') > 0:
                [nph2,nph1] = self.getPhase(member,dam,sire)
                if (nph1.count('x') + nph2.count('x')) < (ph1.count('x') + ph2.count('x')): ph1,ph2 = nph1,nph2
            sr,dr = self.checkDoubles(member,[ph1,ph2],0)
            s2r,d2r = self.checkDoubles(member,[ph1,ph2],1)
            try:
                if testSex in ['1','9']:
                    try:
                        fout.write(member+sep+sire+sep+'P'+sep+str(len(sr)/2)+sep+sep.join([''.join(ph1[i:i+50]) for i in xrange(0,len(ph1),50)])+sep+str(len(s2r)))
                    except TypeError:
                        pass
                        #print member,sire,family,sr,s2r,ph1
                        #sys.exit()
                    fout.write(sep+sep.join([str(el) for el in s2r]))
                    fout.write('\n')
                if testSex in ['0','9']:
                    fout.write(member+sep+dam+sep+'M'+sep+str(len(dr)/2)+sep+sep.join([''.join(ph2[i:i+50]) for i in xrange(0,len(ph2),50)])+sep+str(len(d2r)))
                    if len(d2r) > 0: fout.write(sep+sep.join([str(el) for el in d2r]))
                    fout.write('\n')
            except IOError:
                pass
                
    def writePhaseSafe(self,fout,testAnim='0'):
        """ Writes i/o-phases for the 'safe' haplotypes, testAnim is used to limit output to only children of one sire/dam
        """
        sep = '\t'
        testSex = '9'
        if testAnim != '0': testSex = self.ped.getSex(testAnim)
        for member in self.ped.getAnimals():
            sire,dam,sex = self.ped.getSire(member),self.ped.getDam(member),self.ped.getSex(member)
            try: gsire = self.ped.getSire(sire)
            except KeyError: continue
            if (sire == '0' and dam == '0') or (testAnim != '0' and testAnim not in [sire,dam]) or gsire == '0' or (sire not in self.ped.getAnimals()): continue
            ph1 = self.getSafePhase(member,sire,gsire)
            sr,dr = self.checkDoubles(member,[ph1,ph1],0)
            s2r,d2r = self.checkDoubles(member,[ph1,ph1],1)
            try:
                if testSex in ['1','9']:
                    try:
                        fout.write(member+sep+sire+sep+'P'+sep+str(len(sr)/2)+sep+sep.join([''.join(ph1[i:i+50]) for i in xrange(0,len(ph1),50)])+sep+str(len(s2r)))
                    except TypeError:
                        #pass
                        print member,sire,family,sr,s2r,ph1
                        sys.exit()
                    fout.write(sep+sep.join([str(el) for el in s2r]))
                    fout.write('\n')
            except IOError:
                #pass
                print 'Error'
                
    def writePhaseSafe2(self,fout):
        """ Writes i/o-phases for the 'safe' haplotypes, testAnim is used to limit output to only children of one sire/dam
        """
        sep = '\t'
        for member in self.ped.getAnimals():
            sire = self.ped.getSire(member)
            try: gsire = self.ped.getSire(self.ped.getSire(member))
            except KeyError: continue
            if gsire == '0': continue
            ph1 = self.getSafePhase2(member,sire,gsire)
            sr,dr = self.checkDoubles(member,[ph1,ph1],0)
            s2r,d2r = self.checkDoubles(member,[ph1,ph1],1)
            try:
                try:
                    fout.write(member+sep+sire+sep+'P'+sep+str(len(sr)/2)+sep+sep.join([''.join(ph1[i:i+50]) for i in xrange(0,len(ph1),50)])+sep+str(len(s2r)))
                except TypeError:
                    #pass
                    print member,sire,family,sr,s2r,ph1
                    sys.exit()
                fout.write(sep+sep.join([str(el) for el in s2r]))
                fout.write('\n')
            except IOError:
                #pass
                print 'Error'
                    
    def fixPhase(self):
        """ Corrects double recombinations caused by single markers
        """
        sep = '\t'
        testSex = '9'
        for member in self.ped.getAnimals():
            sire,dam,sex = self.ped.getSire(member),self.ped.getDam(member),self.ped.getSex(member)
            if sire == '0' and dam == '0': continue
            [ph1,ph2] = self.getPhase(member,sire,dam)
            s2r,d2r = self.checkDoubles(member,[ph1,ph2],1)
            for el in s2r:
                self.__setitem__([member,el],['0','0'])
            for el in d2r:
                self.__setitem__([member,el],['0','0'])
                                                                                                                                                                                                
    def cleanFile(self,flog=None):
        """ Removes haplotypes with only one allele and set haplotypes containing '9' (error) to '0'
        """
        if flog: flog.write('#Removing illegal alleles and double recombinants\n')
        for animal in self.ped.getAnimals():
            for marker in self.mark.getMarkers():
                [aI,aO] = self.__getitem__([animal,marker])
                if [aI,aO] == ['0','0']: continue
                if aI in ['9','0'] or aO in ['0','9']:
                    self.__setitem__([animal,marker],['0','0'])
                    if flog: flog.write('%s %s a[%s,%s] x[%s,%s]\n' % (member,marker,aI,aO,'0','0'))
                
    def scrub(self,flog = None):
        """ Finds illegal alleles and double recombinants and remove them
        """
        def oppo(a,m1,m2):
            if a == m1: return m2
            if a == m2: return m1
            return '0'
        
        if flog: flog.write('#Removing illegal alleles and double recombinants\n')
        sep = '\t'
        testSex = '9'
        changed = True
        while changed:
            changed = False
            for member in self.ped.getAnimals():
                sire,dam,sex = self.ped.getSire(member),self.ped.getDam(member),self.ped.getSex(member)
                if sire == '0' and dam == '0': continue
                [ph1,ph2] = self.getPhase(member,sire,dam)
                #if ph1.count('x') > 0 or ph2.count('x') > 0:
                #    [nph2,nph1] = self.getPhase(member,dam,sire)
                #    if (nph1.count('x') + nph2.count('x')) < (ph1.count('x') + ph2.count('x')): ph1,ph2 = nph1,nph2
                s2r,d2r = self.checkDoubles(member,[ph1,ph2],1) # Two lists with positions of double recombs
                #removing genotypes giving double recombinations
                for marker in self.mark.getMarkers():
                    pos = self.mark[marker][1]
                    sph = ph1[pos]
                    dph = ph2[pos]
                    if sph != 'x' and dph != 'x' and pos not in s2r and pos not in d2r: continue
                    m1,m2 = self.getMAlleles(marker)
                    siredouble = damdouble = sirex = damx = False
                    if pos in s2r: siredouble = True
                    if pos in d2r: damdouble = True
                    if ph1[pos] == 'x': sirex = True
                    if ph2[pos] == 'x': damx = True
                    [aI,aO] = self.__getitem__([member,marker])
                    if aI == '9': aI = '0'
                    if aO == '9': aO = '0'
                    [xI,xO] = self.findLegalAlleles(member,sire,dam,marker)
                    if xI != '0' and xO != '0':
                        legal = 0
                        if (xI == aI and not siredouble and not sirex) or (xI != aI and (siredouble or sirex)): legal += 1
                        if (xO == aO and not damdouble and not damx) or (xO != aO and (damdouble or damx)): legal += 1
                        if legal == 2:
                            self.__setitem__([member,marker],[xI,xO])
                            changed = True
                            if flog: flog.write('%s %s a[%s,%s] x[%s,%s]\n' % (member,marker,aI,aO,xI,xO))
                        else:
                            self.__setitem__([member,marker],['0','0'])
                            if flog: flog.write('%s %s a[%s,%s] x[%s,%s]\n' % (member,marker,aI,aO,'0','0'))
                        continue
                    if '0' in aI+aO: print '%s %s a[%s,%s] x[%s,%s]\n' % (member,marker,aI,aO,xI,xO)
                    if (siredouble or sirex) and xI in ['0',oppo(aI,m1,m2)]: xI = oppo(aI,m1,m2)
                    else: xI = aI
                    if (damdouble or damx) and xO in ['0',oppo(aO,m1,m2)]: xO = oppo(aO,m1,m2)
                    else: xO = aO
                    if '0' in xI+xO: print '%s %s a[%s,%s] x[%s,%s] m[%s,%s]\n' % (member,marker,aI,aO,xI,xO,m1,m2)
                    self.__setitem__([member,marker],[xI,xO])
                    if flog: flog.write('%s %s a[%s,%s] x[%s,%s]\n' % (member,marker,aI,aO,xI,xO))
                    changed = True
                    
    def writeBglMark(self,fout):
        def trans(s):
            if s in '1A': return 'A'
            if s in '2C': return 'C'
            if s in '3G': return 'G'
            if s in '4T': return 'T'
            if s in 'M': return 'M'
            return '0'
        
        sep = '\t'
        warned = False
        for marker in self.mark.getMarkers():
            a = [] 
            s = ''
            for animal in self.genotypes:
                try: allele = self.__getitem__([animal,marker])
                except (KeyError,IndexError): continue
                if '0' in allele: continue
                s += allele[0] + allele[1]
                if allele[0] not in a: a.append(allele[0])
                if allele[1] not in a: a.append(allele[1])
            pos = self.mark[marker][3]
            if len(a) < 2 and not warned:
                print "WARNING! Monomorphic markers!"
                warned = True
            if len(a) > 2:
                print "ERROR!",marker,a
                sys.exit(1)
            if len(a) == 0: a.append('M')
            if len(a) == 1: a.append('M')
            sc1 = s.count(a[0])
            sc2 = s.count(a[1])
            if sc1 < sc2: a = [a[1],a[0]]
            fout.write(marker+sep+str(pos)+sep+sep.join([trans(m) for m in a])+sep+str(sc1)+sep+str(sc2)+'\n')
        fout.close()
                        
    def testSpeed(self,fout):
        import time
        t = time.time()
        fout.write(str(self))
        print "Option 1: %.3f" % (time.time()-t)
        t = time.time()
        self.fastPrint(fout)
        print "Option 2: %.3f" % (time.time()-t)
        
    def alleleDist(self,fout):
        fout.write('#animal\tblank\thom1\thet\thom2\twrong\n')
        markinfo = {}
        for animal in self.ped.getAnimals():
            res = [0,0,0,0,0] # blank,hom1,het,hom2,wrong
            for marker in self.mark.getMarkers():
                if marker not in markinfo: markinfo[marker] = [0,0,0,0,0]
                m1,m2,MAF = self.getMAF(marker)
                a1,a2 = self.__getitem__([animal,marker])
                if a1 != a2 and a1 in m1+m2 and a2 in m1+m2: # Het
                    res[2] += 1
                    markinfo[marker][2] += 1
                elif '0' in a1+a2:
                    res[0] += 1
                    markinfo[marker][0] += 1
                elif a1 == m1:
                    res[1] += 1
                    markinfo[marker][1] += 1
                elif a1 == m2:
                    res[3] += 1
                    markinfo[marker][3] += 1
                else:
                    res[4] += 1
                    markinfo[marker][4] += 1
                    sys.stdout.write('(%s,%s):a(%s,%s) m(%s,%s)\n' % (animal,marker,a1,a2,m1,m2) )
            fout.write(animal+'\t'+'\t'.join([str(r) for r in res])+'\n')
        fout.write('#marker\tblank\thom1\thet\thom2\twrong\tHW-test\tMAF\n')
        for marker in self.mark.getMarkers():
            n = markinfo[marker][1]+markinfo[marker][2]+markinfo[marker][3]
            m1,m2,MAF = self.getMAF(marker)
            try:
                p = 1.0*(2*markinfo[marker][1]+markinfo[marker][2]) / (2*n)
                q = 1-p
                EAA,EAa,Eaa = p*p*n,2*p*q*n,q*q*n
                xstat = pow((markinfo[marker][1]-EAA),2)/EAA+pow((markinfo[marker][2]-EAa),2)/EAa+pow((markinfo[marker][3]-Eaa),2)/Eaa
            except ZeroDivisionError:
                xstat = 0
            fout.write(marker+'\t'+'\t'.join([str(r) for r in markinfo[marker]]))
            fout.write('\t%.5f\t%.5f\n' % (xstat,MAF))
    
    def calcStats(self,fout):
        data = {}
        sep = '\t'
        for animal in self.ped.getAnimals(): fout.write(sep+animal)
        fout.write('\n')
        for animal1 in self.ped.getAnimals():
            fout.write(animal1)
            for animal2 in self.ped.getAnimals():
                if animal1 == animal2 or (animal2,animal1) in data:
                    if animal1 != animal2:  fout.write(sep+str(data[animal2,animal1]))
                    else: fout.write(sep+'0')
                    continue
                dist = 0
                for marker in self.mark.getMarkers():
                    a1,a2 = self.__getitem__([animal1,marker])
                    b1,b2 = self.__getitem__([animal2,marker])
                    if '0' in [a1,a2,b1,b2]: continue
                    if a1 == a2 and b1 == b2 and a1 != b1: dist += 2
                    elif a1 == a2 and b1 != b2 or a1 != a2 and b1 == b2: dist += 1
                fout.write(sep+str(dist))
                data[animal1,animal2] = dist
            fout.write('\n')
        fout.close()
        
    def splitmark(self,numsplits,outfile):
        for i in xrange(0,numsplits):
            print i
            
def extractPedigree(infile):
    """ Returns a libPed object """
    #if not self.localMark:
    #    self.markIndex = [[self.mark[str(mark)][1]*2,self.mark[str(mark)][1]*2+1] for mark in self.mark['markers']]
    if infile[-3:] == '.gz': fgeno = gzip.open(infile,'r')
    else: fgeno = open(infile,'r')
    ped = libPed.Ped()
    for line in fgeno:
        if line.startswith('#'): continue
        try: animal,sire,dam,geno = line.strip().split(None,3)
        except ValueError:
            sys.stderr.write('Error: Wrong format of input file\n')
            sys.exit(1)
        ped.addAnimal(animal,dam,sire,'F0','3')
    ped.updateSex()
    fgeno.close()
    return ped
    
def extractMarkers(infile):
    """ Returns a libMark object """
    if infile[-3:] == '.gz': fgeno = gzip.open(infile,'r')
    else: fgeno = open(infile,'r')
    line = fgeno.next().strip().split()
    if line[0] != '#': line = ['#']+['M'+str(i) for i in xrange( (len(line)-3)/2 )]
    mark = libMark.Mark(line[1:])
    fgeno.close()
    return mark
    
def extractBglMark(infile):
    """ Returns a libMark object with markeralleles set """
    mark = extractMarkers(infile)
    if infile[-3:] == '.gz': fgeno = gzip.open(infile,'r')
    else: fgeno = open(infile,'r')
    alldict = {}
    marklist = mark.getMarkers()
    for marker in marklist: alldict[marker] = []
    for line in fgeno:
        if len(alldict) == 0: break
        if line.startswith('#'): continue
        l = line.strip().split()
        for i in xrange(0,len(marklist)):
            if marklist[i] not in alldict: continue
            a1,a2 = l[i*2+3:i*2+5]
            if a1 == '0': continue
            if a1 not in alldict[marklist[i]]: alldict[marklist[i]].append(a1)
            if a2 not in alldict[marklist[i]]: alldict[marklist[i]].append(a2)
            if len(alldict[marklist[i]]) == 2:
                mark['alleles',marklist[i]] = sorted(alldict[marklist[i]])
                alldict.pop(marklist[i])
    print "%d monomorphic or ungenotyped markers" % len(alldict)
    for marker in alldict:
        a = alldict[marker]
        if len(a) == 1:
            mark['alleles',marker] = (a[0],'0')
        elif len(a) == 0:
            mark['alleles',marker] = ('0','0')
        else:
            mark['alleles',marker] = ('NA',''.join(a))
            print "WARNING:",marker
    fgeno.close()
    return mark

def toLett(a):
    if a in 'A1': return 'A'
    elif a in 'C2': return 'C'
    elif a in 'G3': return 'G'
    elif a in 'T4': return 'T'
    elif a in 'D5': return 'D'
    elif a in 'I6': return 'I'
    else: return '0'

def writeBglMark(mark, fout):
    count = 0
    for marker in mark.getMarkers():
        m1,m2 = mark.getMAlleles(marker)
        if m1 == 'NA':
            fout.write('%s\t%d\t[%s]\n' % (marker,count,m2) )
        else:
            fout.write(marker+'\t'+str(count)+'\t'+toLett(m1)+'\t'+toLett(m2)+'\n')
        count += 1

def main():
    try:
        fin = open(sys.argv[1],'r')
        opt = sys.argv[2]
        fout = open(sys.argv[3],'w')
    except IOError:
        print "Error opening files"
        sys.exit(0)
    except:
        print """ ped/mark/switch/stat/del#/bglmark/adist\n"""
        sys.exit(0)
    if opt not in ['ped','mark','bglmark']: cr = Geno(sys.argv[1],'0','0')
    if opt == 'ped': fout.write(str(extractPedigree(sys.argv[1])))
    elif opt == 'mark': fout.write(str(extractMarkers(sys.argv[1])))
    elif opt == 'bglmark':
        mark = extractBglMark(sys.argv[1])
        writeBglMark(mark,fout)
    elif opt == 'switch': 
        cr.switchAlleles()
        cr.fastPrint(fout)
    elif opt[:3] == 'del': cr.removeGenotypes(fout,int(opt[3:]))
    elif opt == 'findmaf': cr.writeBglMark(fout)
    elif opt == 'speed': cr.testSpeed(fout)
    elif opt == 'stat': cr.calcStats(fout)
    elif opt == 'phase': cr.writePhase(fout,'0')
    elif opt == 'phasesafe': cr.writePhaseSafe2(fout)
    elif opt == 'clean':
        cr.cleanFile()
        cr.fastPrint(fout)
    elif opt == 'scrub':
        cr.scrub()
        cr.fastPrint(fout)
    elif opt == 'adist':
        cr.alleleDist(fout)
    elif opt[:9] == 'splitmark':
        cr.splitmark(int(opt[9:]),sys.argv[3])
    else: print "Illegal option"
    fin.close()
    fout.close()
    
class EmptyGeno(object):
    def __init__(self):
        """ if markfile and pedfile are not valid files, markobj and pedobj are generated locally and can be exported
        """
        pass

    def __getitem__(self, key):
        return ['0','0']

    def getMAlleles(self, marker):
        return '-1','-1'
            
    def __str__(self):
        return ''
        
    def getAlleles(self,animal,marker,opt = 'lett'):
        return '0','0'
        
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
