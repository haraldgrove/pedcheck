#!/usr/bin/env python

# snpstat, version 1.0, 2014-03-14

from __future__ import division, print_function
import sys
import argparse
import numpy as np

class SNP(object):

    def __init__(self,ic,ia,informat,verbose):
        """
            ic = number of information columns before genotypes start
            ia = allele representation, 1 = 0/1/2, 2 = 11/13/33, 3 = 1 1/1 3/3 3
                 for 2 and 3, alleles are represented as 1/2/3/4
        """
        self.verbose = verbose
        self.sep = '\t'
        if not informat:
            self.ic = ic
            self.ia = ia
            self.nc = 0
        elif informat.lower() == 'plink':
            self.ic = 6
            self.ia = 3
            self.nc = 1
        elif informat.lower() == 'dmu':
            self.ic = 1
            self.ia = 1
            self.nc = 0
        elif informat.lower() == 'linkage':
            self.ic = 6
            self.ia = 3
            self.nc = 1
        else:
            sys.stderr.write('Unknown input format: "%s"\n' % informat)
            sys.exit(1)

    def tbasea(self,a,mark):
        # Transforms from 0/1/2 to 1/2/3/4
        m1 = self.mark[mark]['alleles'][0]
        m2 = self.mark[mark]['alleles'][1]
        if a == 1: return m1+m2
        if a == 0: return m1+m1
        if a == 2: return m2+m2
        return '00'

    def tbase012(self,a,mark):
        # Transform from 1/2/3/4 to 0/1/2
        if len(a) == 1: a = a+a
        try:
            m1 = self.mark[mark]['alleles'][0]
        except IndexError:
            if a[0] != a[1]:
                self.mark[mark]['alleles'] = [a[0],a[1]]
                m1 = a[0]
            elif a[0] != '0':
                self.mark[mark]['alleles'].append(a[0])
                m1 = a[0]
            else:
                m1 = 'X'
        try:
            m2 = self.mark[mark]['alleles'][1]
        except IndexError:
            if a[0] != '0' and a[0] != m1:
                self.mark[mark]['alleles'].append(a[0])
                m2 = a[0]
            elif a[1] != '0' and a[1] != m1:
                self.mark[mark]['alleles'].append(a[1])
                m2 = a[1]
            else:
                # The second allelel has not been encountered yet.
                m2 = 'X'
        if a[0] != a[1]: return '1'
        if a[0] == m1: return '0'
        if a[0] == m2: return '2'
        if a[0] != '0': sys.stderr.write('ERROR: Marker %s has more than 2 alleles: [%s]\n' % mark,str(a))
        return np.nan

    def readGenos(self,genofile,referencefile=None):
        """
        Reads the genotype file and converts the genotypes into a numpy array as 0/1/2
        The -a parameter will specify if the alleles are given as:
          0/1/2 (-a 1),
          11/13/33 (-a 2),
          1 1/1 3/3 3 (-a 3)
        """
        if referencefile: self.gen = np.zeros((len(self.ped1)+len(self.ped2),len(self.mark)))
        else: self.gen = np.zeros((len(self.ped1),len(self.mark)))
        self.gen[:] = np.nan
        with open(genofile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mlist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                irow = self.ped1[l[self.nc]]['rank']
                for i,mark  in enumerate(mlist):
                    if mark not in self.mark: continue
                    icol = self.mark[mark]['rank']
                    if self.ia == 1:
                        a = l[i+self.ic]
                    elif self.ia == 2: 
                        a = self.tbase012(l[i+self.ic],mark)
                    elif self.ia == 3:
                        a = self.tbase012(l[i*2+self.ic]+l[i*2+1+self.ic],mark)
                    if a not in ['0','1','2']: a = np.nan
                    else: a = int(a)
                    self.gen[irow,icol] = a-1
        if not referencefile: return
        with open(referencefile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mlist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                irow = self.ped2[l[self.nc]]['rank']
                for i,mark  in enumerate(mlist):
                    if mark not in self.mark: continue
                    icol = self.mark[mark]['rank']
                    if self.ia == 1:
                        a = l[i+self.ic]
                    elif self.ia == 2:
                        a = self.tbase012(l[i+self.ic],mark)
                    elif self.ia == 3:
                        a = self.tbase012(l[i*2+self.ic]+l[i*2+1+self.ic],mark)
                    if a not in ['0','1','2']: a = np.nan
                    else: a = int(a)
                    self.gen[irow,icol] = a-1

    def readPedigree(self,pedfile,count=0,real=True):
        """
        Reads a pedigree from either a separate pedigree file or from the given genotype file (real=False)
        The first 3 columns have to be: sample_name, father and mother, regardless of input file
        If the number of information columns are either 1 or 2, it will assume that no parents are
            present or just the father, respectively.
        count is the starting position for the first sample, use when reading from more than one file
        """
        ped = {}
        pedlist = []
        with open(pedfile,'r') as fin:
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                name,father,mother,family,sex = '0','0','0','0','3'
                if len(l) > 0: name = l[self.nc]
                if (real or self.ic > 1) and len(l) > 1: father = l[self.nc+1]
                if (real or self.ic > 2) and len(l) > 2: mother = l[self.nc+2]
                if real and len(l) > 4: family,sex = l[3],l[4]
                if name == '0': continue
                if name not in ped:
                    ped[name] = {'father':father,
                                 'mother':mother,
                                 'rank':count,
                                 'sex':sex,
                                 'children':[]}
                    count += 1
                    pedlist.append(name)
                else:
                    sys.stderr.write('%s present more than once\n' % name)
        self.updatePed(ped)
        return ped,pedlist

    def collectPedigree(self,pedfile,file1,file2):
        """
        Reads a pedigree and lists of the animals in each file
        ped and pedlist are pedigree information for the animals to be checked
        ped1 and names1 are samples from query file
        ped2 and names2 are samples from reference file
        """
        if pedfile:
            self.ped,self.pedlist = self.readPedigree(pedfile)
            self.ped1,self.names1 = self.readPedigree(file1,0,False)
        else:
            self.ped,self.pedlist = self.readPedigree(file1,0,False)
            self.ped1,self.names1 = self.ped,self.pedlist
        if file2: self.ped2,self.names2 = self.readPedigree(file2,len(self.ped1),False)
        else: self.ped2,self.names2 = self.ped1,self.names1


    def updatePed(self,ped):
        # Assign children to parents and set sex of parents
        for n in ped:
            father,mother = ped[n]['father'],ped[n]['mother']
            if father in ped:
                ped[father]['sex'] = '1'
                ped[father]['children'].append(n)
            if mother in ped:
                ped[mother]['sex'] = '0'
                ped[mother]['children'].append(n)

    def readMarkers(self,markerfile):
        """
            Read a marker file, format is Plink map-file
            Will also accept the addition of the marker alleles as two extra columns
            and a marker file containing just one column of marker names.
        """
        self.mark = {}
        self.marklist = []
        with open(markerfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                if len(l) >= 6: chrom,name,distance,position,a1,a2 = l[:6]
                elif len(l) == 4:
                    chrom,name,distance,position = l # Plink
                    a1,a2 = [],[]
                elif len(l) == 1:
                    name = l[0]
                    chrom,position,a1,a2 = '0',count,[],[]
                if name not in self.mark:
                    self.mark[name] = {'chrom':chrom,
                                       'pos':int(position),
                                       'alleles': a1+a2,
                                       'rank':count}
                    count += 1
                    self.marklist.append(name)

    def collectMarkers(self, ingeno):
        """ 
            Reads input file to search for the necessary marker information
            If the markers are not present as a comment line on top of the file,
            it will calculate the number of markers to be the same as the number of alleles
            after the information columns.
        """
        self.mark = {}
        self.marklist = []
        with open(ingeno,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    l = line.strip('#').strip().split()
                    for i,e in enumerate(l):
                        self.mark[e] = {'chrom':'0',
                                       'pos':i,
                                       'alleles': [],
                                       'rank':i}
                        self.marklist.append(e)
                    break
            else:
                l = line.strip().split()
                if self.ia == 3:
                    for i in xrange(0,len(l[self.ic:])//2):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': [],
                                           'rank':i}
                        self.marklist.append(str(i))
                else:
                    for i,e in enumerate(l[self.ic:]):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': [],
                                           'rank':i}
                        self.marklist.append(str(i))

#*****************************************************************************************************

    def findDiscords(self,anim,father,mother,limit=100):
        """ 
            Detects all mendelian discords within the given trio
            Will work even if one of the parents is missing
        """
        results = {}
        if self.verbose: sys.stdout.write('Checking: %s vs. %s and %s' % (anim,father,mother))
        # Father
        try:
            resF = self.gen[self.ped2[father]['rank'],:]*self.gen[self.ped1[anim]['rank'],:]
            wrongF = resF==-1
            sitesF = len(resF) - np.count_nonzero(np.isnan(resF))
        except KeyError:
            wrongF = [False]*len(self.mark)
            sitesF = 0
            resF = wrongF
        if sitesF == 0:
            pct = -1
        else:
            pct = 100*np.count_nonzero(wrongF) / sitesF
        if self.verbose: sys.stdout.write(' %.3f' % pct)
        results[father] = {'disc':np.count_nonzero(wrongF),'sites':sitesF,'pct':pct}
        # Mother
        try:
            resM = self.gen[self.ped2[mother]['rank'],:]*self.gen[self.ped1[anim]['rank'],:]
            wrongM = resM==-1
            sitesM = len(resM) - np.count_nonzero(np.isnan(resM))
        except KeyError:
            wrongM = [False]*len(self.mark)
            sitesM = 0
            resM = wrongM
        if sitesM == 0:
            pct = -1
        else:
            pct = 100*np.count_nonzero(wrongM) / sitesM
        if self.verbose: sys.stdout.write(' %.3f' % pct)
        results[mother] = {'disc':np.count_nonzero(wrongM),'sites':sitesM,'pct':pct}
        # Trios
        try:
            resT = self.gen[self.ped2[father]['rank'],:]*self.gen[self.ped2[mother]['rank'],:] -\
                   (self.gen[self.ped1[anim]['rank'],:]*self.gen[self.ped1[anim]['rank'],:])
            wrongT = resT==1
            sitesT = len(resT) - np.count_nonzero(np.isnan(resT))
        except KeyError:
            wrongT = [False]*len(self.mark)
            sitesT = 0
            resT = wrongT
        wrongTrio = np.logical_or(np.logical_or(wrongF,wrongM),wrongT)
        wrongSites = np.logical_or(np.logical_or(np.isnan(resF),np.isnan(resM)),np.isnan(resT))
        if len(resT)-np.count_nonzero(wrongSites) == 0:
            pct = -1
        else:
            pct = 100*np.count_nonzero(wrongTrio) / (len(resT)-np.count_nonzero(wrongSites))
        if self.verbose: sys.stdout.write(' %.3f\n' % pct)
        if pct > limit or pct < 0: return None
        results[father,mother] = {'disc':np.count_nonzero(wrongTrio),'sites':len(resT)-np.count_nonzero(wrongSites),'pct':pct}
        return results

    def checkPed(self,outfile):
        # Checks the provided pedigree for discords
        if outfile: fout = open(outfile,'w')
        else: fout = sys.stdout
        for sample in self.pedlist:
            father = self.ped[sample]['father']
            mother = self.ped[sample]['mother']
            if father != '0' or mother != '0':
                res = self.findDiscords(sample,father,mother)
                if not res:
                    f1 = f2 = m1 = m2 = t1 = t2 = 0
                    f3 = m3 = t3 = -1
                else:
                    f1,f2,f3 = res[father]['disc'],res[father]['sites'],res[father]['pct']
                    m1,m2,m3 = res[mother]['disc'],res[mother]['sites'],res[mother]['pct']
                    t1,t2,t3 = res[father,mother]['disc'],res[father,mother]['sites'],res[father,mother]['pct']
            else:
                f1 = f2 = m1 = m2 = t1 = t2 = 0
                f3 = m3 = t3 = -1 
            fout.write('%s\t%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n' % (sample,father,f1,f2,f3,mother,m1,m2,m3,t1,t2,t3))
        if outfile: fout.close()

    def findSingleDiscords(self,anim,father,limit):
        """
            Detects all mendelian discords within the given trio
            Will work even if one of the parents is missing
        """
        if self.verbose: sys.stdout.write('Checking: %s vs. %s ' % (anim,father))
        results = {}
        # Father
        try:
            resF = self.gen[self.ped2[father]['rank'],:]*self.gen[self.ped1[anim]['rank'],:]
            wrongF = resF==-1
            sitesF = len(resF) - np.count_nonzero(np.isnan(resF))
        except KeyError:
            wrongF = [False]*len(self.mark)
            sitesF = 0
            resF = wrongF
        if sitesF == 0:
            pct = -1
        else:
            pct = 100*np.count_nonzero(wrongF) / sitesF
        if self.verbose: sys.stdout.write('%d\n' % pct)
        if pct > limit or pct < 0: return None
        results = {'disc':np.count_nonzero(wrongF),'sites':sitesF,'pct':pct}
        return results

    def findParent(self,outfile,limit):
        if outfile: fout = open(outfile,'w')
        else: fout = sys.stdout
        # Checks the provided pedigree for potential parents
        for sample in self.pedlist:
            for parent in self.ped2:
                if sample == parent: continue
                temp  = self.findSingleDiscords(sample,parent,limit)
                if not temp: continue
                f1,f2,f3 = temp['disc'],temp['sites'],temp['pct']
                fout.write('%s\t%s\t%d\t%d\t%.2f\n' % (sample,parent,f1,f2,f3))
        if outfile: fout.close()

    def findPed(self,outfile,limit):
        res = {} # All 1-to-1 pairings under the given limit
        parents = {} # List of potential parents for each sample
        if not outfile: fout = sys.stdout
        else: fout = open(outfile,'w')
        # Locate all single parents
        for sample in self.pedlist:
            parents[sample] = []
            for parent in self.ped2:
                if sample == parent: continue
                if (parent,sample) in res or (sample,parent) in res:
                    parents[sample].append(parent)
                    continue
                temp  = self.findSingleDiscords(sample,parent,limit)
                if not temp: continue
                res[sample,parent] = temp
                parents[sample].append(parent)
        # Search single parents to find couples
        for sample in self.pedlist:
            foundCouple = False
            if len(parents[sample]) == 0: continue
            for i,parent1 in enumerate(parents[sample]):
                if len(parents[sample]) == 1:
                    continue
                for parent2 in parents[sample][i+1:]:
                    sex1,sex2 = self.ped2[parent1]['sex'],self.ped2[parent2]['sex']
                    if sex1 == '1' and sex2 == '0' or sex1 == '1' and sex2 == '3' or sex1 == '3' and sex2 == '0': father,mother = parent1,parent2
                    elif sex1 == '0' and sex2 == '1' or sex1 == '0' and sex2 == '3' or sex1 == '3' and sex2 == '1': father,mother = parent2,parent1
                    elif sex1 == sex2 and sex1 != '3': continue
                    else:
                        father,mother = parent1,parent2
                    temp = self.findDiscords(sample,father,mother,limit*2)
                    if not temp: continue
                    f1,f2,f3 = temp[father]['disc'],temp[father]['sites'],temp[father]['pct']
                    m1,m2,m3 = temp[mother]['disc'],temp[mother]['sites'],temp[mother]['pct']
                    t1,t2,t3 = temp[father,mother]['disc'],temp[father,mother]['sites'],temp[father,mother]['pct']
                    fout.write('%s\t%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n' % (sample,father,f1,f2,f3,mother,m1,m2,m3,t1,t2,t3))
                    foundCouple = True
            if not foundCouple:
                for parent in parents[sample]:
                    try: temp = res[sample,parent]
                    except KeyError: temp = res[parent,sample]
                    if self.ped2[parent]['sex'] == '0':
                        m1,m2,m3 = temp['disc'],temp['sites'],temp['pct']
                        f1 = f2 = t1 = t2 = 0
                        f3 = t3 = -1
                        father,mother = '',parent
                    else:
                        f1,f2,f3 = temp['disc'],temp['sites'],temp['pct']
                        m1 = m2 = t1 = t2 = 0
                        m3 = t3 = -1
                        father,mother = parent,'0'
                    fout.write('%s\t%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n' % (sample,father,f1,f2,f3,mother,m1,m2,m3,t1,t2,t3))
        
    def findDif(self,anim,father,limit):
        """
            Detects all difference within the given trio
            Will work even if one of the parents is missing
        """
        if self.verbose: sys.stdout.write('Checking: %s vs. %s ' % (anim,father))
        results = {}
        # Father
        try:
            resF = self.gen[self.ped2[father]['rank'],:]-self.gen[self.ped[anim]['rank'],:]
            wrongF = resF!=0
            sitesF = len(resF) - np.count_nonzero(np.isnan(resF))
        except KeyError:
            wrongF = [False]*len(self.mark)
            sitesF = 0
            resF = wrongF
        if sitesF == 0:
            pct = -1
        else:
            pct = 100*np.count_nonzero(wrongF) / sitesF
        if self.verbose: sys.stdout.write('%d/%d\n' % (pct,limit))
        if pct < limit: return None
        results = {'disc':np.count_nonzero(wrongF),'sites':sitesF,'pct':pct}
        return results
        
    def findDup(self,outfile,limit):
        if outfile: fout = open(outfile,'w')
        else: fout = sys.stdout
        lim = 100-limit
        # Checks the provided pedigree for potential parents
        for i,sample1 in enumerate(self.pedlist):
            for sample2 in self.pedlist[i+1:]:
                temp  = self.findDif(sample1,sample2,lim)
                if not temp:
                    if self.verbose: sys.stdout.write('%s and %s gave no results\n' % (sample1,sample2))
                    continue
                f1,f2,f3 = temp['disc'],temp['sites'],temp['pct']
                if f3 > lim: fout.write('%s\t%s\t%d\t%d\t%.2f\n' % (sample1,sample2,f1,f2,f3))
        if outfile: fout.close()

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('-r','--repfile',help='Output report file')
    parser.add_argument('-o','--mode',help='Type of operation (check/findparent/findped/dup)',default='check')
    parser.add_argument('-j','--reference',help='File with potential parents')
    parser.add_argument('-p','--pedigree',help='Pedigree file',default=None)
    parser.add_argument('-m','--markers',help='Marker file')
    parser.add_argument('-n','--informat',help='Format of input file (Plink/DMU/Linkage)')
    parser.add_argument('-c',dest='infocol',type=int,help='Non-genotype columns', default=3)
    parser.add_argument('-a',dest='allele',type=int,help='Alleleformat, 1=0/1/2, 2=11/13/33, 3=1 1/1 3/3 3', default=3)
    parser.add_argument('-l',dest='limit',type=float,help='Noise level in percent', default = 100.0)
    parser.add_argument('-v','--verbose',action="store_true",help='Prints runtime info')
    args = parser.parse_args()
    gen = SNP(args.infocol,args.allele,args.informat,args.verbose)
    # Read pedigree information, collect from genotype file if needed
    gen.collectPedigree(args.pedigree,args.ingeno,args.reference)
    # Read marker information, collect from genotype file if needed
    if args.markers:
        gen.readMarkers(args.markers)
    else:
        gen.collectMarkers(args.ingeno)
    gen.readGenos(args.ingeno,args.reference)
    if args.mode.lower() == 'check': gen.checkPed(args.repfile)
    elif args.mode.lower() == 'findparent': gen.findParent(args.repfile,args.limit)
    elif args.mode.lower() == 'findped': gen.findPed(args.repfile,args.limit)
    elif args.mode.lower() == 'dup': gen.findDup(args.repfile,args.limit)
    else: sys.stderr.write('ERROR: Unknown mode [%s]\n' % args.mode)

if __name__ == '__main__':
    main()
