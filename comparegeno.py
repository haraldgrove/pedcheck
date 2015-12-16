#!/usr/bin/env python

# snpstat, version 1.0, 2014-03-14

from __future__ import division, print_function
import sys
import argparse
import gzip
import numpy as np

class SNP(object):

    def __init__(self,informat,verbose):
        """
            ic = number of information columns before genotypes start
            ia = allele representation, 1 = 0/1/2, 2 = 11/13/33, 3 = 1 1/1 3/3 3
                 for 2 and 3, alleles are represented as 1/2/3/4
        """
        self.verbose = verbose
        self.sep = '\t'
        if informat.lower() == 'geno':
            self.ic = 3
            self.ia = 3
            self.nc = 0
        elif informat.lower() == 'plink':
            self.ic = 6
            self.ia = 3
            self.nc = 1
        elif informat.lower() == 'dmu':
            self.ic = 1
            self.ia = 1
            self.nc = 0
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
        if a[0] != '0': sys.stderr.write('ERROR: Marker %s has more than 2 alleles: [%s]\n' % (mark,str(a)))
        return np.nan

    def readGenos(self,genofile,referencefile=None):
        """
        Reads the genotype file and converts the genotypes into a numpy array as 0/1/2
        The -a parameter will specify if the alleles are given as:
          0/1/2 (-a 1),
          11/13/33 (-a 2),
          1 1/1 3/3 3 (-a 3)
        """
        if referencefile:
            self.gen = np.zeros((len(self.ped),len(self.mark)))
            self.gen1 = np.zeros((len(self.ped),len(self.mark)))
            self.gen1[:] = np.nan
        else: self.gen = np.zeros((len(self.ped),len(self.mark)))
        self.gen[:] = np.nan
        if genofile.rsplit('.',1)[1] == 'gz': op = gzip.open
        else: op = open
        with op(genofile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mlist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                try: irow = self.ped[l[self.nc]]['rank']
                except KeyError: continue
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
        if referencefile.rsplit('.',1)[1] == 'gz': op = gzip.open
        else: op = open
        with op(referencefile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mlist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                try: irow = self.ped[l[self.nc]]['rank']
                except KeyError: continue
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
                    self.gen1[irow,icol] = a-1

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
        if pedfile.rsplit('.',1)[1]  == 'gz': op = gzip.open
        else: op = open
        with op(pedfile,'r') as fin:
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
        if file2: self.ped2,self.names2 = self.readPedigree(file2,0,False)
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

    def readMarkers(self,markerfile,chro=None):
        """
            Read a marker file, format is Plink map-file
            Will also accept the addition of the marker alleles as two extra columns
            and a marker file containing just one column of marker names.
        """
        self.mark = {}
        self.marklist = []
        if markerfile.rsplit('.',1)[1]  == 'gz': op = gzip.open
        else: op = open
        with op(markerfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                if len(l) == 6: chrom,name,distance,position,a1,a2 = l
                elif len(l) == 4:
                    chrom,name,distance,position = l # Plink
                    a1,a2 = [],[]
                elif len(l) == 1:
                    name = l[0]
                    chrom,position,a1,a2 = '0',count,[],[]
                if chro and chrom != chro: continue
                if name not in self.mark:
                    self.mark[name] = {'chrom':chrom,
                                       'pos':int(position),
                                       'alleles': a1+a2,
                                       'rank':count}
                    count += 1
                    self.marklist.append(name)

#*****************************************************************************************************

    def findDif(self,anim,father):
        """
            Detects all difference
        """
        if self.verbose: sys.stdout.write('Checking: %s vs. %s ' % (anim,father))
        results = {}
        try:
            s1 = self.gen1[self.ped[father]['rank'],:]
            s2 = self.gen[self.ped[anim]['rank'],:]
            resF = s1-s2
            try:
                corrF = np.corrcoef(s1[np.isfinite(resF)],s2[np.isfinite(resF)])[0,1]
            except IndexError:
                corrF = np.nan
            rightF = resF==0
            sitesF = len(resF) - np.count_nonzero(np.isnan(resF))
        except KeyError:
            rightF = [False]*len(self.mark)
            sitesF = 0
            resF = rightF
            corrF = np.nan
        if sitesF == 0:
            pct = -1
        else:
            pct = 100*(np.count_nonzero(rightF) / sitesF)
        if self.verbose: sys.stdout.write('%d\n' % pct)
        results = {'disc':np.count_nonzero(rightF),'sites':sitesF,'pct':pct,'corr':corrF}
        return results
        
    def findDup(self,outfile):
        if outfile: fout = open(outfile,'w')
        else: fout = sys.stdout
        # Checks the provided pedigree for differences between the files
        for sample1 in self.pedlist:
            temp  = self.findDif(sample1,sample1)
            if not temp: continue
            f1,f2,f3,f4 = temp['disc'],temp['sites'],temp['pct'],temp['corr']
            if np.isnan(f4):
                f4 = 'nan'
            else:
                f4 = '%.4f' % f4
            fout.write('%s\t%d\t%d\t%.2f\t%s\n' % (sample1,f1,f2,f3,f4))
        if outfile: fout.close()

    def findDifMark(self,mark):
        results = {}
        rc = self.mark[mark]['rank']
        s1 = self.gen1[:,rc]
        s2 = self.gen[:,rc]
        resF = s1-s2
        try:
            corrF = np.corrcoef(s1[np.isfinite(resF)],s2[np.isfinite(resF)])[0,1]
        except IndexError:
            corrF = np.nan
        sitesF = len(resF) - np.count_nonzero(np.isnan(resF))
        rightF = resF==0
        if sitesF == 0:
            pct = -1
        else:
            pct = 100*(np.count_nonzero(rightF) / sitesF)
        results = {'disc':np.count_nonzero(rightF),'sites':sitesF,'pct':pct,'corr':corrF}
        return results

    def findMark(self,outfile):
        if outfile: fout = open(outfile,'w')
        else: fout = sys.stdout
        for mark1 in self.marklist:
            temp = self.findDifMark(mark1)
            if not temp: continue
            f1,f2,f3,f4 = temp['disc'],temp['sites'],temp['pct'],temp['corr']
            if np.isnan(f4):
                f4 = 'nan'
            else:
                f4 = '%.4f' % f4
            fout.write('%s\t%d\t%d\t%.2f\t%s\n' % (mark1,f1,f2,f3,f4))
        if outfile: fout.close()

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('-r','--repfile',help='Output report file')
    parser.add_argument('-j','--reference',help='File with the true genotypes')
    parser.add_argument('-p','--pedigree',help='Pedigree file',required=True)
    parser.add_argument('-t','--type',action="store_true",help='Compare markers instead of animals')
    parser.add_argument('-m','--markers',help='Marker file',required=True)
    parser.add_argument('-n','--informat',help='Format of input file (Plink/DMU/Geno)',default='Geno')
    parser.add_argument('-c','--chrom',help='Chromosome to work on')
    parser.add_argument('-v','--verbose',action="store_true",help='Prints runtime info')
    args = parser.parse_args()
    gen = SNP(args.informat,args.verbose)
    gen.collectPedigree(args.pedigree,args.ingeno,args.reference)
    gen.readMarkers(args.markers,args.chrom)
    gen.readGenos(args.ingeno,args.reference)
    if args.type: gen.findMark(args.repfile)
    else: gen.findDup(args.repfile)

if __name__ == '__main__':
    main()
