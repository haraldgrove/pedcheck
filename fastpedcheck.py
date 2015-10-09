#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2012
# This is GNU GPL Software: http://www.gnu.org/
# Written with vim

import sys
import time
import libMark
import libPed
import libGeno
from optparse import OptionParser
import numpy as np
import gzip
import math

def tbase012(a1,a2,A,B):
	res = ''
	if a1 != a2: return 1
	elif a1 == a2 == '0': return 'nan'
	elif a1 == A: return 0
	elif a1 == B: return 2
	else:
		print "Wrong allele",a1,a2,A,B
		sys.exit(1)

def tbasenum(a):
	if a in 'A1': return '1'
	elif a in 'C2': return '2'
	elif a in 'G3': return '3'
	elif a in 'T4': return '4'
	return '0'

def pedcheck(options):
	#if options.genotypefile[-3:0] == '.gz': fin = gzip.open(options.genotypefile,'r')
	#else: fin = open(options.genotypefile,'r')
	#***** Check requirements and read data *****
	if options.pedigree: pedigree = libPed.Ped(options.pedigree)
	else:
		print "Gathering pedigree from data"
		pedigree = libGeno.extractPedigree(options.genofile)
		#sys.stderr.write('Pedigree file needed.\n')
		#sys.exit(1)
	if options.markers: markers = libMark.Mark(options.markers)
	else:
		print "Gathering markers from data"
		markers = libGeno.extractBglMark(options.genofile)
		#sys.stderr.write('Marker file needed.\n')
		#sys.exit(1)
	checkAll,checkOrphans = False,False
	if ',' in options.pedlims: lim = options.pedlims.split(',')
	else: lim = (options.pedlims,options.pedlims)
	oldLim,newLim = float(lim[0]),float(lim[1])
	r = {}
	hits = {}
	count = 0
	for anim in pedigree.getAnimals():
		r[anim] = count
		count += 1
	#****** Reading data, converting if needed
	gen2 = None
	fout = None
	if options.reportfile: fout = open(options.reportfile,'w')
	fin = open(options.genofile,'r')
	rows = len(pedigree)
	columns = len(markers)
	gen = np.zeros((rows,columns))
	marklist = fin.next().strip('#').strip().split()
	newanims = {}
	for line in fin:
		l = line.strip().split()
		icolumn = 0
		irow = r[l[0]]
		newanims[l[0]] = 1
		for i in xrange(0,len(marklist)*2,2):
			a1,a2 = l[i+3:i+5]
			try: m1,m2 = markers.getMAlleles(marklist[icolumn])
			except KeyError:
				if marklist[icolumn] not in markers.getMarkers():
					print "WARNING! Incomplete markerlist"
					a1,a2 = '0','0'
				else:
					print "ERROR! Failed to assign marker alleles",marklist[icolumn]
					sys.exit(1)
			if m1==m2 or '0' in m1+m2: gen[irow,icolumn] = 1
			else: gen[irow,icolumn] = tbase012(a1,a2,tbasenum(m1),tbasenum(m2))
			icolumn += 1
	fin.close()
	if options.genofile2:
		fin = open(options.genofile2,'r')
		for line in fin:
			if line.startswith('#'): continue
			l = line.strip().split()
			irow = r[l[0]]
			icolumn = 0
			for i in xrange(0,len(marklist)*2,2):
				a1,a2 = l[i+3:i+5]
				try: m1,m2 = markers.getMAlleles(marklist[icolumn])
				except KeyError:
					if marklist[icolumn] not in markers.getMarkers():
						print "WARNING! Incomplete markerlist"
						a1,a2 = '0','0'
					else:
						print "ERROR! Failed to assign marker alleles"
						sys.exit(1)
				if m1 == m2 or '0' in m1+m2: gen[irow,icolumn] = 1
				else: gen[irow,icolumn] = tbase012(a1,a2,tbasenum(m1),tbasenum(m2))
				icolumn += 1
			irow += 1
		fin.close()
	out = ''
	sep = '\t'
	if fout: fout.write('#ID\tparent\tdiscords\tinfo_sites\tdiscord%\tcategory_sex\n')
	for anim in pedigree.getAnimals():
		if anim not in newanims: continue
		mismatch = False
		sire,dam = pedigree.getSire(anim),pedigree.getDam(anim)
		if sire != '0' and sire in r and not checkAll:
			res = gen[r[sire],:]-gen[r[anim],:]
			wrong = len(res[res==2]) + len(res[res==-2])
			identical = len(res[res==0])
			info = len(res)-len(res[np.isnan(res)])
			if info == 0: out += anim+sep+sire+sep+str(wrong)+sep+str(info)+sep+'-1'+sep+'1'+'\n'
			else: out += anim+sep+sire+sep+str(wrong)+sep+str(info)+sep+str(wrong*100.0/info)+sep+'1'+'\n'
			if info == 0: pass
			elif wrong*100.0/info > oldLim: 
				mismatch = True
				if fout: fout.write('%s\t%s\t%d\t%d\t%.3f\t%s\n' % (anim,sire,wrong,info,wrong*100.0/info,'W1'))
			hits[anim,sire] = wrong,info
		if dam != '0' and dam in r and not checkAll:
			res = gen[r[dam],:]-gen[r[anim],:]
			wrong = len(res[res==2]) + len(res[res==-2])
			identical = len(res[res==0])
			info = len(res)-len(res[np.isnan(res)])
			if info == 0: out += anim+sep+dam+sep+str(wrong)+sep+str(info)+sep+'-1'+sep+'0'+'\n'
			else: out += anim+sep+dam+sep+str(wrong)+sep+str(info)+sep+str(wrong*100.0/info)+sep+'0'+'\n'
			if info == 0: pass
			elif wrong*100.0/info > oldLim:
				mismatch = True
				if fout: fout.write('%s\t%s\t%d\t%d\t%.3f\t%s\n' % (anim,dam,wrong,info,wrong*100.0/info,'W0'))
			hits[anim,dam] = wrong,info
		if sire == '0' and dam == '0' and len(pedigree.getOffspring(anim)) == 0: mismatch = True
		if mismatch: # Search for better matches
			for anim2 in pedigree.getAnimals():
				if anim2 == anim: continue
				sex = pedigree.getSex(anim2)
				if (anim,anim2) in hits:
					wrong,info = hits[anim,anim2]
					rep = True
				elif (anim2,anim) in hits: 
					wrong,info = hits[anim2,anim]
					rep = True
				else:
					res = gen[r[anim2],:]-gen[r[anim],:]
					wrong = len(res[res==2]) + len(res[res==-2])
					identical = len(res[res==0])
					info = len(res)-len(res[np.isnan(res)])
					rep = False
					hits[anim,anim2] = wrong,info
				if info == 0 or rep: pass
				elif wrong*100.0/info <= newLim:
					if fout: fout.write('%s\t%s\t%d\t%d\t%.3f\t%s\n' % (anim,anim2,wrong,info,wrong*100.0/info,'N'+sex))
	if fout: fout.close()
	if len(out) > 0 and options.reportped:
		fout = open(options.reportped,'w')
		fout.write('#ID\tparent\tdiscords\tinfo_sites\tdiscord%\tparent_sex\n'+out)
		fout.close()

def main():
	usage = """usage: %prog [options]"""
	parser = OptionParser(usage)
	parser.add_option("-i",dest="genofile",help="input file", metavar="FILE")
	parser.add_option("-j",dest="genofile2",help="extra input file", metavar="FILE")
	parser.add_option("-r",dest="reportfile",help="problems report", metavar="FILE")
	parser.add_option("-s",dest="reportped",help="pedigree report",metavar="FILE")
	parser.add_option("-p",dest="pedigree",help="Pedigree", metavar="FILE")
	parser.add_option("-m",dest="markers",help="Markers", metavar="FILE")
	parser.add_option("-a",dest="pedlims",help="Limits for pedcheck", metavar="N[,N]",default='5')
	parser.add_option("-v",dest="verbose",action="store_true",help="prints runtime info",default=False)
	parser.add_option("-w",dest="nowarning",action="store_true",help="does not stop at missing markers/animals",default=False)
	parser.add_option("-G",dest="galaxy",action="store_true",help="Script is being run from galaxy",default=False)
	(options,args) = parser.parse_args()
	t = time.time()
	if options.galaxy:
		if options.genofile == 'None': options.genofile = None
		if options.pedigree == 'None': options.pedigree = None
		if options.markers == 'None': options.markers = None
		if options.genofile2 == 'None': options.genofile2 = None
		if options.reportped == 'None': options.reportped = None
	if not options.genofile:
		parser.error('ERROR: Input file is required.')
		sys.exit(1)
	if not options.reportfile and not options.reportped:
		parser.error('ERROR: Report file is required.')
		sys.exit(1)
	if options.genofile2 and not options.pedigree:
		parser.error('ERROR: Pedigree required when using second genotype file.')
		sys.exit(1)
	pedcheck(options)
	print "Time pedcheck: %.3f" % (time.time()-t)
        
if __name__ == "__main__":
	#import cProfile
	#cProfile.run('main()')
	main()
