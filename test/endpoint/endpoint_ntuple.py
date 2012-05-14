#!/usr/bin/env python

import ROOT
import os,sys
#from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

#f = ROOT.TFile('data/Run2012PlusDCSOnlyMuonsOnly/ana_datamc_data.root')
path = 'data/Run2012PlusDCSOnlyMuonsOnly/ana_datamc_data.root'
tmp_fn = 'micro_ntuple.temp1.txt'
branch_spec = 'vertex_m'
branch_spec = 'vertex_m:lep_pt'
branch_spec = 'lep_id:lep_pt:lep_eta:lep_phi:lep_tk_pt:lep_cocktail_pt'
branch_spec = 'lep_id:lep_tk_pt:lep_eta:lep_phi'
branch_spec = 'lep_id:lep_cocktail_pt:lep_eta:lep_phi'
cut = 'lep_pt'
cut1 = 'OurSelNew'
cut2 = 'OurSelNewNoSign'

f = ROOT.TFile(path)
t = f.SimpleNtupler.Get('t')

c1 = t.GetAlias(cut1)
c2 = t.GetAlias(cut2)
c3 = t.GetAlias('loose_new_0')
c4 = t.GetAlias('loose_new_1')

print c3
c3n = c3.replace('45','15')
c4 = c4.replace('45','15')
print c3n
print c3
t.SetAlias("loose_new_0",c3n)
t.SetAlias("loose_new_1",c4)

cendpoint = cut1+'& vertex_m>60 & lep_id[0]*lep_id[1]<0'

t.GetPlayer().SetScanRedirect(True)
t.GetPlayer().SetScanFileName(tmp_fn)
t.Scan(branch_spec, cendpoint)
t.GetPlayer().SetScanRedirect(False)
f.Close()

import math
iline =0

ftp2=open("arp.txt","wt")
for line in open(tmp_fn).readlines():

	iline +=1
	if iline <4: continue
	if (len(line)<6): continue

	if ' * ' not in line and 'Row' in line: continue
	line = line.split(' *')[2:]
	if (len(line)<4): continue

#	vals = [int(line[0])/(-13),line[1:]]
	vals = [int(line[0])/(-13*float(line[1])),line[2],line[3]]

	for v in vals: ftp2.write("%.5f\t"%float(v))
	ftp2.write('\n')

#	ftp2.write('\n'.join(vals))


#	x =1
#	math.copysign(x,-line[0])
#	spt = math.copysign(-1*line[1],line[0])
#	print line
#	print vals 
	

#	if iline >10: break

f = ROOT.TFile("tmp2.root","recreate")
t = ROOT.TTree('tree','tree')
t.ReadFile('arp.txt', 'k:eta:phi')
f.Write()
f.Close()
	

'''
lines = [line.split(' *')[1:] for line in open(tmp_fn).readlines() if ' * ' in line and 'Row' not in line]
lines = ['\t'.join(y.strip() for y in x) for x in lines]
open(tmp_fn, 'wt').write('\n'.join(lines))
'''


'''
#f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
f = ROOT.TFile("tmp.root","recreate")
t = ROOT.TTree('t','')
t.ReadFile(tmp_fn, branch_spec)
f.Write()
f.Close()

#os.remove(tmp_fn)
'''
