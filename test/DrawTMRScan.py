#!/usr/bin/env python

# This script takes as input a bunch of ROOT files created by
# Zprime2muResolution, created by scanning over the cut value. The
# scanning is done by a script such as the following:
#
# #!/bin/tcsh
#
# foreach cut (`seq 0 0.2 10`)
#   setenv TMRCUT ${cut}
#   cmsRun testZprime2muResolution_cfg.py
#   mv zp2mu_histos.root tmrscanresults/${cut}.root
# end
#
# This also requires in testZprime2muResolution_cfg.py the four lines
# that obtain the cut value from the environment variable TMRCUT to be
# uncommented.
#
# 

# Which rec level to be looking at (TMR).
rec_level = 9
# The path to the ROOT files.
root_file_dir = 'tmrscanresults/'
# Output file for plots.
output_file = 'TMRscan.ps'
# Range of cut values to look at on the plots.
cut_range = (0, 20)
# Proposed cut value (to be tweaked after looking at the plots).
proposed_cut = 3.5

import sys
import glob
from array import array
from PSDrawer import PSDrawer
from ROOT import TFile, TGraph, TLine, kRed

# Make some arrays for TGraphs.
cut = array('d')
stats = {
    # Order of tuple is entries, underflows, overflows, sigmas, rms
    'LeptonInvPtRes':  (array('d'), array('d'), array('d'), array('d'), array('d')),
    'DileptonMassRes': (array('d'), array('d'), array('d'), array('d'), array('d'))
    }

# Loop over every root file and extract the results of the fits.
for fn in glob.glob('%s*.root' % root_file_dir):
    c = float(fn.replace(root_file_dir, '').replace('.root', ''))
    if c < cut_range[0] or c > cut_range[1]: continue
    cut.append(c)

    f = TFile(fn)

    for base_name, (entries, under, over, sigma, rms) in stats.items():
        h = f.Zprime2muResolution.Get('%s%i' % (base_name, rec_level))
        h.Fit('gaus','Q')

        entries.append(h.GetEntries())
        under.append(h.GetBinContent(0))
        over.append(h.GetBinContent(h.GetNbinsX()+1))

        sigma.append(h.GetFunction('gaus').GetParameter(2))
        rms.append(h.GetRMS())

    f.Close()

# Draw the scatterplots.
psd = PSDrawer(output_file)

graphs = []
names = ['Entries', 'Underflows', 'Overflows', 'Sigma (fraction)', 'RMS (fraction)']
for base_name, arrays in stats.items():
    pad = psd.new_page(base_name, (2,3))
    for i, a in enumerate(arrays):
        pad.cd(i+1)
        
        g = TGraph(len(cut), cut, a)
        graphs.append(g) # Need to keep the graph object around.
        
        g.SetMarkerSize(0.6)
        g.GetXaxis().SetTitle('cut')
        g.SetTitle(names[i])
        
        g.Draw('AP')

# Draw a scatterplot of the track probabilities, with the proposed cut
# value superimposed.
f = TFile('%s%.1f.root' % (root_file_dir, proposed_cut))
h = f.Zprime2muHistos.Get('TeVMuonLnProb56')
if h is None:
    raise RuntimeError, 'input file does not contain the histogram!'

pad = psd.new_page('Choice of cut value', (1,2))
pad.cd(1)

h.DrawCopy()

pad.cd(2)

h.GetXaxis().SetRangeUser(0, 10)
h.GetYaxis().SetRangeUser(0, 20)
h.Draw()

maxTK = 10

l = TLine(0, proposed_cut, maxTK, proposed_cut + maxTK)
l.SetLineColor(kRed)
l.Draw()

psd.close()
    
########################################################################

#import glob
#from mymisc import depickle
#
#for f in sorted(glob.glob('tmrscanresults/*.pydata')):
#    d = depickle(f)
#    cut = float(f.replace('tmrscanresults/','').replace('.pydata',''))
#    for t in ('LeptonInvPtRes','DileptonMassRes'):
#        print '%s %f %i %i %i %f %f' % ((t, cut) + d[t][9])

########################################################################

#levels = ['GR','TK','FS','PR','OP','TR']
#stats = {}
#
#def dump_sigmas(base_name):
#    stats[base_name] = {}
#    print '%s:' % base_name
#    print '%8s %8s %8s %8s %8s %8s' % ('level', 'entries', 'under', 'over', 'sigma', 'rms')
#    print '-'*53
#    for rec in xrange(4,10):
#        level = levels[rec-4]
#        h = histos.Get('%s%i' % (base_name, rec))
#        h.Fit('gaus','Q')
#        sigma = h.GetFunction('gaus').GetParameter(2)
#        rms = h.GetRMS()
#        entries = h.GetEntries()
#        under = h.GetBinContent(0)
#        over = h.GetBinContent(h.GetNbinsX()+1)
#        stats[base_name][rec] = (entries, under, over, sigma, rms)
#        print '%8s  %8i %8i %8i %8.5f %8.5f' % (level, entries, under, over, sigma, rms)
#    print
#
#dump_sigmas('LeptonInvPtRes')
#dump_sigmas('DileptonMassRes')
#
#import cPickle as pickle
#open(rootFile.replace('.root','.pydata'), 'wb').write(pickle.dumps(stats))

