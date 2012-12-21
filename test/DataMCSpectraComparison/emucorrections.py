# (py emucorrections.py >! plots/out.emucorrections) && tlp plots/emucorrections

from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)

titles_on_individual = False
if titles_on_individual:
    ROOT.gStyle.SetTitleX(0.25)
    ROOT.gStyle.SetTitleY(0.40)

ps = plot_saver('plots/emucorrections', size=(600,600), pdf=True, pdf_log=True)

bins = array('d', [120, 200, 400, 600])

total_mumu = ROOT.TH1F('mumu', '', len(bins)-1, bins)
total_emu  = ROOT.TH1F('emu',  '', len(bins)-1, bins)

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import ttbar_powheg, ww, tW, tbarW, ztautau, wz, zz
samples = [ttbar_powheg, ww, tW, tbarW, ztautau, wz, zz]
cut = 'OurNew'

def draw(sample):
    for h,c in [(sample.emu, ROOT.kRed), (sample.mumu, ROOT.kBlue)]:
        #h.SetFillColor(c)
        h.SetLineWidth(2)
        h.SetLineColor(c)
        h.SetStats(0)
        title = '%s;dilepton mass (GeV);reconstructed events/bin' % (sample.nice_name if titles_on_individual else '')
        h.SetTitle(title)
        h.GetYaxis().SetTitleOffset(1.1)
        h.GetYaxis().SetLabelSize(0.025)
        h.GetXaxis().SetRangeUser(100, 620)

    h1,h2 = sort_histogram_pair(sample.mumu, sample.emu)
    h1.Draw('e')
    h2.Draw('e same')

    lg = ROOT.TLegend(0.67, 0.79, 0.90, 0.93)
    lg.AddEntry(sample.emu,  '#mu^{+}e^{-}/#mu^{-}e^{+}', 'LE')
    lg.AddEntry(sample.mumu, '#mu^{+}#mu^{-}', 'LE')
    lg.Draw()
    
    ps.save('%s_mumu_emu' % sample.name)

    sample.div = g = binomial_divide(sample.mumu, sample.emu, clopper_pearson_poisson_means, force_lt_1=False)
    g.SetTitle(';dilepton mass (GeV);correction factor n(#mu#mu)/n(e#mu)')
    g.GetYaxis().SetTitleOffset(1.2)
    g.GetYaxis().SetLabelSize(0.03)
    g.Draw('AP')
    ps.save('%s_div' % sample.name, log=False)

for sample in samples:
    f = ROOT.TFile('mc/ana_datamc_%s.root' % sample.name)

    h = f.Get(cut + 'MuonsPlusMuonsMinusHistos').Get('DileptonMass')
    sample.mumu = h.Rebin(len(bins)-1, '%s_mumu' % sample.name, bins)
    h = f.Get(cut + 'MuonsElectronsOppSignHistos').Get('DileptonMass')
    sample.emu = h.Rebin(len(bins)-1, '%s_emu' % sample.name, bins)

    sample.weight = sample.mumu.Clone('%s_weight' % sample.name)
    sample.weight.SetDirectory(0)
    sample.weight.Scale(sample.partial_weight)
    
    total_mumu.Add(sample.mumu, sample.partial_weight)
    total_emu .Add(sample.emu,  sample.partial_weight)

    draw(sample)


# Normalize by the sum of the weights.
lg = ROOT.TLegend(0.73, 0.12, 0.90, 0.445)
ROOT.gStyle.SetPaintTextFormat('.1e')
sum_weights = samples[0].weight.Clone('sum_weights')
for sample in samples[1:]:
    sum_weights.Add(sample.weight)
first = True
for sample in samples:
    sample.weight.Divide(sum_weights)
    sample.weight.SetTitle(';dilepton mass (GeV);relative weight in bin')
    sample.weight.SetStats(0)
    sample.weight.SetLineWidth(2)
    sample.weight.SetLineColor(sample.color)
    sample.weight.GetYaxis().SetLabelSize(0.025)
    sample.weight.GetYaxis().SetRangeUser(1e-3,1)
    for i in xrange(1, sample.weight.GetNbinsX()+1):
        sample.weight.SetBinError(i, 0)
    sample.weight.Draw('hist text00' if first else 'hist text00 same')
    first = False
    lg.AddEntry(sample.weight, sample.nice_name, 'LE')
lg.Draw()
ps.save('weights')


# Now draw the correction factors overlaid.
lg = ROOT.TLegend(0.60, 0.16, 0.87, 0.32)
first = True
for sample in [ttbar_powheg, ww]: # samples:
    g = sample.div
    g.SetLineWidth(2)
    g.SetLineColor(sample.color)
    g.SetMinimum(0.46)
    g.SetMaximum(0.72)
    g.Draw('AP' if first else 'P same')
    first = False
    lg.AddEntry(sample.weight, sample.nice_name, 'LE')

    n = g.GetN()
    x,y = ROOT.Double(), ROOT.Double()
    print sample.name
    for i in xrange(n):
        g.GetPoint(i, x, y)
        exl = g.GetErrorXlow(i)
        exh = g.GetErrorXhigh(i)
        eyl = g.GetErrorYlow(i)
        eyh = g.GetErrorYhigh(i)
        print '%20s%10.4f%20s' % ('%.f-%.f' % (x-exl, x+exh), y, '[%5.4f, %5.4f]' % (y-eyl, y+eyh))
    print

lg.Draw()
ps.save('all')
