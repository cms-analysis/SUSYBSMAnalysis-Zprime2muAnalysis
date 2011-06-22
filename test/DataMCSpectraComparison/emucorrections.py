from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/emucorrections')

bins = array('d', [120, 200, 300, 400, 600])

total_mumu = ROOT.TH1F('mumu', '', len(bins)-1, bins)
total_emu  = ROOT.TH1F('emu',  '', len(bins)-1, bins)

from samples import *
samples = [ttbar, ww, singletop_tW, ztautau, wz, zz]
cut = 'OurNew'

def draw(sample):
    for h,c in [(sample.emu, 9), (sample.mumu, 46)]:
        #h.SetFillColor(c)
        h.SetLineWidth(2)
        h.SetLineColor(c)
        h.SetStats(0)
        h.SetTitle('%s;dilepton mass (GeV);arb. units' % sample.name)

    h1,h2 = sort_histogram_pair(sample.mumu, sample.emu)
    h1.Draw('e')
    h2.Draw('e same')
    ps.save('%s_mumu_emu' % sample.name)

    sample.div = g = binomial_divide(sample.mumu, sample.emu, clopper_pearson_poisson_means, force_lt_1=False)
    g.SetTitle(';dilepton mass (GeV);correction factor n(#mu#mu)/n(e#mu)')
    g.Draw('AP')
    ps.save('%s_div' % sample.name, log=False)

for sample in samples:
    f = ROOT.TFile('/uscms/home/tucker/nobackup/zp2mu_ana_datamc_mc/latest/ana_datamc_%s.root' % sample.name)

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

#draw('total', mumu, emu, simple_div=True)

lg = ROOT.TLegend(0.62, 0.13, 0.81, 0.33)

# Normalize by the sum of the weights.
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
        
first = True
for sample in samples:
    g = sample.div
    g.SetLineWidth(2)
    g.SetLineColor(sample.color)
    g.Draw('AP' if first else 'P same')
    first = False

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

