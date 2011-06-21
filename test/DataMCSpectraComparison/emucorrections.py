from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/emucorrections')

bins = array('d', range(0, 2000, 50))
bins = array('d', [0, 60, 120, 200, 400, 600, 1000, 3000])

mumu = ROOT.TH1F('mumu', '', len(bins)-1, bins)
emu  = ROOT.TH1F('emu',  '', len(bins)-1, bins)

from samples import *
cut = 'OurNew'

def draw(name, this_mumu, this_emu):
    for h,c in [(this_emu, 9), (this_mumu, 46)]:
        h.SetFillColor(c)
        h.SetLineColor(c)
        h.SetStats(0)
        h.SetTitle('%s;dilepton mass (GeV);events/50 GeV' % sample.name)

    h1,h2 = sort_histogram_pair(this_mumu, this_emu)
    h1.Draw('hist')
    h2.Draw('hist same')
    ps.save('%s_mumu_emu' % name)

    h = this_mumu.Clone('h')
    h.Divide(this_emu)
    h.GetXaxis().SetRangeUser(120,2000)
    h.Draw('hist')
    ps.save('%s_div' % name, log=False)

for sample in [ztautau, zz, wz, ww, singletop_tW, ttbar]:
    f = ROOT.TFile('/uscms/home/tucker/nobackup/zp2mu_ana_datamc_mc/latest/ana_datamc_%s.root' % sample.name)

    this_mumu = f.Get(cut + 'MuonsPlusMuonsMinusHistos')  .Get('DileptonMass')
    this_emu  = f.Get(cut + 'MuonsElectronsOppSignHistos').Get('DileptonMass')

    for h,n in [(this_mumu, 'this_mumu'), (this_emu, 'this_emu')]:
        z = h.Rebin(len(bins)-1, n, bins)
        exec '%s = z' % n

    mumu.Add(this_mumu, sample.partial_weight)
    emu .Add(this_emu,  sample.partial_weight)

    draw(sample.name, this_mumu, this_emu)

draw('total', mumu, emu)
