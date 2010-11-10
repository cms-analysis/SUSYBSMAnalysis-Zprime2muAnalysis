from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, plot_saver, set_zp2mu_style
set_zp2mu_style()

ps = plot_saver('plots/emucorrections')

from samples import samples
for sample in samples:
    f = ROOT.TFile('ana_datamc/ana_datamc_%s.root' % sample.name)

    cut = 'VBTF'
    mumu = f.Get(cut + 'MuonsPlusMuonsMinusHistos').Get('DileptonMass')
    emu = f.Get(cut + 'MuonsElectronsOppSignHistos').Get('DileptonMass')

    #for h in (mumu,emu):
    #    h.Rebin(5)

    bins = array('d', range(0, 201, 40) + [2000])
    for h,n in [(mumu, 'mumu'), (emu, 'emu')]:
        z = h.Rebin(len(bins)-1, n, bins)
        exec '%s = z' % n

    mumu.Draw('hist')
    ps.save('%s_mumu' % sample.name)
    emu.Draw('hist')
    ps.save('%s_emu' % sample.name)

    h = mumu.Clone('h')
    h.Divide(emu)
    h.Draw('hist')
    ps.save('%s_div' % sample.name)
