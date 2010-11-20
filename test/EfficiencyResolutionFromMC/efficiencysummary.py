from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/ourvsvbtfeff')

our_f = ROOT.TFile('plots/trigeffvsmassmctruth/RecoWrtAcc_totals.root')
vbtf_f = ROOT.TFile('plots/trigeffvsmassmctruth_vbtf/RecoWrtAcc_totals.root')
vbtf21_f = ROOT.TFile('plots/trigeffvsmassmctruth_vbtf21/RecoWrtAcc_totals.root')

our = our_f.Get('c12').FindObject('Graph')
vbtf = vbtf_f.Get('c12').FindObject('Graph')
vbtf21 = vbtf21_f.Get('c12').FindObject('Graph')

our.SetLineColor(ROOT.kRed)
vbtf.SetLineColor(ROOT.kBlue)
vbtf21.SetLineColor(ROOT.kBlack)
our.SetTitle('')
our.GetXaxis().SetTitle('vector boson pole mass (GeV)')
our.GetYaxis().SetTitle('trigger+reco+selection efficiency wrt acceptance')
our.GetYaxis().SetTitleOffset(1.2)
our.GetYaxis().SetRangeUser(0.65, 0.96)

ox,oy = ROOT.Double(),ROOT.Double()
vx,vy = ROOT.Double(),ROOT.Double()
v21x,v21y = ROOT.Double(),ROOT.Double()
ratio = our.Clone('ratio')
ratio21 = our.Clone('ratio21')
print '%10s%10s%10s%10s%10s%10s' % ('pole mass', 'our', 'vbtf', 'vbtf21', 'ratio', 'rat21')
for i in xrange(our.GetN()):
    our.GetPoint(i, ox, oy)
    vbtf.GetPoint(i, vx, vy)
    vbtf21.GetPoint(i, v21x, v21y)
    print '%10i%10.3f%10.3f%10.3f%10.3f%10.3f' % (int(ox), oy, vy, v21y, vy/oy, v21y/oy)
    ratio.SetPoint(i, ox, vy/oy)
    ratio21.SetPoint(i, ox, v21y/oy)
ratio.SetLineColor(ROOT.kGreen+2)
ratio21.SetLineColor(28)

for i,h in enumerate([our,vbtf,vbtf21,ratio,ratio21]):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(1.0)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
our.Draw('APL')
vbtf.Draw('PL same')
vbtf21.Draw('PL same')

lg = ROOT.TLegend(0.39, 0.13, 0.86, 0.38)
lg.AddEntry(our,  'Using Our selection', 'LPE')
lg.AddEntry(vbtf, 'Using VBTF selection', 'LPE')
lg.AddEntry(vbtf21, 'Using VBTF selection, acc. = |#eta| < 2.1', 'LPE')
lg.Draw()

ps.save('RecoWrtAcc', log=False)

ratio.GetYaxis().SetTitle('ratio VBTF/Ours')
ratio.GetYaxis().SetTitleOffset(1.3)
ratio.GetYaxis().SetRangeUser(0.78, 1)
ratio.Draw('APL')
ratio21.Draw('PL same')

lg = ROOT.TLegend(0.39, 0.13, 0.86, 0.38)
lg.AddEntry(ratio, 'Ratio VBTF/Ours', 'LPE')
lg.AddEntry(ratio21, 'Ratio VBTF(acc. = |  #eta| < 2.1)/Ours', 'LPE')
lg.Draw()

ps.save('RecoWrtAccRatios', log=False)
