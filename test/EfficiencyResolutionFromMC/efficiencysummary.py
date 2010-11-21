from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/ourvsvbtfeff')
cname = 'c12'

################################################################################

acc_f = ROOT.TFile('plots/trigeffvsmassmctruth/Acceptance_totals.root')
acc21_f = ROOT.TFile('plots/trigeffvsmassmctruth_vbtf/Acceptance_totals.root')

acc = acc_f.Get(cname).FindObject('Graph')
acc21 = acc21_f.Get(cname).FindObject('Graph')

acc.SetLineColor(ROOT.kRed)
acc21.SetLineColor(ROOT.kBlue)
acc.SetTitle('')
acc.GetXaxis().SetTitle('vector boson pole mass (GeV)')
acc.GetYaxis().SetTitle('acceptance')
acc.GetYaxis().SetTitleOffset(1.2)
acc.GetYaxis().SetRangeUser(0.4, 0.96)
for i,h in enumerate([acc, acc21]):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
acc.Draw('APL')
acc21.Draw('PL same')

lg = ROOT.TLegend(0.38, 0.143, 0.86, 0.316)
lg.AddEntry(acc,   'both muons | #eta| < 2.4', 'LPE')
lg.AddEntry(acc21, 'both muons | #eta| < 2.1', 'LPE')
lg.Draw()

ps.save('Acceptance', log=False)

################################################################################

our_f = ROOT.TFile('plots/trigeffvsmassmctruth/RecoWrtAcc_totals.root')
vbtf_f = ROOT.TFile('plots/trigeffvsmassmctruth_vbtf21/RecoWrtAcc_totals.root')

our = our_f.Get(cname).FindObject('Graph')
vbtf = vbtf_f.Get(cname).FindObject('Graph')

ox,oy,vx,vy = ROOT.Double(),ROOT.Double(),ROOT.Double(),ROOT.Double()
ratio = our.Clone('ratio')
print '%10s%10s%10s%10s' % ('pole mass', 'our', 'vbtf', 'ratio')
for i in xrange(our.GetN()):
    our.GetPoint(i, ox, oy)
    vbtf.GetPoint(i, vx, vy)
    print '%10i%10.3f%10.3f%10.3f' % (int(ox), oy, vy, vy/oy)
    ratio.SetPoint(i, ox, vy/oy)

our.SetTitle('')
our.GetXaxis().SetTitle('vector boson pole mass (GeV)')
our.GetYaxis().SetTitle('trigger+reco+selection efficiency wrt acceptance')
our.GetYaxis().SetTitleOffset(1.2)
our.GetYaxis().SetRangeUser(0.815, 1)

our.SetLineColor(ROOT.kRed)
vbtf.SetLineColor(ROOT.kBlue)
ratio.SetLineColor(ROOT.kGreen+2)

for i,h in enumerate([our,vbtf,ratio]):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
our.Draw('APL')
vbtf.Draw('PL same')
ratio.Draw('PL same')

lg = ROOT.TLegend(0.25, 0.13, 0.86, 0.29)
lg.AddEntry(our,   'efficiency w/ Our selection (acc:   #left|#eta#right| < 2.4)', 'LPE')
lg.AddEntry(vbtf,  'efficiency w/ VBTF selection (acc:   #left|#eta#right| < 2.1)', 'LPE')
lg.AddEntry(ratio, 'ratio VBTF/Ours', 'LPE')
lg.Draw()

ps.save('RecoWrtAcc', log=False)

################################################################################

def get_z0(g):
    g.GetPoint(0, ox, oy)
    assert(abs(ox - 90.5) < 0.1)
    return float(oy)
our_z0 = get_z0(our)
vbtf_z0 = get_z0(vbtf)

our_ratio_to_z0 = our.Clone('our_ratio_to_z0')
vbtf_ratio_to_z0 = our.Clone('vbtf_ratio_to_z0')
for i in xrange(our.GetN()):
    our.GetPoint(i, ox, oy)
    vbtf.GetPoint(i, vx, vy)
    our_ratio_to_z0.SetPoint(i, ox, our_z0/oy)
    vbtf_ratio_to_z0.SetPoint(i, vx, vbtf_z0/vy)

our_ratio_to_z0.SetLineColor(ROOT.kRed)
vbtf_ratio_to_z0.SetLineColor(ROOT.kBlue)

for i,h in enumerate([our_ratio_to_z0, vbtf_ratio_to_z0]):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

our_ratio_to_z0.GetXaxis().SetTitle('vector boson pole mass (GeV)')
our_ratio_to_z0.GetYaxis().SetTitle("ratio efficiency(Z0)/efficiency(Z')")
our_ratio_to_z0.GetYaxis().SetTitleOffset(1.3)
our_ratio_to_z0.GetYaxis().SetRangeUser(0.87, 1.01)
our_ratio_to_z0.Draw('APL')
vbtf_ratio_to_z0.Draw('PL same')

lg = ROOT.TLegend(0.45, 0.68, 0.86, 0.85)
lg.AddEntry(our_ratio_to_z0,  'Our selection', 'LPE')
lg.AddEntry(vbtf_ratio_to_z0, 'VBTF selection', 'LPE')
lg.Draw()

ps.save('LowToHighMass', log=False)

################################################################################

our_f = ROOT.TFile('plots/trigeffvsmassmctruth/TotalReco_totals.root')
vbtf_f = ROOT.TFile('plots/trigeffvsmassmctruth_vbtf/TotalReco_totals.root')
vbtf_acc_times_eff_at_z0 = get_z0(vbtf_f.Get(cname).FindObject('Graph'))

our_acc_times_eff = our_f.Get(cname).FindObject('Graph')
our_acc_times_eff.SetTitle('')
our_acc_times_eff.GetXaxis().SetTitle('vector boson pole mass (GeV)')
our_acc_times_eff.GetYaxis().SetTitle('acceptance times trigger+reco+selection efficiency')
our_acc_times_eff.GetYaxis().SetTitleOffset(1.2)
our_acc_times_eff.GetYaxis().SetRangeUser(0.4, 0.96)
our_acc_times_eff.SetMarkerStyle(20)
our_acc_times_eff.SetMarkerSize(0.9)

ps.c.cd()
our_acc_times_eff.Draw('APstats')

fcn = ROOT.TF1('fcn', '[0] + [1]/(x + [2])', 0, 2000)
fcn.SetParNames("a", "b", "c");
ROOT.gStyle.SetOptFit(11)
our_acc_times_eff.Fit(fcn, 'LVR')

tl = ROOT.TLatex(0.56, 0.36, 'a + #frac{b}{mass + c}')
tl.SetNDC()
tl.Draw()

ps.c.Update()
s = our_acc_times_eff.GetListOfFunctions().FindObject('stats')
s.SetFitFormat('5.3g')
s.SetX1NDC(0.51)
s.SetY1NDC(0.15)
s.SetX2NDC(0.85)
s.SetY2NDC(0.31)

ps.save('AcceptanceTimesEfficiency', log=False)
