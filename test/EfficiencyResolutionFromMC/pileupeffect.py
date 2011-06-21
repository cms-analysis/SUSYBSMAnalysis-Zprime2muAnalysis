from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/pileupeffect')
cname = 'c0'

across_intime = ['0509', '61109', '122609']
across_late = ['0509', '051026']

def intime_files(pattern):
    return [ROOT.TFile(pattern % x) for x in across_intime]      

def late_files(pattern):
    return [ROOT.TFile(pattern % x) for x in across_late]      

def get_graphs(files):
    return [f.Get(cname).FindObject('Graph') for f in files]

################################################################################

few_late_f, lot_late_f = fs = late_files('plots/trigeffvsmassmctruth/%s/TotalReco_totals.root')
few_late, lot_late = gs = get_graphs(fs)

few_late.SetLineColor(ROOT.kRed)
lot_late.SetLineColor(ROOT.kBlue)
few_late.SetTitle('')
few_late.GetXaxis().SetTitle('vector boson pole mass (GeV)')
few_late.GetYaxis().SetTitle('efficiency')
few_late.GetYaxis().SetTitleOffset(1.2)
few_late.GetXaxis().SetRangeUser(600, 1900)
few_late.GetYaxis().SetRangeUser(0.78, 0.88)
for i,h in enumerate(gs):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
few_late.Draw('APL')
lot_late.Draw('PL same')

lg = ROOT.TLegend(0.38, 0.143, 0.86, 0.316)
lg.AddEntry(few_late, '0-9 1BX late PU')
lg.AddEntry(lot_late, '10-25 1BX late PU')
lg.Draw()

ps.save('LateTotalReco', log=False)

################################################################################

few_late_f, lot_late_f = fs = late_files('plots/trigeffvsmassmctruth/%s/RecoWrtAcc_totals.root')
few_late, lot_late = gs = get_graphs(fs)

few_late.SetLineColor(ROOT.kRed)
lot_late.SetLineColor(ROOT.kBlue)
few_late.SetTitle('')
few_late.GetXaxis().SetTitle('vector boson pole mass (GeV)')
few_late.GetYaxis().SetTitle('efficiency')
few_late.GetYaxis().SetTitleOffset(1.2)
few_late.GetYaxis().SetRangeUser(0.85, 0.95)
for i,h in enumerate(gs):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
few_late.Draw('APL')
lot_late.Draw('PL same')

lg = ROOT.TLegend(0.38, 0.143, 0.86, 0.316)
lg.AddEntry(few_late, '0-9 1BX late PU')
lg.AddEntry(lot_late, '10-25 1BX late PU')
lg.Draw()

ps.save('LateRecoWrtAcc', log=False)

################################################################################


fs = intime_files('plots/trigeffvsmassmctruth/%s/TotalReco_totals.root')
few_intime, mid_intime, lot_intime = gs = get_graphs(fs)

few_intime.SetLineColor(ROOT.kRed)
mid_intime.SetLineColor(ROOT.kGreen+2)
lot_intime.SetLineColor(ROOT.kBlue)
few_intime.SetTitle('')
few_intime.GetXaxis().SetTitle('vector boson pole mass (GeV)')
few_intime.GetYaxis().SetTitle('efficiency')
few_intime.GetYaxis().SetTitleOffset(1.2)
few_intime.GetXaxis().SetRangeUser(600, 1900)
few_intime.GetYaxis().SetRangeUser(0.78, 0.88)
for i,h in enumerate(gs):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
few_intime.Draw('APL')
mid_intime.Draw('PL same')
lot_intime.Draw('PL same')

lg = ROOT.TLegend(0.38, 0.143, 0.86, 0.316)
lg.AddEntry(few_intime, '0-5 1BX intime PU')
lg.AddEntry(mid_intime, '6-11 1BX intime PU')
lg.AddEntry(lot_intime, '12-25 1BX intime PU')
lg.Draw()

ps.save('IntimeTotalReco', log=False)

################################################################################

fs = intime_files('plots/trigeffvsmassmctruth/%s/RecoWrtAcc_totals.root')
few_intime, mid_intime, lot_intime = gs = get_graphs(fs)

few_intime.SetLineColor(ROOT.kRed)
mid_intime.SetLineColor(ROOT.kGreen+2)
lot_intime.SetLineColor(ROOT.kBlue)
few_intime.SetTitle('')
few_intime.GetXaxis().SetTitle('vector boson pole mass (GeV)')
few_intime.GetYaxis().SetTitle('efficiency')
few_intime.GetYaxis().SetTitleOffset(1.2)
few_intime.GetXaxis().SetRangeUser(600, 1900)
few_intime.GetYaxis().SetRangeUser(0.85, 0.95)
for i,h in enumerate(gs):
    h.SetMarkerStyle(20+i)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(h.GetLineColor())

ps.c.cd()
few_intime.Draw('APL')
mid_intime.Draw('PL same')
lot_intime.Draw('PL same')

lg = ROOT.TLegend(0.38, 0.143, 0.86, 0.316)
lg.AddEntry(few_intime, '0-5 1BX intime PU')
lg.AddEntry(mid_intime, '6-11 1BX intime PU')
lg.AddEntry(lot_intime, '12-25 1BX intime PU')
lg.Draw()

ps.save('IntimeRecoWrtAcc', log=False)
