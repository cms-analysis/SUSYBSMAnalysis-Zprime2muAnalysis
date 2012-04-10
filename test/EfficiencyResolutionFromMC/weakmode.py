import os
from math import pi, sin
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetOptStat(112210)

ROOT.gRandom = ROOT.TRandom3(12191982)

ps = plot_saver('plots/weakmode')

def doit(blammo, M, sigma, draw_gaus=False):
    limits = 0, M*2
    nbins = int((limits[1] - limits[0])/5)

    formula = '1./%(sigma)f/sqrt(2*3.14159265)*exp(-(x-%(M)f)*(x-%(M)f)/2/%(sigma)f/%(sigma)f)' % locals()
    f = ROOT.TF1('f', formula, *limits)
    h = ROOT.TH1F('h', '', nbins, *limits)
    hweakmode = ROOT.TH1F('hweakmode', '', nbins, *limits)

    for i in xrange(1000000):
        m = f.GetRandom()
        h.Fill(m)
        phi = ROOT.gRandom.Rndm()*2*pi - pi
        mweakmode = 1./(1./m + blammo/92/92 * sin(phi))
        hweakmode.Fill(mweakmode)

    h.SetTitle('M = %.f, pure resolution %.1f, no weak-mode effect' % (M, sigma))
    hweakmode.SetTitle('M = %.f, pure resolution %.1f, weak-mode effect of %.2f' % (M, sigma, blammo))

    for hh in (h, hweakmode):
        hh.GetYaxis().SetTitle('entries/5 GeV')
        hh.GetXaxis().SetTitle('mass (GeV)')

    if draw_gaus:
        h.Draw()
        ps.save('gaus_%i_%i_%i' % (blammo*100, M, sigma))
    hweakmode.Draw()
    ps.save('gaus_weakmoded_%i_%i_%i' % (blammo*100, M, sigma))
    return hweakmode

if False:
    doit(0.57, 1000, 0.122  * 1000)
    doit(0.57, 1000, 0.0825 * 1000)
    doit(0.73, 1000, 0.0825 * 1000)
    doit(0.73, 1000, 0.071  * 1000)
    doit(0.73,  500, 0.048  *  500)

print '%8s%8s%8s%8s' % ('M', 'm-M', 's/M', 'swm/M')
for which, blammo in [('tkonly', 0.56), ('tunep', 0.73)]:
    print which
    fn = 'plots/resfromzp/%s.root' % which
    f = ROOT.TFile(fn)
    pure_res = f.Get('c0').FindObject('Graph')
    weakmoded_res = pure_res.Clone('weakmoded_%s_res' % which)
    x,y = ROOT.Double(), ROOT.Double()
    for i in xrange(pure_res.GetN()):
        pure_res.GetPoint(i,x,y)
        M = float(x)
        sigma = float(y)*M
        hwm = doit(blammo, M, sigma)
        print '%8.f%8.4f%8.4f%8.4f' % (M, hwm.GetMean() - M, sigma/M, hwm.GetRMS()/M)
        weakmoded_res.SetPoint(i, M, hwm.GetRMS()/M)
        weakmoded_res.SetPointError(i, pure_res.GetErrorX(i), hwm.GetRMSError()/M)
        del hwm
    ps.c.cd()
    weakmoded_res.Draw('AP')
    ps.save('weakmoded_' + which)
