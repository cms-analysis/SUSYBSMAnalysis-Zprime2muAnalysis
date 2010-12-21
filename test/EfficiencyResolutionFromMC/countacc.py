import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

max_eta = 2.1
min_pt = 20

print '%15s%9s%9s%9s%7s%7s%9s%7s%7s' % ('which', 'tot' , 'num', 'den', 'acc', 'err', 'denr', 'accr', 'errr')
for fn in sys.argv[1:]:
    f = ROOT.TFile(fn)
    num = 0
    denr = 0
    den = 0
    t = f.GenNtuple.Get('t')
    tot = t.GetEntries()
    for j,t in ttree_iterator(t):
        if j > 150000:
            break
        #print t.res_mass, t.dil_mass, t.lep_eta[0], t.lep_eta[1], t.lep_pt[0], t.lep_pt[1]
        if t.res_mass > 60 and t.res_mass < 120:
            denr += 1
        dmok = t.dil_mass > 60 and t.dil_mass < 120
        if dmok:
            den += 1
        if dmok and abs(t.lep_eta[0]) < max_eta and abs(t.lep_eta[1]) < max_eta and t.lep_pt[0] > min_pt and t.lep_pt[1] > min_pt:
                num += 1
    acc,l,h = clopper_pearson(num,den)
    err = acc - l
    accr,l,h = clopper_pearson(num,denr)
    errr = accr - l
    print '%15s%9i%9i%9i%7.4f%7.4f%9i%7.4f%7.4f' % (fn.replace('.root', '').replace('genntuple_', ''), tot, num, den, acc, err, denr, accr, errr) 
    
