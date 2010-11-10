import sys
sys.argv.append('-b')
import ROOT

mu = PDG_MASS_Z = 91.1876
gamma = PDG_WIDTH_Z = 2.4952
mu2 = mu*mu
gamma2 = gamma*gamma

def ZPole(x, par):
    mass = x[0]
    mass2 = mass*mass
    mass4 = mass2*mass2

    A, sigma, theta, B, C, kappa = par

    # first part
    poleval = A*ROOT.TMath.Voigt(mass-mu,gamma,sigma)
    pdfTerm = ROOT.TMath.Exp(-1.*theta*mass)

    # interference
    deltam2 = mass2 - mu2
    deltam22 = deltam2*deltam2
    interference = B*mu*deltam22 / (deltam22 + mass4*gamma2/mu2)
    #(mass2-mu2)*(mass2-mu2)/((mass2-mu2)*(mass2-mu2)+mass2*mass2*gamma*gamma/mu2)

    # final expo
    outexpo = C/mass2*ROOT.TMath.Exp(-kappa*mass)
    return (poleval + interference)*pdfTerm + outexpo

# tail fit
def Tail(x, par):
    mass = x[0]
    A, alpha, kappa = par
    return ROOT.TMath.Exp(A+-1.*alpha*ROOT.TMath.Power(mass,kappa))

# mix the two
def MixFuncEx(x, par, t):
    mass = x[0]

    tailA, tailAlpha, tailKappa, sigA, sigSigma, sigTheta, sigB, sigC, sigKappa, turnonWidth, scale = par

    erfval = 0.5 + 0.5*ROOT.TMath.Erf(turnonWidth * (mass-120))

    val1a = ZPole(x, (sigA,sigSigma,sigTheta,sigB,sigC,sigKappa))
    val1b = Tail(x, (tailA,tailAlpha,tailKappa))

    val1 = val1a * (1.-erfval)
    val2 = val1b * erfval

    if t == 1:
        return scale * val1
    elif t == 2:
        return scale * val2
    elif t == 3:
        return scale * val1a
    elif t == 4:
        return scale * val1b
    else:
        return scale*(val1+val2)

def MixFunc(x, par):
    return MixFuncEx(x, par, 0)
def MixFuncZOnly(x, par):
    return MixFuncEx(x, par, 1)
def MixFuncDYOnly(x, par):
    return MixFuncEx(x, par, 2)
def MixFuncZOnlyNoErf(x, par):
    return MixFuncEx(x, par, 3)
def MixFuncDYOnlyNoErf(x, par):
    return MixFuncEx(x, par, 4)
def FracZOnly(x, par):
    return MixFuncZOnly(x,par)/MixFunc(x,par)
def FracDYOnly(x, par):
    return MixFuncDYOnly(x,par)/MixFunc(x,par)

if __name__ == '__main__':
    pars = '''1  p0           2.44298e+01     fixed
   2  p1           6.99626e+00     fixed
   3  p2           2.03121e-01     fixed
   4  p3           9.45741e+06     fixed
   5  p4           2.63045e+00     fixed
   6  p5           1.80690e-02     fixed
   7  p6           5.12014e+01     fixed
   8  p7           0.00000e+00     fixed
   9  p8           0.00000e+00     fixed
   10  p9           2.19753e-02     fixed
   11  p10          5.63909e-03   1.69470e-03'''.split('\n')
    pars = [float(x.split()[2]) for x in pars if x]
    tailA, tailAlpha, tailKappa, sigA, sigSigma, sigTheta, sigB, sigC, sigKappa, turnonWidth, scale = pars

    leg = ROOT.TLegend(0.6, 0.6, 1.0, 1.0)
    for i,x in enumerate(['MixFunc', 'MixFuncZOnly', 'MixFuncDYOnly', 'MixFuncZOnlyNoErf', 'MixFuncDYOnlyNoErf', 'FracZOnly', 'FracDYOnly']):
        f = ROOT.TF1(x, eval(x), 120, 500, 11)
        f.SetLineColor(i+1)
        leg.AddEntry(f, x)
        for i,y in enumerate(pars):
            f.FixParameter(i, y)
        exec 'f%s = f' % x

    fMixFunc.Draw()
    fMixFuncZOnly.Draw('same')
    fMixFuncDYOnly.Draw('same')
    fMixFuncZOnlyNoErf.Draw('same')
    fMixFuncDYOnlyNoErf.Draw('same')

    leg.Draw()
    ROOT.c1.SetLogy(1)
    ROOT.c1.SaveAs('asdf/zlineshape.png')
    ROOT.c1.SaveAs('asdf/zlineshape.root')

    fFracZOnly.Draw()
    fFracDYOnly.Draw('same')
    ROOT.c1.SaveAs('asdf/zlineshapefrac.png')
    ROOT.c1.SaveAs('asdf/zlineshapefrac.root')
            
