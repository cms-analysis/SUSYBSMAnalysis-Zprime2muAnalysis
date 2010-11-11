#!/usr/bin/env python

import sys, os
from array import array

sys.argv.append('-b')     # Start ROOT in batch mode;
import ROOT; ROOT.TCanvas # make sure libGui gets initialized while '-b' is specified;
sys.argv.remove('-b')     # and don't mess up sys.argv.

def apply_hist_commands(hist, hist_cmds=None):
    """With hist_cmds a list of n-tuples, where the first entry of the
    tuple is a function name, call the function name on the histogram
    passing the n-1 remaining arguments to the function.

    E.g. if hist_cmds is [('SetMarkerStyle', 5), ('SetStats', 0)],
    this is equivalent to calling

    hist.SetMarkerStyle(5)
    hist.SetStats(0)

    Default is a no-op (hist_cmds = None).
    """

    if hist_cmds is None: return
    for fn, args in hist_cmds:
        t = type(args)
        if t != type(()) or t != type([]):
            args = (args,)
        getattr(hist, fn)(*args)

def clopper_pearson(n_on, n_tot, alpha=1-0.6827, equal_tailed=True):
    if equal_tailed:
        alpha_min = alpha/2
    else:
        alpha_min = alpha

    lower = 0
    upper = 1

    if n_on > 0:
        lower = ROOT.Math.beta_quantile(alpha_min, n_on, n_tot - n_on + 1)
    if n_tot - n_on > 0:
        upper = ROOT.Math.beta_quantile_c(alpha_min, n_on + 1, n_tot - n_on)

    if n_on == 0 and n_tot == 0:
        return 0, lower, upper
    else:
        return float(n_on)/n_tot, lower, upper

def binomial_divide(h1, h2, confint=clopper_pearson):
    nbins = h1.GetNbinsX()
    xax = h1.GetXaxis()
    if h2.GetNbinsX() != nbins: # or xax2.GetBinLowEdge(1) != xax.GetBinLowEdge(1) or xax2.GetBinLowEdge(nbins) != xax.GetBinLowEdge(nbins):
        raise ValueError, 'incompatible histograms to divide'
    x = []
    y = []
    exl = []
    exh = []
    eyl = []
    eyh = []
    xax = h1.GetXaxis()
    for ibin in xrange(1, nbins+1):
        s,t = h1.GetBinContent(ibin), h2.GetBinContent(ibin)
        if t == 0:
            assert(s == 0)
            continue

        p_hat = float(s)/t
        if s > t:
            print 'warning: bin %i has p_hat > 1, in interval forcing p_hat = 1' % ibin
            s = t
        rat, a,b = confint(s,t)
        #print ibin, s, t, a, b

        _x  = xax.GetBinCenter(ibin)
        _xw = xax.GetBinWidth(ibin)/2
        
        x.append(_x)
        exl.append(_xw)
        exh.append(_xw)

        y.append(p_hat)
        eyl.append(p_hat - a)
        eyh.append(b - p_hat)
    eff = ROOT.TGraphAsymmErrors(len(x), *[array('d', obj) for obj in (x,y,exl,exh,eyl,eyh)])
    return eff

def core_gaussian(hist, factor, i=[0]):
    core_mean  = hist.GetMean()
    core_width = factor*hist.GetRMS()
    f = ROOT.TF1('core%i' % i[0], 'gaus', core_mean - core_width, core_mean + core_width)
    i[0] += 1
    return f

def fit_gaussian(hist, factor=None, draw=False, cache=[]):
    """Fit a Gaussian to the histogram, and return a dict with fitted
    parameters and errors. If factor is supplied, fit only to range in
    hist.mean +/- factor * hist.rms.
    """

    if draw:
        opt = 'qr'
    else:
        opt = 'qr0'

    if factor is not None:
        fcn = core_gaussian(hist, factor)
        cache.append(fcn)
        hist.Fit(fcn, opt)
    else:
        hist.Fit('gaus', opt)
        fcn = hist.GetFunction('gaus')
        
    return {
        'constant': (fcn.GetParameter(0), fcn.GetParError(0)),
        'mu':       (fcn.GetParameter(1), fcn.GetParError(1)),
        'sigma':    (fcn.GetParameter(2), fcn.GetParError(2))
        }
    
def get_bin_content_error(hist, value):
    """For the given histogram, find the bin corresponding to the
    value and return its contents and associated
    error. Multi-dimensional histograms are supported; value may be a
    tuple in those cases.
    """

    if type(value) != type(()):
        value = (value,)
    bin = hist.FindBin(*value)
    return (hist.GetBinContent(bin), hist.GetBinError(bin))

def get_integral(hist, xlo, xhi):
    """For the given histogram, return the integral of the bins
    corresponding to the values xlo to xhi along with its error.
    """
    
    binlo = hist.FindBin(xlo)
    binhi = hist.FindBin(xhi)
    integral = hist.Integral(binlo, binhi)
    wsq = 0
    for i in xrange(binlo, binhi+1):
        wsq += hist.GetBinError(i)**2
    return integral, wsq**0.5

def get_hist_stats(hist, factor=None):
    """For the given histogram, return a five-tuple of the number of
    entries, the underflow and overflow counts, the fitted sigma
    (using the function specified by fcnname, which must be an
    already-made ROOT.TF1 whose parameter(2) is the value used), and the
    RMS.
    """

    results = fit_gaussian(hist, factor)
    results.update({
        'entries': hist.GetEntries(),
        'under':   hist.GetBinContent(0),
        'over':    hist.GetBinContent(hist.GetNbinsX()+1),
        'mean':    (hist.GetMean(), hist.GetMeanError()),
        'rms':     (hist.GetRMS(), hist.GetRMSError())
        })
    return results

def make_rms_hist(prof, name='', bins=None, cache={}):
    """Takes an input TProfile and produces a histogram whose bin contents are
    the RMS of the bins of the profile. Caches the histogram so that it doesn't
    get deleted by python before it gets finalized onto a TCanvas.

    If bins is a list of bin lower edges + last bin high edge,
    rebinning is done before making the RMS histogram. Due to a bug in
    ROOT's TProfile in versions less than 5.22 (?), rebinning is done
    manually here.
    """
    
    nbins = prof.GetNbinsX()
    if name == '':
        name = 'RMS' + prof.GetName()
    title = 'RMS ' + prof.GetTitle()
    old_axis = prof.GetXaxis()

    # Play nice with same-name histograms that were OK because they
    # were originally in different directories.
    while cache.has_key(name):
        name += '1'

    # Format of contents list: [(new_bin, (new_bin_content, new_bin_error)), ...]
    contents = []
    
    if bins:
        if type(bins) == type([]):
            bins = array('f', bins)
        new_hist = ROOT.TH1F(name, title, len(bins)-1, bins)
        new_axis = new_hist.GetXaxis()
        new_bins = {}
        for old_bin in xrange(1, nbins+1):
            new_bin = new_axis.FindBin(old_axis.GetBinLowEdge(old_bin))
            if not new_bins.has_key(new_bin):
                new_bins[new_bin] = [0., 0.]
            N = prof.GetBinEntries(old_bin)
            new_bins[new_bin][0] += N*prof.GetBinContent(old_bin)
            new_bins[new_bin][1] += N
        for val in new_bins.values():
            if val[1] > 0:
                val[0] /= val[1]
        contents = new_bins.items()
    else:
        new_hist = ROOT.TH1F(name, title, nbins, old_axis.GetXmin(), old_axis.GetXmax())
        for old_bin in xrange(1, nbins+1):
            f_bin = float(prof.GetBinContent(old_bin))
            ent_bin = float(prof.GetBinEntries(old_bin))
            contents.append((old_bin, (f_bin, ent_bin)))

    for new_bin, (f_bin, ent_bin) in contents:
        if f_bin > 0:
            f_bin = f_bin**0.5
        else:
            f_bin = 0
            
        if ent_bin > 0:
            err_bin = f_bin/(2.*ent_bin)**0.5
        else:
            err_bin = 0

        new_hist.SetBinContent(new_bin, f_bin)
        new_hist.SetBinError(new_bin, err_bin)
        
    cache[name] = new_hist
    return new_hist

def move_above_into_bin(h,a):
    """Given the TH1 h, add the contents of the bins above the one
    corresponding to a into that bin, and zero the bins above."""
    assert(h.Class().GetName().startswith('TH1')) # i bet there's a better way to do this...
    b = h.FindBin(a)
    bc = h.GetBinContent(b)
    bcv = h.GetBinError(b)**2
    for nb in xrange(b+1, h.GetNbinsX()+2):
        bc += h.GetBinContent(nb)
        bcv += h.GetBinError(nb)**2
        h.SetBinContent(nb, 0)
        h.SetBinError(nb, 0)
    h.SetBinContent(b, bc)
    h.SetBinError(b, bcv**0.5)

def move_overflow_into_last_bin(h):
    """Given the TH1 h, Add the contents of the overflow bin into the
    last bin, and zero the overflow bin."""
    assert(h.Class().GetName().startswith('TH1')) # i bet there's a better way to do this...
    nb = h.GetNbinsX()
    h.SetBinContent(nb, h.GetBinContent(nb) + h.GetBinContent(nb+1))
    h.SetBinError(nb, (h.GetBinError(nb)**2 + h.GetBinError(nb+1)**2)**0.5)
    h.SetBinContent(nb+1, 0)
    h.SetBinError(nb+1, 0)

class plot_saver:
    i = 0
    
    def __init__(self, plot_dir=None, html=True):
        self.c = ROOT.TCanvas('c%i' % plot_saver.i, '', 820, 630)
        plot_saver.i += 1
        self.saved = []
        self.html = html
        self.set_plot_dir(plot_dir)

    def __del__(self):
        self.write_index()
        
    def write_index(self):
        if not self.saved or not self.html:
            return
        html = open(os.path.join(self.plot_dir, 'index.html'), 'wt')
        html.write('<html><body><pre>\n')
        for i, (fn, log, root) in enumerate(self.saved):
            bn = os.path.basename(fn)
            html.write('%10i ' % i)
            if log:
                html.write(' <a href="%s">log</a>' % os.path.basename(log))
            if root:
                html.write(' <a href="%s">root</a>' % os.path.basename(root))
            html.write('  <a href="%s">%s</a>' % (bn, bn))
            html.write('\n')
        html.write('<br><br>')
        for i, (fn, log, root) in enumerate(self.saved):
            bn = os.path.basename(fn)
            html.write('%s<br>\n' % bn.replace('.png', ''))
            if log:
                html.write('<img src="%s"><img src="%s"><br><br>\n' % (bn, os.path.basename(log)))
            else:
                html.write('<img src="%s"><br><br>\n' % bn)
        html.write('</pre></body></html>\n')
        
    def set_plot_dir(self, plot_dir):
        self.write_index()
        self.saved = []
        self.plot_dir = plot_dir
        if plot_dir is not None:
            os.system('mkdir -p %s' % self.plot_dir)

    def save(self, n, log=True, root=True):
        if self.plot_dir is None:
            raise ValueError('save called before plot_dir set!')
        self.c.SetLogy(0)
        fn = os.path.join(self.plot_dir, n + '.png')
        self.c.SaveAs(fn)
        if root:
            root = os.path.join(self.plot_dir, n + '.root')
            self.c.SaveAs(root)
        if log:
            self.c.SetLogy(1)
            log = os.path.join(self.plot_dir, n + '.log.png')
            self.c.SaveAs(log)
            self.c.SetLogy(0)
        self.saved.append((fn, log, root))

def rainbow_palette(num_colors=500):
    """Make a rainbow palette with the specified number of
    colors. Also call SetNumberContours so it actually gets used when
    drawing COLZ."""
    
    r = array('d', [0, 0, 0, 1, 1])
    g = array('d', [0, 1, 1, 1, 0])
    b = array('d', [1, 1, 0, 0, 0])
    stops = array('d', [float(i)/4 for i in xrange(5)])
    ROOT.TColor.CreateGradientColorTable(5, stops, r, g, b, num_colors)
    ROOT.gStyle.SetNumberContours(num_colors)

def real_hist_max(h, return_bin=False, user_range=None, use_error_bars=True):
    """Find the real maximum value of the histogram, taking into
    account the error bars and/or the specified range."""

    m_ibin = None
    m = 0

    if user_range is None:
        b1, b2 = 1, h.GetNbinsX() + 1
    else:
        b1, b2 = h.FindBin(user_range[0]), h.FindBin(user_range[1])+1
    
    for ibin in xrange(b1, b2):
        if use_error_bars:
            v = h.GetBinContent(ibin) + h.GetBinError(ibin)
        else:
            v = h.GetBinContent(ibin)
        if v > m:
            m = v
            m_ibin = ibin
    if return_bin:
        return m_ibin, m
    else:
        return m

def real_hist_min(h, return_bin=False, user_range=None):
    """Find the real minimum value of the histogram, ignoring empty
    bins, and taking into account the specified range."""

    m_ibin = None
    m = 99e99

    if user_range is None:
        b1, b2 = 1, h.GetNbinsX() + 1
    else:
        b1, b2 = h.FindBin(user_range[0]), h.FindBin(user_range[1])+1
    
    for ibin in xrange(b1, b2):
        v = h.GetBinContent(ibin)
        if v > 0 and v < m:
            m = v
            m_ibin = ibin
    if return_bin:
        return m_ibin, m
    else:
        return m

def set_zp2mu_style(date_pages=False):
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetFillColor(0)
    if date_pages:
        ROOT.gStyle.SetOptDate()
    ROOT.gStyle.SetOptStat(111111)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetMarkerSize(.1)
    ROOT.gStyle.SetMarkerStyle(8)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetPaperSize(ROOT.TStyle.kA4)
    ROOT.gStyle.SetStatW(0.25)
    ROOT.gStyle.SetStatFormat('6.4g')
    #ROOT.gStyle.SetTitleFont(52, 'XY')
    #ROOT.gStyle.SetLabelFont(52, 'XY')
    #ROOT.gStyle.SetStatFont(52)
    ROOT.gErrorIgnoreLevel = 1001 # Suppress TCanvas::SaveAs messages.

def ttree_iterator(tree, return_tree=True):
    for jentry in xrange(tree.GetEntriesFast()):
        if tree.LoadTree(jentry) < 0: break
        if tree.GetEntry(jentry) <= 0: continue
        if return_tree:
            yield jentry, tree
        else:
            yield jentry

__all__ = [
    'apply_hist_commands',
    'binomial_divide',
    'clopper_pearson',
    'core_gaussian',
    'fit_gaussian',
    'get_bin_content_error',
    'get_integral',
    'get_hist_stats',
    'make_rms_hist',
    'move_above_into_bin',
    'move_overflow_into_last_bin',
    'plot_saver',
    'rainbow_palette',
    'real_hist_max',
    'real_hist_min',
    'set_zp2mu_style',
    'ttree_iterator',
    'ROOT',
    ]
