#!/usr/bin/env python

import sys
from array import array

sys.argv.append('-b') # start ROOT in batch mode
from ROOT import TH1F, TVirtualFitter
sys.argv.remove('-b') # and don't mess up sys.argv

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

def get_hist_stats(hist, extended=False, fcnname='gaus'):
    """For the given histogram, return a five-tuple of the number of
    entries, the underflow and overflow counts, the fitted sigma
    (using the function specified by fcnname, which must be an
    already-made TF1 whose parameter(2) is the value used), and the
    RMS.
    """
    
    if hist.GetFunction(fcnname) is None:
        hist.Fit(fcnname,'Q')
    
    entries = hist.GetEntries()
    under = hist.GetBinContent(0)
    over = hist.GetBinContent(hist.GetNbinsX()+1)
    sigma = hist.GetFunction(fcnname).GetParameter(2)
    rms = hist.GetRMS()
    if extended:
        sigma = (sigma, hist.GetFunction(fcnname).GetParError(2))
        rms = (rms, hist.GetRMSError())
    
    return entries, under, over, sigma, rms

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
        new_hist = TH1F(name, title, len(bins)-1, bins)
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
        new_hist = TH1F(name, title, nbins, old_axis.GetXmin(), old_axis.GetXmax())
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

__all__ = ['get_bin_content_error', 'get_integral', 'get_hist_stats', 'make_rms_hist']

