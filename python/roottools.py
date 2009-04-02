#!/usr/bin/env python

import sys
sys.argv.append('-b') # start ROOT in batch mode
from ROOT import TH1F
sys.argv.remove('-b') # and don't mess up sys.argv

def get_bin_content_error(hist, value):
    """For the given histogram, find the bin corresponding to the
    value and return its contents and associated
    error. Multi-dimensional histograms are supported; value may be a
    tuple in those cases."""

    if type(value) != type(()):
        value = (value,)
    bin = hist.FindBin(*value)
    return (hist.GetBinContent(bin), hist.GetBinError(bin))

def get_integral(hist, xlo, xhi):
    """For the given histogram, return the integral of the bins
    corresponding to the values xlo to xhi along with its error."""
    
    binlo = hist.FindBin(xlo)
    binhi = hist.FindBin(xhi)
    integral = hist.Integral(binlo, binhi)
    wsq = 0
    for i in xrange(binlo, binhi+1):
        wsq += hist.GetBinError(i)**2
    return integral, wsq**0.5

def get_hist_stats(hist, fcnname='gaus'):
    """For the given histogram, return a five-tuple of the number of
    entries, the underflow and overflow counts, the fitted sigma
    (using the function specified by fcnname, which must be an
    already-made TF1 whose parameter(2) is the value used), and the
    RMS."""
    
    if hist.GetFunction(fcnname) is None:
        hist.Fit(fcnname,'Q')
    
    entries = hist.GetEntries()
    under = hist.GetBinContent(0)
    over = hist.GetBinContent(hist.GetNbinsX()+1)
    sigma = hist.GetFunction(fcnname).GetParameter(2)
    rms = hist.GetRMS()
    
    return entries, under, over, sigma, rms

def make_rms_hist(prof, cache={}):
    """Takes an input TProfile and produces a histogram whose bin contents are
    the RMS of the bins of the profile. Caches the histogram so that it doesn't
    get deleted by python before it gets finalized onto a TCanvas."""
    
    nbins = prof.GetNbinsX()
    name = 'RMS' + prof.GetName()
    # Play nice with same-name histograms that were OK because they
    # were originally in different directories.
    while cache.has_key(name):
        name += '1'

    h = TH1F(name, 'RMS ' + prof.GetTitle(), nbins, prof.GetXaxis().GetXmin(), prof.GetXaxis().GetXmax())
    for ibin in xrange(1, nbins+1):
        f_bin = float(prof.GetBinContent(ibin))
        ent_bin = float(prof.GetBinEntries(ibin))

        if f_bin > 0:
            f_bin = f_bin**0.5
        else:
            f_bin = 0
            
        if ent_bin > 0:
            err_bin = f_bin/(2*ent_bin)**0.5
        else:
            err_bin = 0

        h.SetBinContent(ibin, f_bin)
        h.SetBinError(ibin, err_bin)
        
    cache[name] = h
    return h

__all__ = ['get_bin_content_error', 'get_integral', 'get_hist_stats', 'make_rms_hist']

