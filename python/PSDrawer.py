#!/usr/bin/env python

import os
import sys

sys.argv.append('-b') # start ROOT in batch mode
from ROOT import *
sys.argv.remove('-b') # and don't mess up sys.argv

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import rec_levels, rec_level_code
recLevelDict, recLevels = rec_levels()

from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import apply_hist_commands, make_rms_hist, set_zp2mu_style

class PSDrawer:
    GEN = recLevels.index('GN')
    REC_START = recLevels.index('GN') + 1
    TRIG_START = recLevels.index('L1')
    OFFLINE_START = recLevels.index('GR')
    COCKTAIL_START = recLevels.index('OP')
    MAX_LEVELS = len(recLevels)

    divs = {
        'all':      (2,5),
        'offline':  (2,3),
        'no_gen':   (2,5),
        'no_trig':  (2,4),
        'cocktail': (1,2)
        }

    levels = {
        'all':      range(MAX_LEVELS),
        'offline':  range(OFFLINE_START, MAX_LEVELS),
        'no_gen':   range(REC_START, MAX_LEVELS),
        'no_trig':  [GEN] + range(OFFLINE_START, MAX_LEVELS),
        'cocktail': range(COCKTAIL_START, MAX_LEVELS)
        }        
        
    def __init__(self, filename, date_pages=False, as_pdf=False):
        set_zp2mu_style(date_pages)
        self.canvas = TCanvas('c1','',0,0,500,640)
        self.ps = TPostScript(filename, 111)
        self.filename = filename
        self.page = 1    
        self.t = TText()
        self.t.SetTextFont(32)
        self.t.SetTextSize(0.025)
        
        self.as_pdf = as_pdf
        self.closed = False

    def __del__(self):
        self.close()

    def new_page(self, title, div=(1,2)):
        self.canvas.Update()
        self.ps.NewPage()
        self.canvas.Clear()
        self.canvas.cd(0)
        self.title = TPaveLabel(0.1, 0.94, 0.9, 0.98, title)
        self.title.SetFillColor(10)
        self.title.Draw()
        self.t.DrawText(0.9, 0.02, '- %d -' % self.page)
        self.pad = TPad('', '', 0.05, 0.05, 0.93, 0.93)
        self.pad.Draw()
        self.pad.Divide(*div)
        self.page += 1
        return self.pad

    def draw_if(self, histos, name, draw_opt=''):
        h = histos.Get(name)
        if h: h.Draw(draw_opt)
        return h

    def div_levels(self, page_type):
        div = self.divs[page_type]
        levels = self.levels[page_type]
        if div[0]*div[1] < len(levels):
            raise RuntimeError, 'not enough divisions (%i) for number of levels (%i)' % (div[0]*div[1], len(levels))
        return div, levels

    def rec_level_page(self, histos, page_type, histo_base_name, page_title, draw_opt='', log_scale=False, fit_gaus=False, hist_cmds=None, prof2rms=False, fit_other=None):
        div, levels = self.div_levels(page_type)
        pad = self.new_page(page_title, div)
        subpads = []
        for i, level in enumerate(levels):
            subpad = pad.cd(i+1)
            subpads.append(subpad)
            h = histos.Get('%s%s' % (histo_base_name, rec_level_code(level)))
            if h is not None:
                if prof2rms: # and type(h) == type(TProfile) should check this, but not so straightforward
                    h = make_rms_hist(h)
                if log_scale and h.GetEntries() > 0: subpad.SetLogy(1)
                apply_hist_commands(h, hist_cmds)
                h.Draw(draw_opt)
                if fit_gaus and fit_other is None:
                    h.Fit('gaus', 'Q')
                if fit_other is not None:
                    h.Fit(*fit_other)
        return subpads

    def close(self):
        if self.closed:
            return
        self.canvas.Update()
        self.ps.Close()
        if self.as_pdf:
            print 'Converting to PDF...'
            os.system('ps2pdf %s' % self.filename)
            os.system('rm %s' % self.filename)
        else:
            # New ROOT TPostScript breaks gv page number titles.
            print 'PSDrawer through with file, sedding titles...'
            os.system("sed --in-place -e 's/Page: (number /Page: (/g' %s" % self.filename)
        self.closed = True

class PSDrawerIterator:
    def __init__(self, psd, title, div=(2,2)):
        self.psd = psd
        self.title = title
        self.div = div
        self.pagesize = div[0]*div[1]
        self.cd = 0
        self.pad = None
        self.pagecount = 0

    def next(self):
        if self.pad:
            self.pad.Update()
            self.psd.canvas.Update()
        self.cd += 1
        if self.cd > self.pagesize or self.pad == None:
            self.pagecount += 1
            t = self.title
            if self.pagecount > 1:
                t += ' (%i)' % self.pagecount
            self.pad = self.psd.new_page(t, self.div)
            self.cd = 1
        return self.pad.cd(self.cd)

if len(recLevels) > 10:
    # Should reorganize using PSDrawerIterator, but for now settle
    # for a hack.
    PSDrawer.divs['all']      = (3,6)
    PSDrawer.divs['offline']  = (3,5)
    PSDrawer.divs['no_gen']   = (4,4)
    PSDrawer.divs['no_trig']  = (2,7)
    PSDrawer.divs['cocktail'] = (2,4)

__all__ = ['PSDrawer', 'PSDrawerIterator']

if __name__ == '__main__':
    psd = PSDrawer('test_psdrawer.ps')
    it = PSDrawerIterator(psd, 'testing')
    h = TH1F('test','test',100,-5,5)
    h.FillRandom('gaus', 10000)
    h2 = TH1F('test2','test2',100,-1,1)
    h2.FillRandom('gaus', 10000)
    for i in xrange(7):
        it.next()
        if i % 2:
            h.Draw()
        else:
            h2.Draw()
    psd.close()
