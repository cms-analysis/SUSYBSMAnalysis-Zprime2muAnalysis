import os, sys

sys.argv.append('-b') # start ROOT in batch mode
from ROOT import *
sys.argv.remove('-b') # and don't mess up sys.argv

class PSDrawer:
    def __init__(self, filename, datePages=False):
        gROOT.SetStyle("Plain");
        gStyle.SetFillColor(0);
        if datePages:
            gStyle.SetOptDate();
        gStyle.SetOptStat(111111);
        gStyle.SetOptFit(1111);
        gStyle.SetPadTickX(1);
        gStyle.SetPadTickY(1);
        gStyle.SetMarkerSize(.1);
        gStyle.SetMarkerStyle(8);
        gStyle.SetGridStyle(3);
        gStyle.SetPaperSize(TStyle.kA4);
        gStyle.SetStatW(0.25);        # width of statistics box; default is 0.19
        gStyle.SetStatFormat("6.4g"); # leave default format for now
        gStyle.SetTitleFont(52,"XY"); # italic font for axis
        gStyle.SetLabelFont(52,"XY"); # italic font for axis labels
        gStyle.SetStatFont(52);       # italic font for stat. box
  
        self.canvas = TCanvas('c1','',0,0,500,640)
        self.ps = TPostScript(filename, 111)
        self.filename = filename
        self.page = 1    
        self.t = TText()
        self.t.SetTextFont(32)
        self.t.SetTextSize(0.025)

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

    def rec_level_page(self, histos, page_type, histo_base_name, page_title, draw_opt='', log_scale=False, fit_gaus=False):
        if page_type == 'all':
            div = (2,5)
            levels = xrange(9)
        elif page_type == 'offline':
            div = (2,3)
            levels = xrange(4,9)
        elif page_type == 'no_gen':
            div = (2,4)
            levels = xrange(1,9)
        elif page_type == 'no_trig':
            div = (2,3)
            levels = [0] + range(4,9)
        else:
            raise ValueError, 'page_type %s not recognized' % page_type

        pad = self.new_page(page_title, div)
        for i, level in enumerate(levels):
            subpad = pad.cd(i+1)
            if log_scale: subpad.SetLogy(1)
            h = histos.Get('%s%i' % (histo_base_name, level))
            if h is not None:
                h.Draw(draw_opt)
                if fit_gaus:
                    h.Fit('gaus', 'Q')

    def close(self):
        self.canvas.Update()
        self.ps.Close()
        # New ROOT TPostScript breaks gv page number titles.
        os.system("sed --in-place -e 's/Page: (number /Page: (/g' %s" % self.filename)
        
__all__ = ['PSDrawer']
