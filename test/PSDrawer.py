import os
from ROOT import *

class PSDrawer:
    def __init__(self, filename, datePages=False):
        gROOT.SetStyle("Plain");
        #gStyle.SetFillColor(0);
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

    def close(self):
        self.canvas.Update()
        self.ps.Close()
        # New ROOT TPostScript breaks gv page number titles.
        os.system("sed --in-place -e 's/Page: (number /Page: (/g' %s" % self.filename)
        
__all__ = ['PSDrawer']
