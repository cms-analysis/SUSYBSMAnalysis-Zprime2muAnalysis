#!/usr/bin/env python
import sys, os, os.path, re
from math import *

import ROOT
ROOT.gROOT.SetBatch(True)

##    ____        __ _                          _   _                 
##   |  _ \  ___ / _(_)_ __   ___    ___  _ __ | |_(_) ___  _ __  ___ 
##   | | | |/ _ \ |_| | '_ \ / _ \  / _ \| '_ \| __| |/ _ \| '_ \/ __|
##   | |_| |  __/  _| | | | |  __/ | (_) | |_) | |_| | (_) | | | \__ \
##   |____/ \___|_| |_|_| |_|\___|  \___/| .__/ \__|_|\___/|_| |_|___/
##                                       |_|                          
##   
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] file directory")
parser.add_option("-o", "--out",       dest="out", type="string", default="plots", 
                        metavar="DIR", help="write output plots to DIR")
parser.add_option("-r", "--reference-file", dest="ref", type="string", 
                        metavar="FILE", help="take reference plots from FILE")
parser.add_option("--rd", "--reference-dir", dest="refdir", type="string", 
                        metavar="DIR", help="take reference plots from directory DIR")
parser.add_option("-n", "--normalize", dest="norm", type="string", default="integral",
                        help="How to normalize the plots: "
                             "'none' = no normalization, "+
                             "'integral' = normalize each plot to the same integral (default), "+
                             "'external' = normalize to the input specified in the step 1 (usually number of events), +"
                             "'manual,<value>' = rescale MC by <value> (computed in whatever way you want)")
parser.add_option("-t", "--titles", dest="titles", type="string", default="inclusiveMuonPlots_titles.txt",
                        metavar="TXTFILE", help="read axis and plot titles from TXTFILE")
parser.add_option("-s", "--select", dest="select", type="string", action="append",
                        metavar="REGEXP", help="select only plots whose name matches REGEXP; can be used more than once (will print plots that match ANY of them)")
parser.add_option("-x", "--exclude", dest="exclude", type="string",  action="append",
                        metavar="REGEXP", help="exclude plots whose name matches REGEXP; applied AFTER the 'select'; can be used more than once (will veto plots that match ANY of them)")
parser.add_option("--sg", "--select-group", dest="selectGroup", type="string", action="append",
                        metavar="GROUP[,GROUP2,..]", help="select only plots in group GROUP (you can specify this option multiple times")
parser.add_option("--xg", "--exclude-group", dest="excludeGroup", type="string",  action="append",
                        metavar="GROUP[,GROUP2,..]", help="exclude all plots in group GROUP (you can specify this option multiple times)")
parser.add_option("-S", "--stat", dest="showStat",   action="store_true", help="add statistics box to the plots")
parser.add_option("-R", "--ratio", dest="plotRatio", action="store_true", help="add a plot of the ratio this/reference")
parser.add_option("-C", "--cut", dest="plotCut", action="store_true", help="add a plot of the cut efficiency vs cut value")
parser.add_option("-O", "--overflow", dest="showOverflow", action="store_true", help="add overflows and underflows to the two outermost bins")
parser.add_option("-c", "--composite", dest="composite", type="string", help="compose reference histogram by stacking up different subhistograms")

##    ____                        ___                _ _     _       _                     _   _                 
##   |  _ \ __ _ _ __ ___  ___   ( _ )   __   ____ _| (_) __| | __ _| |_ ___    ___  _ __ | |_(_) ___  _ __  ___ 
##   | |_) / _` | '__/ __|/ _ \  / _ \/\ \ \ / / _` | | |/ _` |/ _` | __/ _ \  / _ \| '_ \| __| |/ _ \| '_ \/ __|
##   |  __/ (_| | |  \__ \  __/ | (_>  <  \ V / (_| | | | (_| | (_| | ||  __/ | (_) | |_) | |_| | (_) | | | \__ \
##   |_|   \__,_|_|  |___/\___|  \___/\/   \_/ \__,_|_|_|\__,_|\__,_|\__\___|  \___/| .__/ \__|_|\___/|_| |_|___/
##                                                                                  |_|                          
##   
(options, args) = parser.parse_args()
if len(args)<2:
    parser.print_usage() 
    print "  Use option -h or --help to get a list of the available options";
    sys.exit(2)
if options.norm.startswith("manual,"):
    try:
        options.norm_value = float(options.norm[len("manual,"):])
    except ValueError:
        print "When using --normalize=manual,<value>, value must be a valid floating point number ('%s' is not)" % options.norm
        sys.exit(2)
## For options that take multiple values, split using comma and join again
if options.selectGroup:  options.selectGroup  = sum([i.split(",") for i in options.selectGroup],  [])
if options.excludeGroup: options.excludeGroup = sum([i.split(",") for i in options.excludeGroup], [])
print ""
if options.selectGroup:
    for M in options.selectGroup: print "SelectGroup %s " %  M


##    ___       _ _   _       _ _           ____   ___   ___ _____                   _    __ _ _           
##   |_ _|_ __ (_) |_(_) __ _| (_)_______  |  _ \ / _ \ / _ \_   _|   __ _ _ __   __| |  / _(_) | ___  ___ 
##    | || '_ \| | __| |/ _` | | |_  / _ \ | |_) | | | | | | || |    / _` | '_ \ / _` | | |_| | |/ _ \/ __|
##    | || | | | | |_| | (_| | | |/ /  __/ |  _ <| |_| | |_| || |   | (_| | | | | (_| | |  _| | |  __/\__ \
##   |___|_| |_|_|\__|_|\__,_|_|_/___\___| |_| \_\\___/ \___/ |_|    \__,_|_| |_|\__,_| |_| |_|_|\___||___/
##                                                                                                         
## === GLOBAL VARIABLES ===
fileIn = ROOT.TFile(args[0])
dirIn  = fileIn.Get(args[1])
print "The directory: %s" % (args[1])
## Reference
fileRef = None; dirRef = None;
composite = []
## Information to be added to HTML page (e.g. numbers)
info   = [] 
## Titles, labels, plot groups
titles = {}
groups = {"other":[]}; groupTitles = {"other":"Other variables"}; groupToPlot = {}
index = []
## Things for plotting: canvases, lines, ..
c1 = None  ## Canvas
line = ROOT.TLine(0.,0.,10.,10.);
line.SetLineColor(2);
line.SetLineWidth(3);


## === Open references ===
if options.ref != None:
    fileRef = ROOT.TFile(options.ref)
    if options.refdir == None: options.refdir = args[1]
    dirRef  = fileRef.Get(options.refdir)
    if dirRef == None: raise RuntimeError, "Reference directory %s not found in reference file %s" % (options.refdir, options.ref)

## === Open individual components of references ===
if options.composite:
    compPattern = re.compile("([A-Z]\w*)(\(\d+\))?")
    for Name,Col in compPattern.findall(options.composite):
        color = 1;
        if Col != '': color = int(Col[1:-1])
        roocol = ROOT.gROOT.GetColor(color)
        if roocol == None:
            print "Color %d is not a valid ROOT color " % color
        compDir = fileRef.Get(options.refdir+Name)
        if compDir == None:
            print "Directory %s not found in reference file." % (options.refdir+Name)
            continue
        composite.append( (Name, compDir, color, "rgb(%d,%d,%d)" % (255*roocol.GetRed(), 255*roocol.GetGreen(), 255*roocol.GetBlue())) )

##    ____  _         _                     _       _           _    __                  _   _                 
##   / ___|| |_ _   _| | ___       _ __ ___| | __ _| |_ ___  __| |  / _|_   _ _ __   ___| |_(_) ___  _ __  ___ 
##   \___ \| __| | | | |/ _ \_____| '__/ _ \ |/ _` | __/ _ \/ _` | | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
##    ___) | |_| |_| | |  __/_____| | |  __/ | (_| | ||  __/ (_| | |  _| |_| | | | | (__| |_| | (_) | | | \__ \
##   |____/ \__|\__, |_|\___|     |_|  \___|_|\__,_|\__\___|\__,_| |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
##              |___/                                                                                          
##   
def dataStyle(histo):    
    histo.SetLineWidth(2)
    histo.SetLineColor(1)
    histo.SetMarkerStyle(20)
    histo.SetMarkerSize(1.0)
    histo.SetMarkerColor(1)
def tdrStyle():
  ROOT.gStyle.SetCanvasBorderMode(0);
  ROOT.gStyle.SetCanvasColor(0);
  ROOT.gStyle.SetCanvasDefH(600); #Height of canvas
  ROOT.gStyle.SetCanvasDefW(600); #Width of canvas
  ROOT.gStyle.SetCanvasDefX(0);   #POsition on screen
  ROOT.gStyle.SetCanvasDefY(0);
  # For the Pad:
  ROOT.gStyle.SetPadBorderMode(0);
  # ROOT.gStyle.SetPadBorderSize(Width_t size = 1);
  ROOT.gStyle.SetPadColor(0);
  ROOT.gStyle.SetPadGridX(False);
  ROOT.gStyle.SetPadGridY(False);
  ROOT.gStyle.SetGridColor(0);
  ROOT.gStyle.SetGridStyle(3);
  ROOT.gStyle.SetGridWidth(1);
  # For the frame:
  ROOT.gStyle.SetFrameBorderMode(0);
  ROOT.gStyle.SetFrameBorderSize(1);
  ROOT.gStyle.SetFrameFillColor(0);
  ROOT.gStyle.SetFrameFillStyle(0);
  ROOT.gStyle.SetFrameLineColor(1);
  ROOT.gStyle.SetFrameLineStyle(1);
  ROOT.gStyle.SetFrameLineWidth(1);
  # For the histo:
  ROOT.gStyle.SetHistLineColor(1);
  ROOT.gStyle.SetHistLineStyle(0);
  ROOT.gStyle.SetHistLineWidth(1);
  ROOT.gStyle.SetEndErrorSize(2);
  ROOT.gStyle.SetMarkerStyle(20);
  #For the date:
  ROOT.gStyle.SetOptDate(0);
  # For the statistics box:
  ROOT.gStyle.SetOptFile(0);
  #ROOT.gStyle.SetOptStat(0);
  ROOT.gStyle.SetOptStat("eou");
  ROOT.gStyle.SetStatColor(0);
  ROOT.gStyle.SetStatFont(42);
  ROOT.gStyle.SetStatFontSize(0.04);#/---> ROOT.gStyle.SetStatFontSize(0.025);
  ROOT.gStyle.SetStatTextColor(1);
  ROOT.gStyle.SetStatFormat("6.4g");
  ROOT.gStyle.SetStatBorderSize(1);
  ROOT.gStyle.SetStatH(0.15);
  ROOT.gStyle.SetStatW(0.3);#/---> ROOT.gStyle.SetStatW(0.15);
  # ROOT.gStyle.SetStatStyle(Style_t style = 1001);
  # ROOT.gStyle.SetStatX(Float_t x = 0);
  # ROOT.gStyle.SetStatY(Float_t y = 0);
  # Margins:
  ROOT.gStyle.SetPadTopMargin(0.05);
  ROOT.gStyle.SetPadBottomMargin(0.13);
  ROOT.gStyle.SetPadLeftMargin(0.16);
  ROOT.gStyle.SetPadRightMargin(0.04);
  # For the Global title:
  ROOT.gStyle.SetOptTitle(0);
  # For the axis titles:
  ROOT.gStyle.SetTitleColor(1, "XYZ");
  ROOT.gStyle.SetTitleFont(42, "XYZ");
  ROOT.gStyle.SetTitleSize(0.06, "XYZ");
  ROOT.gStyle.SetTitleXOffset(0.9);
  ROOT.gStyle.SetTitleYOffset(1.25);
  # For the axis labels:
  ROOT.gStyle.SetLabelColor(1, "XYZ");
  ROOT.gStyle.SetLabelFont(42, "XYZ");
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
  ROOT.gStyle.SetLabelSize(0.05, "XYZ");
  # For the axis:
  ROOT.gStyle.SetAxisColor(1, "XYZ");
  ROOT.gStyle.SetStripDecimals(True);
  ROOT.gStyle.SetTickLength(0.03, "XYZ");
  ROOT.gStyle.SetNdivisions(510, "XYZ");
  ROOT.gStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
  ROOT.gStyle.SetPadTickY(1);
  # Postscript options:
  ROOT.gStyle.SetPaperSize(20.,20.);
  ## OVERRIDES
  ROOT.gStyle.SetHistMinimumZero(1)
  ROOT.gStyle.SetErrorX(0.5); 
  if options.showStat == None:
    ROOT.gStyle.SetOptStat(False)
  # Done
  ROOT.gROOT.ForceStyle();

## === Convert ROOT's text specials in HTML ===
def root2html(title):
    htm = "%s" % title
    for s,r in [ ("<", "&lt;"), (">", "&gt;") ]: 
        htm = htm.replace(s,r)
    htm = re.sub(r"_{(.*?)}",  r"<sub>\1</sub>", htm) # fix special characters
    htm = re.sub(r"\^{(.*?)}", r"<sup>\1</sup>", htm) # fix special characters
    htm = re.sub(r"#(\w+)", r"&\1;", htm) # fix special characters
    htm = htm.replace("&DeltaR;", "&Delta;R") # fix deltaR
    return htm
def htmlInline(html):
    return re.sub(r"</?\w+>", "", html)

##    _   _ _     _                                                          _             _       _   _             
##   | | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___    _ __ ___   __ _ _ __ (_)_ __  _   _| | __ _| |_(_) ___  _ __  
##   | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \  | '_ ` _ \ / _` | '_ \| | '_ \| | | | |/ _` | __| |/ _ \| '_ \ 
##   |  _  | \__ \ || (_) | (_| | | | (_| | | | | | | | | | | | | (_| | | | | | |_) | |_| | | (_| | |_| | (_) | | | |
##   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_| |_| |_| |_|\__,_|_| |_|_| .__/ \__,_|_|\__,_|\__|_|\___/|_| |_|
##                         |___/                                              |_|                                    
##   
## === Plot, possibly with references ===
def plot(name, histo):
    dataStyle(histo)
    refs = getrefs(histo,name)
    maybeOverflow(histo,refs)
    c1.SetLogy(0)
    minmax(histo, refs, False) # if necessary, set ranges so that both plots are in the schema
    histo.Draw("E")
    if refs != None: stack(histo, refs)
    printHisto(name,histo.GetTitle(),"linear")
    c1.SetLogy(1)
    minmax(histo, refs, True) # if necessary, set ranges so that both plots are in the schema
    histo.Draw("E")
    if refs != None: stack(histo, refs)
    printHisto(name,histo.GetTitle(),"log")
    if options.plotCut:
        c1.SetLogy(0)
        hgtcutS = gtCut(histo,"signal")
        hgtcutS.Draw("E")
        if refs != None:
            hgtcutB = gtCut(refs[0],"bkgd")
            hgtcutB.SetLineColor(2)
            hgtcutB.Draw("E SAME")
        printHisto(name,hgtcutS.GetTitle(),"gtCut_lin")
        c1.SetLogy(0)
        hltcutS = ltCut(histo,"signal")
        hltcutS.Draw("E")
        if refs != None:
            hltcutB = ltCut(refs[0],"bkgd")
            hltcutB.SetLineColor(2)
            hltcutB.Draw("E SAME")
        printHisto(name,hgtcutS.GetTitle(),"ltCut_lin")
    c1.SetLogy(0)
    if refs != None and options.plotRatio:
        ratio(histo,refs)
        printHisto(name,histo.GetTitle(),"ratio")
## === Handle Overlflow ===
def addOverflowHist(histo):
    n = histo.GetNbinsX();
    if (n > 2):
        under = histo.GetBinContent(0)
        over  = histo.GetBinContent(n+1)
        histo.SetBinContent(1, histo.GetBinContent(1) + under);
        histo.SetBinContent(n, histo.GetBinContent(n) + over );
def maybeOverflow(histo, refs):
    if options.showOverflow:
        addOverflowHist(histo)
        if refs != None:
            for h in refs[0:2]: addOverflowHist(h)
            if (len(refs)>=5):
                for hi in refs[4]: addOverflowHist(hi)

## === Get maximum and minimum === 
def minmax(h,refs,logscale):
    if refs != None:
        href = refs[1]
        max = 0; 
        for b in range(1, h.GetNbinsX()+1):
            val = h.GetBinContent(b)    + h.GetBinError(b)
            ref = href.GetBinContent(b) + href.GetBinError(b)
            if max < val: max = val
            if max < ref: max = ref
        if logscale:
            h.GetYaxis().SetRangeUser(0.8, max*2)
        else:
            h.GetYaxis().SetRangeUser(0, max*1.2)
## === Stack plots and references === 
def stack(histo, refs):
    if len(refs) == 3:
        (ref, refup, refdn) = refs
        ref.Draw("E2 SAME")
        refup.Draw("H SAME")
        refdn.Draw("H SAME")
        histo.Draw("E SAME")
    elif len(refs) == 5:
        (ref, refup, refdn, stack, comps)  = refs
        stack.Draw("HF SAME")
        ref.Draw("E2 SAME")
        refup.Draw("H SAME")
        refdn.Draw("H SAME")
        histo.Draw("E SAME")

def printHisto(name, title, subname):
    for e in ["png"]: #options.exts:
        c1.Print("%s/%s_%s.%s" % (options.out, name, subname, e)) 
    index.append([name,subname,root2html(title)])


##    _   _                            _ _          _   _             
##   | \ | | ___  _ __ _ __ ___   __ _| (_)______ _| |_(_) ___  _ __  
##   |  \| |/ _ \| '__| '_ ` _ \ / _` | | |_  / _` | __| |/ _ \| '_ \ 
##   | |\  | (_) | |  | | | | | | (_| | | |/ / (_| | |_| | (_) | | | |
##   |_| \_|\___/|_|  |_| |_| |_|\__,_|_|_/___\__,_|\__|_|\___/|_| |_|
##                                                                    
##   
def normalize(hist,hdata):
    hist.Sumw2();
    if options.norm == "integral":
        if (hist.Integral() != 0):  
            scale = hdata.Integral()/hist.Integral()
            hist.Scale(scale)
            return scale
    elif options.norm == "external":
        #if externalNorm == None:
        normData = dirIn.Get("normalization")
        normRef  = dirRef.Get("normalization")
        if normData == None: print "Missing normalization for main sample."
        if normRef  == None: print "Missing normalization for reference sample."
        scale = normData.GetBinContent(1)/normRef.GetBinContent(1)
        hist.Scale(scale)
        return scale
    elif options.norm.startswith("manual,"):
        hist.Scale(options.norm_value)
        return options.norm_value
    return 1.0

def ratio(histo, refs):
    ref = refs[0]
    histo.GetYaxis().SetTitle("data/mc ratio")
    histo.Divide(ref);
    histo.GetYaxis().UnZoom();
    min = 0; max = 2.0; 
    for b in range(1, histo.GetNbinsX()+1):
        valup = histo.GetBinContent(b) + histo.GetBinError(b)
        #valdn = histo.GetBinContent(b) + histo.GetBinError(b)
        if valup > max and valup < 4: max = valup
    histo.GetYaxis().SetRangeUser(0, max*1.2)
    histo.Draw("E");
    line.DrawLine(histo.GetXaxis().GetXmin(),1,histo.GetXaxis().GetXmax(),1);
    histo.Draw("E SAME");

def gtCut(histo,name):
    hgt = histo.Clone(histo.GetName()+"_gtCut_"+name)
    hgt.Sumw2()
    hgt.GetYaxis().SetTitle("Efficiency (> cut)")
    for b in range(0,histo.GetNbinsX()+1):
        hgt.SetBinContent(b,histo.Integral(b,histo.GetNbinsX()+1))
    integ = histo.Integral(0,histo.GetNbinsX()+1)
    if integ > 1e-6: hgt.Scale(1.0/integ)
    return (hgt)

def ltCut(histo,name):
    hlt = histo.Clone(histo.GetName()+"_ltCut_"+name)
    hlt.Sumw2()
    hlt.GetYaxis().SetTitle("Efficiency (< cut)")
    for b in range(0,histo.GetNbinsX()+1):
        hlt.SetBinContent(b,histo.Integral(0,b))
    integ = histo.Integral(0,histo.GetNbinsX()+1)
    if integ > 1e-6: hlt.Scale(1.0/integ)
    return (hlt)

def getrefs(hdata, name):
    if dirRef != None:
        hist = dirRef.Get(name)
        if hist == None: raise RuntimeError, "Reference plot %s not found in reference file %s, dir %s" % (name, options.ref, options.refdir)
        if hist != None:
            scale = 1
            if options.norm != None: scale = normalize(hist,hdata)
            hup = hist.Clone(name+"_up")
            hdn = hist.Clone(name+"_dn")
            for b in range(1, hist.GetNbinsX()+1):
                hup.SetBinContent(b, hist.GetBinContent(b) + hist.GetBinError(b))
                hdn.SetBinContent(b, hist.GetBinContent(b) - hist.GetBinError(b))
                hup.SetBinError(b, 0)
                hdn.SetBinError(b, 0)
            ## FIXME read real options
            hist.SetMarkerStyle(0)
            if options.composite:
                components = []
                stack = ROOT.THStack(name+"_stk", name+"_stk");
                for (compName,compDir,compCol,compColName) in composite:
                    hi = compDir.Get(name);
                    if hi == None: raise RuntimeError, "Reference plot %s not found in reference file %s, dir %s%s" % (name, options.ref, options.refdir, compName)
                    hi.Scale(scale)
                    for b in range(1, hist.GetNbinsX()+1): hi.SetBinError(b, 0)
                    hi.SetFillColor(compCol)
                    components.append(hi)
                    stack.Add(hi)
                hup.SetLineColor(12)
                hdn.SetLineColor(12)
                hist.SetFillColor(12)
                hist.SetFillStyle(3013)
                return (hist, hup, hdn, stack, components)
            else:
                hup.SetLineColor(2)
                hdn.SetLineColor(2)
                hist.SetFillColor(208)
                return (hist, hup, hdn)
    return None
   
def printStats(name, histo):
    global info;
    ndata = histo.GetEntries();
    info += [ "Muons: %.0f +/- %.0f" % (histo.GetEntries(), sqrt(histo.GetEntries())) ]
    refs = getrefs(histo,name)
    #if refs != None : info += [ "refs %s" % (refs[0].GetName()) ]
    if refs != None and options.norm != "integral":
        scale = 1
        if options.norm == "external":
            normData = dirIn.Get("normalization")
            normRef  = dirRef.Get("normalization")
            if normData == None: print "Missing normalization for main sample."
            if normRef  == None: print "Missing normalization for reference sample."
            scale = normData.GetBinContent(1)/normRef.GetBinContent(1)
            info += [ "Scale: data %.0f, mc %f, ratio %.4f" % (normData.GetBinContent(1), normRef.GetBinContent(1), scale) ]
        elif options.norm.startswith("manual,"):
            scale = options.norm_value
            info += [ "Scale: %.4f (by hand)" % scale ]
        nmc = refs[0].GetEntries();
        try:
            ratio  = ndata/(scale*nmc);
            dratio = ratio * sqrt(1.0/ndata + 1.0/nmc);
            #info += [  "nmc * scale = %.0f * %.0f = %.0f" % (scale, nmc, scale*nmc) ]
            info += [  "Normalization:  %s %.0f +/- %.0f,  %s %.0f +/- %.0f, ratio %.4f +/- %.4f" % ( args[1], ndata, sqrt(ndata), options.refdir, scale*nmc, scale*sqrt(nmc), ratio, dratio ) ]
        except ValueError:
            info += [  "Normalization:_ data %.0f +/- %.0f, mc %.0f +/- %.0f" % ( ndata, sqrt(ndata), scale*nmc, scale*sqrt(nmc) ) ]
    if refs != None and options.composite:
        fracts = []
        compHistos = refs[4]
        for i,hi in enumerate(compHistos):
            try:
                fract = hi.Integral()*100.0/refs[0].Integral()
                fracts.append("%s %.1f%%" % (composite[i][0], fract))
            except ValueError:
                fracts.append("%s N/A" % composite[i][0])
        fracts.reverse()
        info += [ "Composition: " + (", ".join(fracts)) ]

##    _____ _ _   _           
##   |_   _(_) |_| | ___  ___ 
##     | | | | __| |/ _ \/ __|
##     | | | | |_| |  __/\__ \
##     |_| |_|\__|_|\___||___/
##                            
##   
def readTitles():
    group = "other"
    file = open(options.titles, "r")
    for line in file:
        gm = re.search(r"\[(\w+)\s*:\s*(.*)\]", line)
        if gm:
            group = gm.group(1); groups[group] = []
            groupTitles[group] = gm.group(2)
            continue
        fields = re.split(r"\s*:\s*", line.strip())
        if len(fields) == 3:
            titles[fields[0]] = (fields[1], fields[2])
            groups[group].append(fields[0]); groupToPlot[fields[0]] = group

## === Define Axis Labels and Titles ===
def axesAndTitles(name, histo):
    if titles.has_key(name):
        (xtitle, ptitle) = titles[name]
        histo.GetXaxis().SetTitle(xtitle)
        histo.SetTitle(ptitle)
    else:
        if options.titles:
            sys.stderr.write("Missing title for plot %s in title file %s\n" % (name, options.titles))
        histo.GetXaxis().SetTitle("muon "+name)
        histo.SetTitle(name)
    if options.norm == "integral":
        histo.GetYaxis().SetTitle("muons (entry norm.)")
    elif options.norm == "external":
        histo.GetYaxis().SetTitle("muons (event norm.)")
    else:
        histo.GetYaxis().SetTitle("muons")

##    __  __    _    ___ _   _    ____ ___  ____  _____ 
##   |  \/  |  / \  |_ _| \ | |  / ___/ _ \|  _ \| ____|
##   | |\/| | / _ \  | ||  \| | | |  | | | | | | |  _|  
##   | |  | |/ ___ \ | || |\  | | |__| |_| | |_| | |___ 
##   |_|  |_/_/   \_\___|_| \_|  \____\___/|____/|_____|
##                                                      
##   
if __name__ == "__main__":
    if not os.path.isdir(options.out): os.mkdir(options.out)
    if options.titles: readTitles()
    ##    ____       _       _                    _            _       _        
    ##   |  _ \ _ __(_)_ __ | |_   _ __ ___   ___| |_ __ _  __| | __ _| |_ __ _ 
    ##   | |_) | '__| | '_ \| __| | '_ ` _ \ / _ \ __/ _` |/ _` |/ _` | __/ _` |
    ##   |  __/| |  | | | | | |_  | | | | | |  __/ || (_| | (_| | (_| | || (_| |
    ##   |_|   |_|  |_|_| |_|\__| |_| |_| |_|\___|\__\__,_|\__,_|\__,_|\__\__,_|
    ##                                                                          
    ##   
    if dirIn.Get("metadata"):
        md = dirIn.Get("metadata")
        info.append("Muons: collection %s, selection cut %s" % (md.Get("muons").GetString(),  root2html(md.Get("selection").GetString())));
        if dirRef and dirRef.Get("metadata"):
            rmd = dirRef.Get("metadata")
            if (md.Get("muons").GetString() != rmd.Get("muons").GetString()) or (md.Get("selection").GetString() != md.Get("selection").GetString()):
                info.append("Reference: collection %s, selection cut %s" % (rmd.Get("muons").GetString(),  root2html(rmd.Get("selection").GetString())));
            if options.composite:
                compInfo = "Reference subdivided by muon classification:<ul>"
                for (compName, compDir, compRooCol, compHtmlCol) in composite:
                    compMetaDir = compDir.Get("metadata")
                    compInfo += "<li><b style=\"color: %s;\">%s</b>: collection %s, selection cut %s</li>" % (compHtmlCol, compName, compMetaDir.Get("muons").GetString(), root2html(compMetaDir.Get("selection").GetString()))
                compInfo += "</ul>"
                info.append(compInfo)
    tdrStyle()
    c1 = ROOT.TCanvas("c1","c1")
    ##    ____       _       _                     _             _       _   
    ##   |  _ \ _ __(_)_ __ | |_    ___  __ _  ___| |__    _ __ | | ___ | |_ 
    ##   | |_) | '__| | '_ \| __|  / _ \/ _` |/ __| '_ \  | '_ \| |/ _ \| __|
    ##   |  __/| |  | | | | | |_  |  __/ (_| | (__| | | | | |_) | | (_) | |_ 
    ##   |_|   |_|  |_|_| |_|\__|  \___|\__,_|\___|_| |_| | .__/|_|\___/ \__|
    ##                                                    |_|                
    ##   
    first = True
    for k in dirIn.GetListOfKeys():
        if k.GetClassName() != "TH1D": continue
        if k.GetName() == "normalization": continue
        if (len(groups) != 1) and options.selectGroup:
            if not (groupToPlot[k.GetName()] in options.selectGroup): continue
        if (len(groups) != 1) and options.excludeGroup:
            if groupToPlot[k.GetName()] in options.excludeGroup: continue
        if options.select:
            if len([x for x in options.select  if re.search(x, k.GetName())]) == 0: continue
        if options.exclude:
            if len([x for x in options.exclude if re.search(x, k.GetName())]) == 0: continue
        obj = dirIn.Get(k.GetName())
        ## Set up axes and titles
        axesAndTitles(k.GetName(), obj)
        ## Plot the histogram and possibly the background
        plot(k.GetName(), obj)
        if first: 
            printStats(k.GetName(), obj)
            first = False
    ##   __        __    _ _         _   _ _____ __  __ _     
    ##   \ \      / / __(_) |_ ___  | | | |_   _|  \/  | |    
    ##    \ \ /\ / / '__| | __/ _ \ | |_| | | | | |\/| | |    
    ##     \ V  V /| |  | | ||  __/ |  _  | | | | |  | | |___ 
    ##      \_/\_/ |_|  |_|\__\___| |_| |_| |_| |_|  |_|_____|
    ##                                                        
    ##   
    htm = open(options.out+"/index.html", "w")
    htm.write("""
<html>
<head>
  <title>Inclusive Muon Plots</title>
   <script src="http://ajax.googleapis.com/ajax/libs/prototype/1.6.0.3/prototype.js" type="text/javascript"></script>
   <script src="http://ajax.googleapis.com/ajax/libs/scriptaculous/1.8.2/scriptaculous.js?load=effects" type="text/javascript"></script>
<style type='text/css'>
body { font-family: "Candara", sans-serif; }
div.plots {
    display: block;
    float: left;
    border: 1px solid gray;
    margin: 3px;
}
div.pic { 
    display: block;
    float: left;
    margin: 3px;
    width: 330px;
}
div.pic img { border: none; width: 300px; }
h3 { text-align: center; 
     padding: 0 0.2em; }
a       { text-decoration: none; color: navy; }
a:hover { text-decoration: underline; color: rgb(120,0,0); }
div#menu  {
    background-color: rgb(200,200,255);
    border: 1px solid silver;
    position:fixed;
    top:0;
    right: 10em;
    padding: 2px 1em;
    font-size: smaller;
    
}
div#menu { text-align: right; }
div#menu ul { padding: 0 }
div#menu li { text-align: left;  }
div#menu li, div#menu ul { padding: 0; margin: 0; }
div#menu a  { color: rgb(0,80,0); }
div#menu > a  { font-size: larger; font-weight: bold;  }
div#menu li { list-style: none; }
ul { padding: 0 2em; margin-top: -0.8em;}
</style>
</head>
<body>
""")
    htm.write("<div class=\"info\">\n");
    for i in info: 
        htm.write("<p>%s</p>\n" % i);
        print i
    htm.write("</div>");
    oldkey = None 
    for name, subname, title in index:
        if title != oldkey:
            if oldkey != None: htm.write("</div>\n");
            htm.write('<div class="plots"><a name="%s"><h3>%s</h3></a>\n' % (name, title) )
            oldkey = title
        fname = "%s_%s.png" % (name, subname)
        htm.write('\t<div class="pic">')
        htm.write('<a href="%s" title="%s">' % (fname,htmlInline(title)))
        htm.write('<img src="%s" alt="%s"/>' % (fname,htmlInline(title)))
        htm.write('</a></div>\n')
    htm.write("</div>")
    htm.write("""
<div id='menu'>
<a href="#" onclick="Effect.toggle('list', 'blind', {duration:0.5}); return false; ">Select Plot</a>
<div id="list" style="display: none;"><div><ul>
""")
    oldkey = None
    for name, subname, title in index:
        if title != oldkey:
            htm.write("<li><a class=\"menu\" href=\"#%s\" onclick=\"Effect.BlindUp('list', {duration:0.25});\">%s</a></li>" % (name,title));
            oldkey = title
    htm.write("""
</ul></div></div>
""")
    htm.write("""
</body>
</html>
""")
