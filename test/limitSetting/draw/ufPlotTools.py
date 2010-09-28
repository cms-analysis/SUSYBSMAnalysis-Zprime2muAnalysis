###############################################################################
# UF plotting tools for TDR-like style and more
###############################################################################

from ctypes import *

from array import array

# import ROOT

from ROOT import *

from setTdrStyle import *

#------------------------------------------------------------------------------
# save a canvas in triplicate - eps, gif and pdf 
###############################################################################
def tripleSave ( canv, baseName ):
   canv.SaveAs('gif/' + baseName + '.gif' )
   canv.SaveAs('eps/' + baseName + '.eps' )
   canv.SaveAs('pdf/' + baseName + '.pdf' )
###############################################################################

#------------------------------------------------------------------------------
# set the properties of a default display pad
###############################################################################
def setupPad ( pad ):
  global __pad
  __pad = pad
  __pad.SetTopMargin   (0.15)
  __pad.SetBottomMargin(0.15)
  __pad.SetLeftMargin  (0.15)
  __pad.SetRightMargin (0.15)

###############################################################################


#------------------------------------------------------------------------------
# set up the default Latex
###############################################################################
def setupLatex():
    global latex__

    latex__.SetTextColor(kBlack)
    latex__.SetTextFont(42)
    latex__.SetTextSize(0.06)
    latex__.SetTextAlign(22)
    latex__.SetNDC()

###############################################################################
    

#------------------------------------------------------------------------------
# create a new canvas with auto-displacement
###############################################################################
def newCanvas( name, xWidth=700, yWidth=700 ):
    global canvasX__
    global canvasY__

    tmpCanv = TCanvas ( name, name, canvasX__, canvasY__, xWidth, yWidth )

    canvasX__+=20
    canvasY__+=20

    if (canvasY__ > 350):
        canvasY__ = 20

    if (canvasX__ > 700):
        canvasX__ = 20
    
    tmpCanv.SetRightMargin  (0.15)
    tmpCanv.SetTopMargin    (0.15)
    tmpCanv.SetLeftMargin   (0.20)
    tmpCanv.SetBottomMargin (0.20)

    tmpCanv.GetFrame().DrawClone()
    tmpCanv.cd()

    return tmpCanv

###############################################################################

#------------------------------------------------------------------------------
# setup nice gradient for displaying greyscale 2d plots
###############################################################################
def setupGreyScale():

    NRGBs = 5
    NCont = 255

    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [0.95, 0.30, 0.20, 0.10, 0.00]
    green = [0.95, 0.30, 0.20, 0.10, 0.00]
    blue  = [0.95, 0.30, 0.20, 0.10, 0.00]

    s = array('d',stops)
    r = array('d',red  )
    g = array('d',green)
    b = array('d',blue )
  
    TColor.CreateGradientColorTable(NRGBs, s, r, g, b, NCont)

    global gStyle
  
    gStyle.SetNumberContours(NCont)

###############################################################################

#------------------------------------------------------------------------------
# setup nice gradient for displaying greyscale 2d plots, z axis log scale
###############################################################################
def setupLogGreyScale():

    NRGBs = 5
    NCont = 255

    stops = [ 0.4096, 0.512, 0.64, 0.8, 1.00 ]
    red   = [ 0.95  , 0.85, 0.60, 0.45, 0.00 ]
    green = [ 0.95  , 0.85, 0.60, 0.45, 0.00 ]
    blue  = [ 0.95  , 0.85, 0.60, 0.45, 0.00 ]

    s = array('d',stops)
    r = array('d',red  )
    g = array('d',green)
    b = array('d',blue )
  
    TColor.CreateGradientColorTable(NRGBs, s, r, g, b, NCont)

    global gStyle
  
    gStyle.SetNumberContours(NCont)

#------------------------------------------------------------------------------
# setup nice color gradient for displaying 2D plots
###############################################################################
def setupColors():

    NRGBs = 5
    NCont = 255

    stops = [ 0.00, 0.0625, 0.25, 0.5625, 1.00 ]
    red   = [ 0.00, 0.00  , 0.87, 1.00  , 0.51 ]
    green = [ 0.00, 0.81  , 1.00, 0.20  , 0.00 ]
    blue  = [ 0.51, 1.00  , 0.12, 0.00  , 0.00 ]

    s = array('d',stops)
    r = array('d',red  )
    g = array('d',green)
    b = array('d',blue )
  
    TColor.CreateGradientColorTable(NRGBs, s, r, g, b, NCont)

    global gStyle
  
    gStyle.SetNumberContours(NCont)

###############################################################################
def getAxes( h ):

    global __xAxis
    global __yAxis
    global __zAxis

    __xAxis = h.GetXAxis()
    __yAxis = h.GetYAxis()
    __zAxis = h.GetZAxis()

###############################################################################
def ufPlotTools():

    print
    print '  **************************'
    print '  ** UF PyRoot Plot Tools **'
    print '  ** ikf v1.0, June 2010  **'
    print '  **************************'
    print
    
    setupLatex()

    setupColors()

    global gStyle

    gStyle.SetNdivisions(505,"X");

    global __random

    __random = TRandom3(13)

##---------------------------------------------------------------------------- 
## this part gets executed at startup
###############################################################################

setTDRStyle()

canvasX__ = 20
canvasY__ = 20

latex__ = TLatex()

ufPlotTools()

