# This is our selection as used for the 2011 EPS results. The only
# change since then (as implemented in OurSelectionNew) is switching
# from tracker hits > 10 to tracker layers > 8. We simply undo that
# change here.

from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff import loose_cut, tight_cut, allDimuons, dimuons

loose_cut = loose_cut.replace('globalTrack.hitPattern.trackerLayersWithMeasurement > 8',
                              'globalTrack.hitPattern.numberOfValidTrackerHits > 10')
                              
allDimuons = allDimuons.clone(loose_cut=loose_cut)
dimuons = dimuons.clone()
