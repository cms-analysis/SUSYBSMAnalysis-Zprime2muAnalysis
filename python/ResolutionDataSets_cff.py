import FWCore.ParameterSet.Config as cms

dataSets = cms.PSet(

SinglePt10 = cms.PSet(
  peakMass = cms.double(30),
  lowerMassWin = cms.double(0.0),
  upperMassWin = cms.double(100.0),
  binSize = cms.int32(2),
  maxTrigMass = cms.double(0.1) # in TeV
),

SinglePt100 = cms.PSet(
  peakMass = cms.double(200),
  lowerMassWin = cms.double(40.0),
  upperMassWin = cms.double(500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.3) # in TeV
),

SinglePt1000 = cms.PSet(
  peakMass = cms.double(2000),
  lowerMassWin = cms.double(40.0),
  upperMassWin = cms.double(500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.3) # in TeV
),

SinglePtFlat = cms.PSet(
  peakMass = cms.double(2000),
  lowerMassWin = cms.double(40.0),
  upperMassWin = cms.double(500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.3) # in TeV
),

DY40 = cms.PSet(
  peakMass = cms.double(91.188),
  lowerMassWin = cms.double(40.0),
  upperMassWin = cms.double(500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.3) # in TeV
),

DY200 = cms.PSet(
  peakMass = cms.double(200.0),
  lowerMassWin = cms.double(100.0),
  upperMassWin = cms.double(800.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.8),
),

DY500 = cms.PSet(
  peakMass = cms.double(500.0),
  lowerMassWin = cms.double(300.0),
  upperMassWin = cms.double(1200.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(1.2),
),

DY1000 = cms.PSet(
  peakMass = cms.double(1000.0),
  lowerMassWin = cms.double(700.0),
  upperMassWin = cms.double(2000.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(2.0),
),

DY2000 = cms.PSet(
  peakMass = cms.double(2000.0),
  lowerMassWin = cms.double(1200.0),
  upperMassWin = cms.double(3500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(3.5),
),

Zp1000 = cms.PSet(
  peakMass = cms.double(1000.0),
  lowerMassWin = cms.double(400.0),
  upperMassWin = cms.double(1600.0),
  binSize = cms.int32(24),
  maxTrigMass = cms.double(5.4),
),

Zp1200 = cms.PSet(
  peakMass = cms.double(1200.0),
  lowerMassWin = cms.double(600.0),
  upperMassWin = cms.double(1800.0),
  binSize = cms.int32(24),
  maxTrigMass = cms.double(5.4),
),

Zp1300 = cms.PSet(
  peakMass = cms.double(1300.0),
  lowerMassWin = cms.double(600.0),
  upperMassWin = cms.double(1800.0),
  binSize = cms.int32(24),
  maxTrigMass = cms.double(5.4),
),

Zp1500 = cms.PSet(
  peakMass = cms.double(1500.0),
  lowerMassWin = cms.double(1000.0),
  upperMassWin = cms.double(2000.0),
  binSize = cms.int32(20),
  maxTrigMass = cms.double(5.4),
),

Zp2000 = cms.PSet(
  peakMass = cms.double(2000.0),
  lowerMassWin = cms.double(1000.0),
  upperMassWin = cms.double(3000.0),
  binSize = cms.int32(20),
  maxTrigMass = cms.double(5.4),
),

Zp2500 = cms.PSet(
  peakMass = cms.double(2500.0),
  lowerMassWin = cms.double(1500.0),
  upperMassWin = cms.double(4000.0),
  binSize = cms.int32(25),
  maxTrigMass = cms.double(5.4),
),

Zp3000 = cms.PSet(
  peakMass = cms.double(3000.0),
  lowerMassWin = cms.double(1500.0),
  upperMassWin = cms.double(4500.0),
  binSize = cms.int32(30),
  maxTrigMass = cms.double(5.4),
),

Zp5000 = cms.PSet(
  peakMass = cms.double(5000.0),
  lowerMassWin = cms.double(3000.0),
  upperMassWin = cms.double(7000.0),
  binSize = cms.int32(40),
  maxTrigMass = cms.double(5.4),
),

G1500 = cms.PSet(
  peakMass = cms.double(1500.0),
  lowerMassWin = cms.double(800.0),
  upperMassWin = cms.double(2200.0),
  binSize = cms.int32(24),
  maxTrigMass = cms.double(5.4),
)

)
