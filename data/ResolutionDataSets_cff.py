import FWCore.ParameterSet.Config as cms

dataSets = cms.PSet(

DY40 = cms.PSet(
  peakMass = cms.double(91.188),
  lowerMassWin = cms.double(40.0),
  upperMassWin = cms.double(500.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(0.3) # in TeV
),

DY200 = cms.PSet(
  peakMass = cms.double(200.0),
  lowerMassWin = cms.double(0.0),
  upperMassWin = cms.double(1000.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(5.4),
),

DY500 = cms.PSet(
  peakMass = cms.double(500.0),
  lowerMassWin = cms.double(0.0),
  upperMassWin = cms.double(1000.0),
  binSize = cms.int32(50),
  maxTrigMass = cms.double(5.4),
),

Zp1000 = cms.PSet(
  peakMass = cms.double(1000.0),
  lowerMassWin = cms.double(400.0),
  upperMassWin = cms.double(1600.0),
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
