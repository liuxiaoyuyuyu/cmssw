import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_cfi import *
from RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_cfi import *
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import *
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
pixeltrackerlocalrecoTask = cms.Task(siPixelClusters,siPixelRecHits)
striptrackerlocalrecoTask = cms.Task(siStripZeroSuppression,siStripClusters,siStripMatchedRecHits)
trackerlocalrecoTask = cms.Task(pixeltrackerlocalrecoTask,striptrackerlocalrecoTask)

from RecoLocalTracker.SiPhase2Clusterizer.phase2TrackerClusterizer_cfi import *
from RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEGeometricESProducer_cfi import *
from RecoLocalTracker.SiPhase2VectorHitBuilder.siPhase2RecHitMatcher_cfi import *
from RecoLocalTracker.SiPhase2VectorHitBuilder.siPhase2VectorHits_cfi import *
from RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi import *

_pixeltrackerlocalrecoTask_phase2 = pixeltrackerlocalrecoTask.copy()
_pixeltrackerlocalrecoTask_phase2.add(siPhase2Clusters)
_pixeltrackerlocalrecoTask_phase2.add(siPhase2RecHits)
#_pixeltrackerlocalrecoTask_phase2.add(siPhase2VectorHits)
phase2_tracker.toReplaceWith(pixeltrackerlocalrecoTask, _pixeltrackerlocalrecoTask_phase2)
phase2_tracker.toReplaceWith(trackerlocalrecoTask, trackerlocalrecoTask.copyAndExclude([striptrackerlocalrecoTask]))

trackerlocalreco = cms.Sequence(trackerlocalrecoTask)

