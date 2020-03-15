import pytest
import numpy as np
import sys
import logging
import importlib.util
import os.path
import pysam

homedir = os.path.expanduser('~')
lowerLimit = 10
upperLimit = 90

# setup path to bam and fasta reference files
PROJ_DIR = 'Documents/lab/DrugResistance/splitStrains'
PATH_TO_BAM = f'{homedir}/{PROJ_DIR}/data/mixed_synth_samples/aligned/trimmed-70.sorted.bam'
PATH_TO_REF = f'{homedir}/{PROJ_DIR}/refs/tuberculosis-tmp.fna'
PATH_TO_GFF = f'{homedir}/{PROJ_DIR}/refs/tuberculosis-tmp.gff'
PATH_TO_SPLITSTRAINS = f'{homedir}/{PROJ_DIR}/splitStrains.py'

# import splitStrains.py
spec = importlib.util.spec_from_file_location("splitStrains", PATH_TO_SPLITSTRAINS)
splitStrains = importlib.util.module_from_spec(spec)
spec.loader.exec_module(splitStrains)


@pytest.fixture
def loadData():
    fileDir = 'freqVec.csv'
    # many entries in freqVec1
    freqVec1 = np.loadtxt(open(f'{homedir}/{PROJ_DIR}/tests/test-data/freqVec1.csv', 'r'), delimiter=',')
    # 2 entries in freqVec2
    freqVec2 = np.loadtxt(open(f'{homedir}/{PROJ_DIR}/tests/test-data/freqVec2.csv', 'r'), delimiter=',')
    return freqVec1, freqVec2

@pytest.fixture
def buildGMM(loadData):
    freqVec1, freqVec2 = loadData
    components = 2
    gmm = splitStrains.fitData(freqVec1, components)
    return gmm

@pytest.fixture
def loadSam():
    samfile = pysam.AlignmentFile(PATH_TO_BAM)     # read bam file
    refFile = pysam.FastaFile(PATH_TO_REF)      # read fasta reference file
    refName = samfile.references[0]
    refLength = samfile.lengths[0]
    return (samfile, refFile)


# @pytest.mark.skip
def test_convolveVec(loadData):
    freqVec1, freqVec2 = loadData
    freqVecFlat1 = freqVec1[:,:-2].flatten()
    freqVecFlat2 = freqVec2[:,:-2].flatten()
    emptyVec = np.array([10,10,10,10, 10, 20, 30, 30, 30, 30, 30, 30, 25, 25, 25, 25])
    proprtionCountThresh = 2
    boxPoints = 4
    assert splitStrains.convolveVec(emptyVec, proprtionCountThresh, boxPoints).size == 0
    assert splitStrains.convolveVec(freqVecFlat1, proprtionCountThresh, boxPoints).size > 0
    assert splitStrains.convolveVec(freqVecFlat2, proprtionCountThresh, boxPoints).size == 0


@pytest.mark.parametrize("entropy_thresh, entropy_window, sum_of_proportions", [
    (1, 1, 0), (2, 80, 0), (0.1, 80, 1), (1,100, 0), (10,1000, 0)])

# @pytest.mark.skip
def test_entropyFilter(loadData, entropy_thresh, entropy_window, sum_of_proportions):
    freqVec1, freqVec2 = loadData
    freqVec2, entropyVec2 = splitStrains.entropyFilter(freqVec2, entropy_thresh, lowerLimit, upperLimit, entropy_window)
    freqVec2 = np.array(freqVec2)
    print(freqVec2)
    if sum_of_proportions == 0:
        for proportionVec in freqVec2[:,:-2]:
            assert sum(proportionVec) > 0
    else:
        for proportionVec in freqVec2[:,:-2]:
            assert sum(proportionVec) < 0


@pytest.mark.parametrize("depthThreshold, entropy_thresh, entropy_window, lower_bound_size", [
    (1, 1, 80, 1), (10, 1, 80, 1), (100, 0.1, 80, 1), (200, 1,100, 0), (300, 1,1000, 0)])

# @pytest.mark.skip
def test_filtervec(loadData, depthThreshold, entropy_thresh, entropy_window, lower_bound_size):
    freqVec1, freqVec2 = loadData
    freqVec1Filtered, beforeEntropyFrecVec1, entropyVec1 = splitStrains.filterVec(freqVec1, depthThreshold, entropy_thresh, entropy_window, lowerLimit, upperLimit)
    freqVec2Filtered, beforeEntropyFrecVec2, entropyVec2 = splitStrains.filterVec(freqVec2, depthThreshold, entropy_thresh, entropy_window, lowerLimit, upperLimit)

    if lower_bound_size == 1:
        assert len(freqVec1Filtered) >= lower_bound_size
    else:
        assert len(freqVec2Filtered) == lower_bound_size


@pytest.mark.parametrize('outputDir, figureName, result', [('/', 'tmp', 1), (f'{homedir}/Desktop', 'tmp', 0), (f'{homedir}/tmp', 'tmp', 1)])

# @pytest.mark.skip
def test_plotHist(loadData, buildGMM, outputDir, figureName, result):
    freqVec1, freqVec2 = loadData
    freqVecFlat1 = freqVec1[:,:-2].flatten()
    freqVecFlat2 = freqVec2[:,:-2].flatten()
    assert splitStrains.plotHist(outputDir, freqVecFlat1,  freqVecFlat1, buildGMM, figureName) == result


@pytest.mark.parametrize('baseQuality, mapQuality, startRegion, endRegion, resultSize', [
    (20, 20, 0, 10000, 10),
    (20, 20, 0, 10, 0),
    (0, 0, 0, 1, 0),
    (100, 100, 0, 10000, 0)
    ])

def test_computeDataFromSam(loadSam, baseQuality, mapQuality, startRegion, endRegion, resultSize):
    samfile, refFile = loadSam
    freqVec = []
    splitStrains.computeDataFromSam(freqVec, samfile, refFile, baseQuality, mapQuality, startRegion, endRegion)

    if resultSize == 10:
        assert len(freqVec) > resultSize
    else:
        assert len(freqVec) == resultSize

@pytest.mark.parametrize('gffDir, start, end, resultVec, result', [
    (PATH_TO_GFF, 0, 20000, [1-1, 2052-1, 3280-1], 0),
    ('wrong_path', 0, 20000, [], 1),
    (PATH_TO_GFF, 20000,20000, [], 0),
    (PATH_TO_GFF, 1,2, [], 0)])

def test_runFuncOverGff(gffDir, start, end, resultVec, result):
    vec = []

    # define simple test function
    def func(vec, regionStart=None, regionEnd=None):
        vec.append(regionStart)

    argv = [vec]

    success = splitStrains.runFuncOverGFF(gffDir,start, end, func, argv)
    assert vec == resultVec
    assert success == result
