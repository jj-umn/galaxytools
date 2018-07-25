#!/usr/bin/env python

from __future__ import print_function

import argparse
import re
import sys

create_opts = [
    'outputFileName',
    'useProbTable',
    'useProteinList',
    'printMRMTable',
    'remark',
    'binaryFormat',
    'writeDtaFiles',
    'writeMgfFile',
    'writePAIdent',
    'removeDecoyProteins',
    'plotSpectra',
    'minimumProbabilityToInclude',
    'maximumFDRToInclude',
    'datasetName',
    'setFragmentation',
    'setDeamidatedNXST',
    'addMzXMLFileToDatasetName',
    'centroidPeaks',
    'rawSpectraNoiseThreshold',
    'rawSpectraMaxDynamicRange',
    'minimumNumAAToInclude',
    'minimumNumPeaksToInclude',
    'skipRawAnnotation',
    'minimumDeltaCnToInclude',
    'maximumMassDiffToInclude',
    'bracketSpectra',
    'mergeBracket',
    'filterCriteria',
    'combineAction',
    'buildAction',
    'refreshDatabase',
    'reduceSpectra',
    'refreshDeleteUnmapped',
    'refreshDeleteMultimapped',
    'reannotatePeaks',
    'minimumMRMQ3MZ',
    'maximumMRMQ3MZ',
    'refreshTrypticOnly',
    'minimumNumReplicates',
    'removeDissimilarReplicates',
    'peakQuorum',
    'maximumNumPeaksUsed',
    'maximumNumReplicates',
    'maximumNumPeaksKept',
    'replicateWeight',
    'recordRawSpectra',
    'minimumNumReplicates',
    'qualityLevelRemove,',
    'qualityPenalizeSingletons',
    'qualityImmuneProbThreshold',
    'qualityImmuneMultipleEngines',
    'useBayesianDenoiser',
    'trainBayesianDenoiser',
    'denoiserMinimumSignalProb',
    'denoiserParamFile',
    'decoyConcatenate',
    'decoySizeRatio',
    'decoyPrecursorSwap',
    'normalizeRTWithLandmarks',
    'normalizeRTLinearRegression',
    'unidentifiedClusterIndividualRun',
    'unidentifiedClusterMinimumDot',
    'unidentifiedRemoveSinglyCharged',
    'unidentifiedMinimumNumPeaksToInclude',
    'unidentifiedSingletonXreaThreshold',
    'allowableModTokens'
]

filter_opts = [
    'outputFileName',
    'useProbTable',
    'useProteinList',
    'printMRMTable',
    'remark',
    'binaryFormat',
    'writeDtaFiles',
    'writeMgfFile',
    'writePAIdent',
    'removeDecoyProteins',
    'plotSpectra',
    'minimumProbabilityToInclude',
    'maximumFDRToInclude',
    'datasetName',
    'setFragmentation',
    'setDeamidatedNXST',
    'addMzXMLFileToDatasetName',
    'centroidPeaks',
    'rawSpectraNoiseThreshold',
    'rawSpectraMaxDynamicRange',
    'minimumNumAAToInclude',
    'minimumNumPeaksToInclude',
    'skipRawAnnotation',
    'minimumDeltaCnToInclude',
    'maximumMassDiffToInclude',
    'bracketSpectra',
    'mergeBracket',
    'filterCriteria',
    'combineAction',
    'buildAction',
    'refreshDatabase',
    'reduceSpectra',
    'refreshDeleteUnmapped',
    'refreshDeleteMultimapped',
    'reannotatePeaks',
    'minimumMRMQ3MZ',
    'maximumMRMQ3MZ',
    'refreshTrypticOnly',
    'minimumNumReplicates',
    'removeDissimilarReplicates',
    'peakQuorum',
    'maximumNumPeaksUsed',
    'maximumNumReplicates',
    'maximumNumPeaksKept',
    'replicateWeight',
    'recordRawSpectra',
    'minimumNumReplicates',
    'qualityLevelRemove,',
    'qualityPenalizeSingletons',
    'qualityImmuneProbThreshold',
    'qualityImmuneMultipleEngines',
    'useBayesianDenoiser',
    'trainBayesianDenoiser',
    'denoiserMinimumSignalProb',
    'denoiserParamFile',
    'decoyConcatenate',
    'decoySizeRatio',
    'decoyPrecursorSwap',
    'normalizeRTWithLandmarks',
    'normalizeRTLinearRegression',
    'unidentifiedClusterIndividualRun',
    'unidentifiedClusterMinimumDot',
    'unidentifiedRemoveSinglyCharged',
    'unidentifiedMinimumNumPeaksToInclude',
    'unidentifiedSingletonXreaThreshold',
    'allowableModTokens'
]

search_opts = [
    'libraryFile',
    'databaseFile',
    'databaseType',
    'indexCacheAll',
    # 'filterSelectedListFileName',
    'precursorMzTolerance',
    'precursorMzUseAverage',
    'searchAllCharges',
    'detectHomologs',
    'fvalFractionDelta',
    'useSp4Scoring',
    'fvalUseDotBias',
    'usePValue',
    'useTierwiseOpenModSearch',
    # 'expectedCysteineMod',
    # 'ignoreSpectraWithUnmodCysteine',
    # 'ignoreChargeOneLibSpectra',
    # 'ignoreAbnormalSpectra',
    'outputExtension',
    'outputDirectory',
    'hitListTopHitFvalThreshold',
    'hitListLowerHitsFvalThreshold',
    'hitListShowHomologs',
    'hitListShowMaxRank',
    'hitListOnlyTopHit',
    'hitListExcludeNoMatch',
    'enzymeForPepXMLOutput',
    'printFingerprintingSummary',
    'filterMinPeakCount',
    'filterAllPeaksBelowMz',
    'filterMaxIntensityBelow',
    'filterMinMzRange',
    'filterCountPeakIntensityThreshold',
    'filterRemovePeakIntensityThreshold',
    'filterMaxPeaksUsed',
    'filterMaxDynamicRange',
    'peakScalingMzPower',
    'peakScalingIntensityPower',
    'peakScalingUnassignedPeaks',
    'peakNoBinning',
    'peakBinningNumBinsPerMzUnit',
    'peakBinningFractionToNeighbor',
    'filterLibMaxPeaksUsed',
    'filterLightIonsMzThreshold',
    'filterITRAQReporterPeaks',
    'filterTMTReporterPeaks',
    # 'filterRemoveHuge515Threshold',
]


def __main__():
    parser = argparse.ArgumentParser(
        description='Parse SpectraST search.params files' +
                    ' to create an updated search.params')
    parser.add_argument(
        'param_files', nargs='*',
        help='A SpectraST search.params files')
    parser.add_argument(
        '-m', '--mode', choices=['search','create','filter'],
        help='')
    parser.add_argument(
        '-o', '--output',
        help='Output file  (-) for stdout')
    args = parser.parse_args()

    output_wtr = open(args.output, 'w')\
        if args.output and args.output != '-' else sys.stdout

    optpat = re.compile('^([a-z]\w+)\s*[=:]\s*([^=]+)$')

    valid_opts = search_opts if args.mode == 'search' else create_opts if args.mode == 'create' else filter_opts
    valid_params = dict()

    # Collect all valid_params
    def parse_params(param_file, fh, valid_opts):
        for i, line in enumerate(fh):
            try:
                m = optpat.match(line.rstrip())
                if m:
                    k, v = m.groups()
                    if k in valid_opts:
                        valid_params[k] = v
            except Exception, e:
                print('%s(%d): %s %s' % (param_file, i, line, e),
                      file=sys.stderr)

    if args.param_files:
        for param_file in args.param_files:
            try:
                with open(param_file, 'r') as fh:
                    parse_params(param_file, fh, valid_opts)
            except Exception, e:
                print('parse_params: %s' % e, file=sys.stderr)
    else:
        try:
            parse_params('stdin', sys.stdin)
        except Exception, e:
            print('parse_params: %s' % e, file=sys.stderr)

    # Write valid_params
    for valid_opt in valid_opts:
        if valid_opt in valid_params:
            print('%s = %s' % (valid_opt, valid_params[valid_opt]), file=output_wtr)


if __name__ == "__main__":
    __main__()

