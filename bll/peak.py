"""Module Description

Copyright (c) 2010,2011 Zhenqing Ye <iamyezhenqing@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Zhenqing Ye
@contact: iamyezhenqing@gmail.com
"""

import os
import networkx
import numpy
from bll import util 

from bll.core import BioTasker
from bll.core import TaskDatar
from bll.worm import InchWorm
import collections
     
class PeakScanner(BioTasker):
    '''A class implemented for Peak-calling, with "threading" feature in mind.'''
    
    def __init__(self, pargs, plog):
        self.pargs   = pargs
        self.log = plog 

        # get those parameters
        self.R = self.pargs.rscan
        self.s = self.pargs.minstep
        self.S = self.pargs.maxstep
        self.pvalue = self.pargs.pvalue

        # preparing the output fold
        if (not os.path.exists(self.pargs.output + '/Peak')):
            tmp = os.makedirs(self.pargs.output + '/Peak', exist_ok=True)
        self.outDir = self.pargs.output + '/Peak'

    def configureTask(self, chromNode):
        chromDatar = chromNode.data 
        taskDatar = self.configurePairMode(chromDatar)
        '''
        if(self.pargs.mode == 'pair'):
            taskDatar = self.configurePairMode(chromDatar)
        else:
            taskDatar = self.configureSingleMode(chromDatar)
        '''
        return taskDatar


    def configurePairMode(self, chromDatar):
        taskDatar = TaskDatar()
        lociDiGraph = chromDatar.lociDiGraph
        sortedPoints = chromDatar.sortedPoints

        # for pair-mode, we need separate the positive strand
        # and negative strand from the above data-structure,
        # because initially we have "absolute" the R2 positions.
        pRefDiGraph = networkx.DiGraph()
        pHead = 'PHEAD'
        pRefDiGraph.add_node(pHead)
        nRefDiGraph = networkx.DiGraph()
        nHead = 'NHEAD'
        nRefDiGraph.add_node(nHead)
        for R2 in sortedPoints:
            toR1sDict = lociDiGraph.nodes[R2]['toR1s']
            pcountDict = collections.defaultdict(int)
            ncountDict = collections.defaultdict(int)
            for R1 in toR1sDict.keys():
                if (R1 > 0): # R1 is positive, so R2 is negative
                    ncountDict[R1] += toR1sDict[R1]
                else: # now R2 is positive
                    pcountDict[R1] += toR1sDict[R1]
            pdegree = sum(pcountDict.values())
            if (pdegree > 0):
                pRefDiGraph.add_node(R2, toR1s=pcountDict, depth=pdegree)
                pRefDiGraph.add_edge(pHead, R2)
                pHead = R2
            ndegree = sum(ncountDict.values())
            if (ndegree > 0):
                nRefDiGraph.add_node(R2, toR1s=ncountDict, depth=ndegree)
                nRefDiGraph.add_edge(nHead, R2)
                nHead = R2

        pHead = list(pRefDiGraph.successors('PHEAD'))[0]
        pRefDiGraph.remove_node('PHEAD')
        nHead = list(nRefDiGraph.successors('NHEAD'))[0]
        nRefDiGraph.remove_node('NHEAD')
        # data-structure has been transformed and attached to taskDatar object
        # it's time to clear these unnecessary data members now
        taskDatar.pRefDiGraph = pRefDiGraph
        taskDatar.pHead = pHead
        taskDatar.nRefDiGraph = nRefDiGraph
        taskDatar.nHead = nHead
        delattr(chromDatar, 'lociDiGraph')
        delattr(chromDatar, 'sortedPoints')
        return taskDatar

    def configureSingleMode(self, chromDatar):
        # wait for implementation in future
        taskDatar = TaskDatar()
        return taskDatar

    def calcEntropy(self, group, refDiGraph):
        degrees = [refDiGraph.nodes[n]['depth'] for n in group]
        degsum = sum(degrees)
        problist = [1.0*val/degsum for val in degrees ]
        entropy = sum([-prob*numpy.log2(prob) for prob in problist])
        return entropy

    def peakScanningProc(self, taskDatar, chrid):
        pRefDiGraph = taskDatar.pRefDiGraph 
        nRefDiGraph = taskDatar.nRefDiGraph
        measure = lambda e: abs(e[1] - e[0])

        pComponents = self.breakingIntoComponents(pRefDiGraph, taskDatar.pHead, measure, k=3.0)
        nComponents = self.breakingIntoComponents(nRefDiGraph, taskDatar.nHead, measure, k=3.0)

        #here we filtering those small components
        degreepf = lambda group: sum([pRefDiGraph.nodes[n]['depth'] for n in group])
        pComponents = [item for item in pComponents if degreepf(item) >= self.R]
        degreenf = lambda group: sum([nRefDiGraph.nodes[n]['depth'] for n in group])
        nComponents = [item for item in nComponents if degreenf(item) >= self.R]
        '''
        #maybe we also need filter the components by entropy further
        pEntropyList = [ self.calcEntropy(item, pRefDiGraph) for item in pComponents ]
        pMean, pStd = numpy.mean(pEntropyList), numpy.std(pEntropyList)
        pComponents = [ pComponents[i] for i in range(0, len(pComponents)) if pEntropyList[i]>pMean-2.0*pStd ]
        nEntropyList = [ self.calcEntropy(item, nRefDiGraph) for item in nComponents ]
        nMean, nStd = numpy.mean(nEntropyList), numpy.std(nEntropyList)
        nComponents = [ nComponents[i] for i in range(0, len(nComponents)) if nEntropyList[i]>nMean-2.0*nStd ]
        '''

        pPeakList, pCompSerials = [], []
        compid = 0
        for component in pComponents:
            pPeaks = self.processComponent(component, pRefDiGraph)
            pPeakList.extend(pPeaks)
            pCompSerials.extend([compid for i in range(0, len(pPeaks))])
            compid += 1
        nPeakList, nCompSerials = [], []
        compid = 0
        for component in nComponents:
            nPeaks = self.processComponent(component, nRefDiGraph)
            nPeakList.extend(nPeaks)            
            nCompSerials.extend([compid for i in range(0, len(nPeaks))])
            compid += 1
        #print("raw pPeakList: %d" % len(pPeakList))
        #print("raw nPeakList: %d" % len(nPeakList))

        #compSerials for remember which component these peaks come from
        #and this information would be used in peakMerging, only those
        #peaks from the same component would take part in the merging activity
        pmergedPeakDict = self.peakMerging(pPeakList, pCompSerials, pRefDiGraph)
        nmergedPeakDict = self.peakMerging(nPeakList, nCompSerials, nRefDiGraph)
        return pmergedPeakDict, nmergedPeakDict

    def processComponent(self, component, refDiGraph):
        #the members in the component already sorted
        #cL, cN = self.getRScanLocalContext(component, refDiGraph)
        cL, cN = self.getRScanLocalContext(component[0], component[-1], refDiGraph)
        s, S = self.pargs.minstep, self.pargs.maxstep
        R, pvalue = self.pargs.rscan, self.pargs.pvalue
        smartWorm = InchWorm(refDiGraph, component, cL, cN, R, s, S, pvalue)
        smartWorm.configure()
        smartWorm.run()
        return smartWorm.peaks

    def breakingIntoComponents(self, refDiGraph, head, measure, k=3.0):
        edgeValList = []
        current = head
        while (True):
            succ = list(refDiGraph.successors(current))
            if (succ):
                value = measure((current, succ[0]))
                edgeValList.append(value)
                current = succ[0]
            else:
                break  
        obslist, enumber = edgeValList, len(edgeValList)
        edgeFlagList = [1 for i in range(0, enumber)]
        while (True):
            valarray = numpy.array(obslist)
            mu, sigma = numpy.mean(valarray), numpy.std(valarray)
            edgeFlagList = [0 if val-mu > k*sigma else 1 for val in edgeValList]
            obslist = [val for i, val in enumerate(edgeValList) if edgeFlagList[i]]
            onumber = len(obslist)
            if (onumber < enumber):
                enumber = onumber
            else:
                break

        components = [[head]]
        current, cursor = head, 0
        while (True):
            succ = list(refDiGraph.successors(current))
            if (succ):
                if (edgeFlagList[cursor]):
                    components[-1].append(succ[0])
                else:
                    components.append([succ[0]])
                current = succ[0]
                cursor += 1
            else:
                break
        return components

    #def getRScanLocalContext(self, component, refDiGraph, extend=10000):
    def getRScanLocalContext(self, start, end, refDiGraph, extend=10000):
        #local background is based on the whole component rather than
        #each individual peak candidate. It means every peak residents
        #in this component would have the same background control.
        #start, end = component[0], component[-1]
        middle = int(0.5 * (start + end))
        span = int(0.5 * extend)
        L = end - start + 1
        #N = sum([refDiGraph.nodes[node]['depth'] for node in component])
        poslist, vallist = self.nodeDepthFlowing(start, end, refDiGraph)
        N = sum(vallist)
        while (middle-start < span):
        #if the whole component is smaller than required extension 10kb
            head = list(refDiGraph.predecessors(start))
            if (head):
                if ((middle-head[0]) <= span):
                    N += refDiGraph.nodes[head[0]]['depth']
                    L += start - head[0]
                    start = head[0]
                else:
                    break
            else:
                break
        while (end-middle < span):
            tail = list(refDiGraph.successors(end))
            if (tail):
                if ((tail[0]-middle) <= span): 
                    N += refDiGraph.nodes[tail[0]]['depth']
                    L += tail[0] - end
                    end = tail[0]
                else:
                    break
            else:
                break
        return L, N

    def outlierEdgeCutting(self, paramGraph, measure, k=3.0):
        paramEdges = list(paramGraph.edges())
        edgeValList = [measure(edge) for edge in paramEdges]
        edgeFlagList = [1 for _ in range(0, len(edgeValList))]
        obslist, enumber = edgeValList, len(edgeValList)
        while (True):
            valarray = numpy.array(obslist)
            mu, sigma = numpy.mean(valarray), numpy.std(valarray)
            edgeFlagList = [0 if val-mu > k*sigma else 1 for val in edgeValList]
            obslist = [val for i, val in enumerate(edgeValList) if edgeFlagList[i]]
            onumber = len(obslist)
            if (onumber < enumber):
                enumber = onumber
            else:
                break
        eraseEdges = [paramEdges[i] for i,flag in enumerate(edgeFlagList) if not flag]
        paramGraph.remove_edges_from(eraseEdges)
        return paramGraph

    def peakMerging(self, peakList, compSerials, refDiGraph):
        # peakList is already sorted, and each peak would have the content
        # looks like : Peak = namedtuple('Peak', 'start, end, w, r, score, pvalue')
        # to be consistent with the output from SmartWorm
        peakMgGraph = networkx.Graph()
        for i in range(0, len(peakList)-1):
            prevPeak, prevCompId = peakList[i], compSerials[i]
            nextPeak, nextCompId = peakList[i+1], compSerials[i+1]
            prevStart, prevEnd = prevPeak.start, prevPeak.end
            nextStart, nextEnd = nextPeak.start, nextPeak.end
            dist = abs(nextStart - prevEnd)
            if(prevCompId==nextCompId): #peak from the same component
                peakMgGraph.add_edge(i, i+1, dist=dist)
            else:
                peakMgGraph.add_nodes_from([i, i+1])
        measure = lambda e: peakMgGraph.get_edge_data(*e)['dist']
        #1 - norm.cdf(3.5) = 0.00023
        #here we tend to keep those peaks in the same component to be linked
        #so pick the k with a little bit higher 
        peakMgGraph = self.outlierEdgeCutting(peakMgGraph, measure, k=3.0)

        components = networkx.connected_components(peakMgGraph)
        MergedPeak = collections.namedtuple('MergedPeak', 'chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue')
        mergedPeakDict = collections.defaultdict(None)
        for comp in components:
            comp = sorted(comp)
            chromStart = peakList[comp[0]].start
            chromEnd = peakList[comp[-1]].end #included
            posList, depthList = self.nodeDepthFlowing(chromStart, chromEnd, refDiGraph)
            valList = [1.0*pos*val for pos,val in zip(posList, depthList)]
            uPoint = int(sum(valList)/sum(depthList)) # the "average" point based on density distribution
            blockCount = len(comp)
            blockSizes = [peakList[index].end - peakList[index].start + 1 for index in comp]
            blockStarts = [ peakList[index].start-chromStart for index in comp ]
            name = str(chromStart) #+ '-' + str(chromEnd) 
            #name = str(uPoint)
            cL, cN = self.getRScanLocalContext(chromStart, chromEnd, refDiGraph)
            w, r = chromEnd - chromStart + 1, sum(valList) + 1
            score, pValue = self.calculateScanScore(r, w, cL, cN)
            mgPeak = MergedPeak(chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue)
            mergedPeakDict[chromStart] = mgPeak
        return mergedPeakDict


    def nodeDepthFlowing(self, start, end, refDiGraph):
        # start always less than end
        '''
        pos = start
        while(pos <= end):
            depth = refDiGraph.nodes[pos]['depth']
            yield pos, depth
            pos = refDiGraph.successors(pos)[0]
        '''
        posList, valList = [], []
        pos = start
        while (pos <= end):
            depth = refDiGraph.nodes[pos]['depth']
            posList.append(pos)
            valList.append(depth)
            succs = list(refDiGraph.successors(pos))
            if (succs):
                pos = succs[0]
            else:
                break  # to the end
        return posList, valList

    def calculateScanScore(self, r, w, L, N):
        logx = numpy.log(w) - numpy.log(L) + (1.0 + 1.0/r) * numpy.log(N)
        score = r*logx - util.fastLogFactorial(r)
        lda = numpy.exp(score)
        prob = 1.0 - numpy.exp(-lda)
        return score, prob

    def peakPairing(self, pmPeakDict, nmPeakDict):
        #MergedPeak(chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue)
        #Peaks on the "same" strand would never overlapped
        pairDiGraph = networkx.DiGraph()
        pPoints, nPoints = sorted(pmPeakDict.keys()), sorted(nmPeakDict.keys())
        pN, nN = len(pPoints), len(nPoints)
        print(pN, nN)
        #cursor for positive and negative strands
        pc, nc = 0, 0
        #find the nearest peak on the other strand
        for pPos in pPoints:
            pUPoint = pmPeakDict[pPos].uPoint
            nPos = nPoints[nc]
            nUPoint  = nmPeakDict[nPos].uPoint
            while (nUPoint <= pUPoint and nc < nN-1):
                nc += 1
                nPos = nPoints[nc]
                nUPoint = nmPeakDict[nPos].uPoint
            GAPSIZE = abs(pUPoint - nUPoint)
            pairDiGraph.add_edge('P-' + str(pPos), 'N-' + str(nPos), dist=GAPSIZE)
        #now treat the negative strand
        for nPos in nPoints:
            nUPoint = nmPeakDict[nPos].uPoint
            pPos = pPoints[pc]
            pUPoint  = pmPeakDict[pPos].uPoint
            if (pUPoint >= nUPoint): #ignoring the first items
                continue
            while (pUPoint < nUPoint and pc < pN-1):
                pc += 1
                pPos = pPoints[pc]
                pUPoint = pmPeakDict[pPos].uPoint
            #back one step
            pc -= 1
            pPos = pPoints[pc]
            pUPoint = pmPeakDict[pPos].uPoint
            GAPSIZE = abs(nUPoint - pUPoint)
            pairDiGraph.add_edge('N-' + str(nPos), 'P-' + str(pPos), dist=GAPSIZE)
        #removing those outlier edges by inspection of their distances
        ef = lambda e: pairDiGraph.get_edge_data(*e)['dist']
        pairDiGraph = self.outlierEdgeCutting(pairDiGraph, ef, k=3.5)
        return pairDiGraph

    def goodPeakPairs(self, peakPairDiGraph):
        PNPairList = set()
        for edge in peakPairDiGraph.edges():
            source, target = edge[0], edge[1]
            if (source in list(peakPairDiGraph.successors(target)) \
              and target in list(peakPairDiGraph.successors(source))):
                pair = (source, target) if source.startswith('P-') else (target, source)
                PNPairList.add(pair)
        return list(PNPairList)

    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        #self.log.info("%s: Starting PeakCalling-task on %s..." % (time.ctime(), chrid))
        self.log.info("Starting PeakCalling-task on %s..." % (chrid))
        taskDatar = self.configureTask(chromNode)
        pmergedPeakDict, nmergedPeakDict = self.peakScanningProc(taskDatar, chrid)
        '''
        outfile = self.outDir + '/peakcall_on_' + chrid + '.bed'
        handle = open(outfile, 'w')
        for pid in pmergedPeakDict:
            mpPeak = pmergedPeakDict[pid]
            linestr = self.buildOneLineStr(chrid, '+', mpPeak, 'pp')
            handle.write(linestr)

        for pid in nmergedPeakDict:
            mnPeak = nmergedPeakDict[pid]
            linestr = self.buildOneLineStr(chrid, '-', mnPeak, 'np')
            handle.write(linestr)
        handle.close()  
        '''
        pairDiGraph = self.peakPairing(pmergedPeakDict, nmergedPeakDict)
        pairs = self.goodPeakPairs(pairDiGraph)

        unpairNodes = set(pairDiGraph.nodes())
        pNodes = set([pair[0] for pair in pairs])
        nNodes = set([pair[1] for pair in pairs])

        unpairNodes = unpairNodes.difference(pNodes)
        unpairNodes = unpairNodes.difference(nNodes)
        #print('Paired Peaks number: %d' % len(pairs))
        #print('Unpaired Peaks number: %d' % len(unpairNodes))
        #gappedPeak format (http://genome.ucsc.edu/FAQ/FAQformat.html#format14)
        #---------------------------------------
        # output section
        #---------------------------------------
        ##MergedPeak(chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue)
        outfile = self.outDir + '/peakcall_on_' + chrid + '.paired.bed'
        handle = open(outfile, 'w')
        for pair in pairs:
            mpPeakId = int(pair[0].split('-')[1])
            mpPeak = pmergedPeakDict[mpPeakId]
            linestr = self.buildOneLineStr(chrid, '+', mpPeak, 'pp')
            handle.write(linestr)
            
            mnPeakId = int(pair[1].split('-')[1])
            mnPeak = nmergedPeakDict[mnPeakId]
            linestr = self.buildOneLineStr(chrid, '-', mnPeak, 'np')
            handle.write(linestr)
        '''    
        for item in unpairNodes:
            if(item.startswith('P-')):
                upPeakId = int(item.split('-')[1])
                upPeak = pmergedPeakDict[upPeakId]
                linestr = self.buildOneLineStr(chrid, '+', upPeak, 'up')
                handle.write(linestr)
            else:
                unPeakId = int(item.split('-')[1])
                unPeak = nmergedPeakDict[unPeakId]
                linestr = self.buildOneLineStr(chrid, '-', unPeak, 'un')
                handle.write(linestr)
        '''        
        handle.close()
        #---------------------------------------
        chromDatar = chromNode.data
        chromDatar.pRefDiGraph = taskDatar.pRefDiGraph
        chromDatar.pHead = taskDatar.pHead
        chromDatar.pmergedPeakDict = pmergedPeakDict
        
        chromDatar.nRefDiGraph = taskDatar.nRefDiGraph
        chromDatar.nHead = taskDatar.nHead
        chromDatar.nmergedPeakDict = nmergedPeakDict
        chromDatar.peakPairs = pairs

        #self.log.info("%s: PeakCalling-task is finished on %s..." % (time.ctime(), chrid) ) 
        self.log.info("PeakCalling-task is finished on %s..." % (chrid))

    def buildOneLineStr(self, chrid, strand, peak, flag):
        #peak.chromEnd+1 because it's closed end here, for bed format, it's opened end.
        linestr = chrid + '\t' + str(peak.chromStart) + '\t' + str(peak.chromEnd+1) + '\t'
        if (strand == '+'):
            name = 'PP-' + peak.name
        else:
            name = 'PN-' + peak.name
        score = 0
        linestr += name + '\t' + str(score) + '\t' + strand + '\t'
        thickStart, thickEnd = str(peak.chromStart), str(peak.chromEnd+1)
        linestr += thickStart + '\t' + thickEnd + '\t'
        #flag: pp means positive strand and paired
        #      np means negative strand and paired
        #      pu means positive strand and unpaired
        #      nu means negative strand and unpaired
        if (flag == 'pp'):
            itemRgb = "255,0,0"
        elif (flag == 'np'):
            itemRgb = "0,0,255"
        else:
            itemRgb = "0,255,0"
        linestr += itemRgb + '\t'
        linestr += str(peak.blockCount) + '\t'
        blockSizes = ','.join(map(str, peak.blockSizes))
        linestr += blockSizes + '\t'
        blockStarts = ','.join(map(str, peak.blockStarts))
        linestr += blockStarts + '\n'
        return linestr
