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
import random
import networkx
import numpy
from scipy import stats
from networkx.readwrite import json_graph
from collections import defaultdict
from numpy.linalg.linalg import LinAlgError

from bll.core import BioTasker
from bll.core import TaskDatar

class BorderScanner(BioTasker):
    '''A class implemented for Border-positionate'''
    
    def __init__(self, pargs, plog):
        #========================================================
	    #  Basic data-structure setting #
	    #========================================================
        self.pargs   = pargs
        self.log = plog 
        self.c = self.pargs.chernoff
        self.k = self.pargs.outlier        
        
        # preparing the output fold
        if(not os.path.exists(self.pargs.output + '/Border')):
            tmp = os.mkdir(self.pargs.output + '/Border') 
        self.outDir = self.pargs.output + '/Border'
    
    #========================================================
	#  Real data preparation
	#========================================================
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
        
        pRefDiGraph = chromDatar.pRefDiGraph  # R2
        pHead = chromDatar.pHead
        pmergedPeakDict = chromDatar.pmergedPeakDict
        
        nRefDiGraph = chromDatar.nRefDiGraph
        nHead = chromDatar.nHead
        nmergedPeakDict = chromDatar.nmergedPeakDict
        
        peakPairs = chromDatar.peakPairs
        
        # data-structure has been transformed and attached to taskDatar object
        # it's time to clear these unnecessary data members now
        taskDatar.pRefDiGraph = pRefDiGraph
        taskDatar.pHead = pHead
        taskDatar.pmergedPeakDict = pmergedPeakDict
        taskDatar.nRefDiGraph = nRefDiGraph
        taskDatar.nHead = nHead
        taskDatar.nmergedPeakDict = nmergedPeakDict
        taskDatar.peakPairs = peakPairs
        
        delattr(chromDatar, 'pRefDiGraph') 
        delattr(chromDatar, 'pHead')
        delattr(chromDatar, 'pmergedPeakDict')
        delattr(chromDatar, 'nRefDiGraph')
        delattr(chromDatar, 'nHead')
        delattr(chromDatar, 'nmergedPeakDict')
        delattr(chromDatar, 'peakPairs')
        
        return taskDatar

    def retriveR1SignalsForward(self, startR2, endR2, refR2DiGraph):
        #startR2 must be a node in refR2DiGraph
        #endR2 may not
        posSignalDict = defaultdict(int)
        current = startR2
        while(current and current <= endR2):
            toR1sDict = refR2DiGraph.node[current]['toR1s']
            for R1 in toR1sDict.keys():
                if(R1 < 0 ): # In the forward case, R2 is positive
                    posSignalDict[abs(R1)] += toR1sDict[R1]
            #current = refR2DiGraph.successors(current)[0]
            succs = refR2DiGraph.successors(current)
            if(succs):
                current = succs[0]
            else:
                break
        return posSignalDict  

    def retriveR2BackgroundForward(self, startR2, endR2, refR2DiGraph):
        #startR2 must be a node in refR2DiGraph
        #endR2 may not
        posSignalDict = defaultdict(int)
        current = startR2
        while(current and current <= endR2):
            degree = refR2DiGraph.node[current]['depth']
            posSignalDict[current] = degree
            #current = refR2DiGraph.successors(current)[0]
            succs = refR2DiGraph.successors(current)
            if(succs):
                current = succs[0]
            else:
                break
        return posSignalDict        
    
    def retriveR1SignalsBackward(self, startR2, endR2, refR2DiGraph):
        #endR2 must be a node in refR2DiGraph in this backward case
        #startR2 may not
        posSignalDict = defaultdict(int)
        current = endR2
        while(current and current >= startR2):
            toR1sDict = refR2DiGraph.node[current]['toR1s']
            for R1 in toR1sDict.keys():
                if(R1 > 0 ): # In the backward case, R2 is negative
                    posSignalDict[abs(R1)] += toR1sDict[R1]
            #current = refR2DiGraph.predecessors(current)[0]
            preds = refR2DiGraph.predecessors(current)
            if(preds):
                current = preds[0]
            else:
                break
        return posSignalDict
        
    def retriveR2BackgroundBackward(self, startR2, endR2, refR2DiGraph):
        #endR2 must be a node in refR2DiGraph in this backward case
        #startR2 may not
        posSignalDict = defaultdict(int)
        current = endR2
        while(current and current >= startR2):
            degree = refR2DiGraph.node[current]['depth']
            posSignalDict[current] = degree
            #current = refR2DiGraph.predecessors(current)[0]
            preds = refR2DiGraph.predecessors(current)
            if(preds):
                current = preds[0]
            else:
                break
        return posSignalDict
        
    def runBorderScanning(self, taskDatar):
        # for pair-mode, we need separate the positive strand 
        # and negative strand from the above data-structure of R1
        pPeakBorderDict, nPeakBorderDict = {}, {}
        pPeakStrengthDict, nPeakStrengthDict = {}, {}
        for pPeak, nPeak in taskDatar.peakPairs:
            pStart = int(pPeak.split('-')[1])
            nStart = int(nPeak.split('-')[1])
            pEnd = taskDatar.pmergedPeakDict[pStart].chromEnd
            nEnd = taskDatar.nmergedPeakDict[nStart].chromEnd
            
            nR1PosSignalDict = self.retriveR1SignalsForward(pStart, nEnd, taskDatar.pRefDiGraph)
            nR1Strength = sum(nR1PosSignalDict.values())
            nPeakStrengthDict[nPeak] = nR1Strength
            nR2PosBackgroundDict = self.retriveR2BackgroundBackward(pStart, nEnd, taskDatar.nRefDiGraph)
            nR1BorderSites = self.estimateBorder(nR1PosSignalDict, nR2PosBackgroundDict)
            #nPeakBorderDict[nPeak] = nR1BorderSites
            nPeakBorderDict[nPeak] = self.borderMerged(nR1BorderSites)
            
            #retrive the positive R1 depth list by scanning the negative R2 segment
            #and the corresponding positive R2 as background distribution
            pR1PosSignalDict = self.retriveR1SignalsBackward(pStart, nEnd, taskDatar.nRefDiGraph)
            pR1Strength = sum(pR1PosSignalDict.values())
            pPeakStrengthDict[pPeak] = pR1Strength
            pR2PosBackgroundDict = self.retriveR2BackgroundForward(pStart, nEnd, taskDatar.pRefDiGraph)
            pR1BorderSites = self.estimateBorder(pR1PosSignalDict, pR2PosBackgroundDict)   
            #pPeakBorderDict[pPeak] = pR1BorderSites 
            pPeakBorderDict[pPeak] = self.borderMerged(pR1BorderSites)

        taskDatar.pPeakBorderDict = pPeakBorderDict
        taskDatar.nPeakBorderDict = nPeakBorderDict
        taskDatar.pPeakStrengthDict = pPeakStrengthDict
        taskDatar.nPeakStrengthDict = nPeakStrengthDict        
        return taskDatar   
    
    
    def estimateBorder(self, posSignalDict, backgroundDict):
        borderSites = []
        
        oPosi = [ pos for pos in posSignalDict.keys() ]
        oPosi = sorted(oPosi)
        oDepth = [ posSignalDict[pos] for pos in oPosi ]
        ou = numpy.mean(oDepth)
        
        #Chernoff bound: http://en.wikipedia.org/wiki/Chernoff_bounds
        df = lambda c: 1.0*c/ou - 1
        pf = lambda d: numpy.power(numpy.exp(d)/numpy.power(1+d, 1+d), ou)
 
        bPosi  = [ pos for pos in backgroundDict.keys() ]
        bPosi = sorted(bPosi)
        bDepth = [ backgroundDict[pos] for pos in bPosi ]
        bu = numpy.mean(bDepth)
       
        #bRugSignals = [ 1.0*pos  for pos, val in backgroundDict.items() for j in range(0, val) ]
        bRugSignals = [ random.uniform(pos-0.5, pos+0.5) for pos, val in backgroundDict.items() for j in range(0, val) ]
        try:
            brugsig_kde = stats.gaussian_kde(bRugSignals)  
        except LinAlgError:
            return borderSites
           
        N = sum(oDepth)
        for i in range(0, len(oPosi)):
            pos, depth = oPosi[i], oDepth[i]
            delta = df(depth)
            if(delta > 0):
                p_chernoff = pf(delta)
                ##if(p_chernoff < 0.05):
                if(p_chernoff < self.c): #
                    p_rug = brugsig_kde.integrate_box_1d(pos-0.5, pos+0.5)  
                    prob = 1.0 - stats.binom.cdf(depth, N, p_rug)
                    if(prob < 0.00001):
                        borderSites.append((pos, pos+1, depth, p_chernoff, prob))                        
        return borderSites
    
    def borderMerged(self, borderSites): 
        #notes: all these borders are from the same strand within the same peak
        #because only these borders are possible overlapped.
        #borderSites already sorted by position ##
        #(pos, pos+1, depth, p_chernoff, prob)
        mergedBorders = []
        for item in borderSites:
            curr_start, curr_end, curr_depth = item[0], item[1], item[2]
            if(len(mergedBorders) == 0):
                mergedBorders.append( item )
                continue
            else:
                    prevBorder = mergedBorders[-1]
                    prev_start, prev_end, prev_depth = prevBorder[0], prevBorder[1], prevBorder[2]
                    if( abs(curr_start - prev_end) <= 3 ):
                        if(curr_depth > prev_depth):
                            mergedBorders[-1] = item 
                    else:
                        mergedBorders.append(item) 
        ##furthermore, to filter those borders with very lower depth
        ##compared to neighbours borders, within the same peak region    
        ##it means the complex is majority dominated by the borders with higher depth
        maxDepth = -numpy.inf
        for border in mergedBorders:
            if(border[2] > maxDepth):
                maxDepth = border[2]
        marked = [] 
        for border in mergedBorders:
            if(1.0*maxDepth/border[2] >= 4.0 or border[2] < 5): # 80% percent
                marked.append(border)
        mergedBorders = [item for item in mergedBorders if item not in marked]                 
        return mergedBorders  
        
    def borderPairing(self, taskDatar):
        pPeakBorderDict = taskDatar.pPeakBorderDict
        nPeakBorderDict = taskDatar.nPeakBorderDict
        # be aware of these borders are located within the same peak-pair region
        # at least one of borderSites should not be none
        borderPairDiGraph = networkx.DiGraph()
        for pPeak, nPeak in taskDatar.peakPairs:
            #[pos, pos+1, depth, p_z, prob]
            #pBorderSites = sorted(pPeakBorderDict[pPeak])
            #nBorderSites = sorted(nPeakBorderDict[nPeak])
            pBorderSites = sorted( [ item[0] for item in pPeakBorderDict[pPeak] ] )
            nBorderSites = sorted( [ item[0] for item in nPeakBorderDict[nPeak] ] )
            pN, nN = len(pBorderSites), len(nBorderSites)
            if(max(pN, nN) == 0):
                continue
            #cursor for positive and negative strands
            pc, nc = 0, 0 
            #find the nearest border on the other strand regarding of orientation
            for pBorder in pBorderSites:
                borderPairDiGraph.add_node('P-' + str(pBorder))
                #---------------------------------------
                #connecting borders on the same plus strand
                index = pBorderSites.index(pBorder)
                if(index >= 1):
                    currNode = 'P-' + str(pBorder)
                    prevNode = 'P-' + str(pBorderSites[index-1])
                    paradist = abs(pBorder - pBorderSites[index-1])
                    borderPairDiGraph.add_edge(currNode, prevNode, dist=paradist, type='P-P')
                #---------------------------------------
                minDist = numpy.inf
                partner = None
                while(True):
                    if(nN == 0):
                        break
                    nBorder = nBorderSites[nc]
                    dist = nBorder - pBorder  # the 0-th is the summit of border
                    if(dist > 0 and dist < minDist):
                        minDist = dist
                        partner = nBorder
                    #else:
                    #    break 
                    #------------------------
                    else:
                        if(partner):
                            break    
                    #------------------------
                    if(nc < nN - 1):                  
                        nc += 1
                    else:
                        break                    
                if(nc != 0): # go back one step
                    nc -= 1            
                if(partner):   
                    source = 'P-'+str(pBorder)
                    target = 'N-'+str(partner)               
                    borderPairDiGraph.add_edge(source, target, dist=minDist, type='P-N')
            #for the borders on the negative strand        
            for nBorder in nBorderSites:
                borderPairDiGraph.add_node('N-'+str(nBorder))
                #-----------------------------------
                #connecting borders on the same minus strand
                index = nBorderSites.index(nBorder)
                if(index >= 1):
                    currNode = 'N-' + str(nBorder)
                    prevNode = 'N-' + str(nBorderSites[index-1])
                    paradist = abs(nBorder - nBorderSites[index-1])
                    borderPairDiGraph.add_edge(currNode, prevNode, dist=paradist, type='N-N')
                #-----------------------------------
                minDist = numpy.inf
                partner = None
                while(True):
                    if(pN == 0):
                        break
                    pBorder = pBorderSites[pc]
                    dist = nBorder - pBorder
                    if(dist > 0 and dist < minDist):
                        minDist = dist
                        partner = pBorder
                    else:
                        break  # here is different from the pBorderSites
                    if(pc < pN - 1):
                        pc += 1
                    else:
                        break                    
                if(pc != 0): # go back one step ##
                    pc -= 1            
                if(partner):
                    source = 'N-'+str(nBorder)
                    target = 'P-'+str(partner)           
                    borderPairDiGraph.add_edge(source, target, dist=minDist, type='N-P')
        return borderPairDiGraph                      
     
    def outlierEdgeCutting(self, paramGraph, measure, k=2.0): 
        paramEdges = paramGraph.edges()
        edgeValList = [ measure(edge) for edge in paramEdges ]             
        edgeFlagList = [1 for i in range(0, len(edgeValList))]
        obslist, enumber = edgeValList, len(edgeValList)
        while(True):   
            valarray = numpy.array(obslist)
            mu, sigma = numpy.mean(valarray), numpy.std(valarray)
            edgeFlagList = [ 0 if val-mu > k*sigma else 1 for val in edgeValList ] 
            obslist = [ val for i,val in enumerate(edgeValList) if edgeFlagList[i] ]
            onumber = len(obslist)
            if(onumber < enumber):
                enumber = onumber
            else:
                break
        eraseEdges = [paramEdges[i] for i,flag in enumerate(edgeFlagList) if not flag]
        paramGraph.remove_edges_from(eraseEdges)        
        return paramGraph
    
    
    def postBorderProc(self, taskDatar):
        #pPeakBorderDict = taskDatar.pPeakBorderDict
        #nPeakBorderDict = taskDatar.nPeakBorderDict
        borderPairDiGraph = self.borderPairing(taskDatar)
        ef = lambda e: borderPairDiGraph.get_edge_data(*e)['dist']        
        borderPairDiGraph = self.outlierEdgeCutting(borderPairDiGraph, ef, k=self.k)
        taskDatar.borderPairDiGraph = borderPairDiGraph
        return taskDatar
        
    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        self.log.info("Starting Border Scanning-task on %s..." % (chrid) )
        taskDatar = self.configureTask(chromNode)
        taskDatar = self.runBorderScanning(taskDatar)
        
        taskDatar = self.postBorderProc(taskDatar)
        
        #----------------------------------------------
        #output border BED file
        #----------------------------------------------
        borderInfoDict = {}
        for peak in taskDatar.pPeakBorderDict.keys():
            borderSites = taskDatar.pPeakBorderDict[peak]
            peakname = 'P' + peak
            for i in range(0, len(borderSites)):
                 start, end, depth, p_chernoff, prob = borderSites[i]
                 bordername = 'BP-' + str(start)
                 borderInfoDict[bordername] = (start, end, depth, p_chernoff, peakname)
        for peak in taskDatar.nPeakBorderDict.keys():
            borderSites = taskDatar.nPeakBorderDict[peak]
            peakname = 'P' + peak
            for i in range(0, len(borderSites)):
                start, end, depth, p_chernoff, prob = borderSites[i]
                bordername = 'BN-' + str(start)
                borderInfoDict[bordername] = (start, end, depth, p_chernoff, peakname)
        # output border BED file
        # chrid    start    end    bordername    depth    strand   chernoff    peakid    compid    pnb
        # noted: pnb means numbers of plus borders and minus borders, as well as backbones
        #----------------------------------------------
        outfile = self.outDir + '/borders_on_' + chrid + '.bed'
        handle = open(outfile, 'w')
        borderPairDiGraph = taskDatar.borderPairDiGraph
        components = networkx.connected_components(borderPairDiGraph.to_undirected())
        for i in range(0, len(components)):
            comp = components[i]
            compid = chrid + '_' + 'Comp_' + str(i+1)
            pnodes = [ node for node in comp if node.__contains__('P-') ]
            nnodes = [ node for node in comp if node.__contains__('N-') ]
            p, n, b = self.getPNBackboneNum(comp, borderPairDiGraph)
            for node in comp:
                bordername = 'B' + node
                start, end, depth, p_chernoff, peakname = borderInfoDict[bordername]
                linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t'
                linestr+= bordername + '\t' + str(depth) + '\t'
                if(node.startswith('P')):
                    linestr += '+\t'
                else:
                    linestr += '-\t'                
                linestr += str(p_chernoff) + '\t' + peakname + '\t'
                linestr+=  compid + '\t' + str(p) + ',' + str(n) + ',' + str(b) + '\n'
                handle.write(linestr)
        handle.close()

        #----------------------------------------------
        #output border-graphs
        #----------------------------------------------
        borderPairDiGraph = taskDatar.borderPairDiGraph
        out = open(self.outDir + '/borderGraphs_on_' + chrid + '.json', 'w')
        json_graph.dump(borderPairDiGraph, out)
        out.close()
        '''
        testG = readwrite.json_graph.load(open('D:\\test.json')) ##
        components = networkx.connected_components(testG.to_undirected())
        singleComponents = [ item for item in components if len(item) == 1 ]
        pairedComponents = [ item for item in components if len(item) == 2 ]
        multipComponents = [ item for item in components if len(item) >= 3 ]
        # be aware of those backbones for further analysis
        '''
        
        self.log.info("Border Scanning-task has finished on %s..." % ( chrid) ) 

    def getPNBackboneNum(self, comp, hostGraph):
        pnodes = [ node for node in comp if node.__contains__('P-') ]
        nnodes = [ node for node in comp if node.__contains__('N-') ]
        backbone = 0
        for pn in pnodes:
            for nn in nnodes:
                e1 = hostGraph.get_edge_data(pn, nn)
                e2 = hostGraph.get_edge_data(nn, pn)
                if(e1 and e2):
                    backbone += 1
        return len(pnodes), len(nnodes), backbone        