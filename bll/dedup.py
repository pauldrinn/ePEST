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

import re
import os
import numpy
import random
import networkx
import collections
from scipy import stats

from bll.core import BioTasker
from bll.core import TaskDatar

class DeduplicateWorker(BioTasker):

    def __init__(self, pargs, plog):
        self.pargs = pargs
        self.log = plog

        if (self.pargs.regex==None):
            self.REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*"
        else:
            self.REGEX = self.pargs.regex

        self.optDistance = self.pargs.dist
        self.qvalue = self.pargs.qvalue
        
        output = self.pargs.output + '/Dedup'
        if (not os.path.exists(self.pargs.output + '/Dedup')):
            tmp = os.makedirs(output, exist_ok=True)

    def configureTask(self, chromNode):
        taskDatar = TaskDatar()    
        chromDatar = chromNode.data
        fragMapDict = chromDatar.map
        #directed graph, nodes are abs(R2)
        lociDiGraph = networkx.DiGraph()
        for fid in fragMapDict.keys():
            R1 = fragMapDict[fid]['R1'] #to +/-
            R2 = abs(fragMapDict[fid]['R2'])
            if (not lociDiGraph.has_node(R2)):
                lociDiGraph.add_node(R2, toR1s=collections.defaultdict(list))
            lociDiGraph.nodes[R2]['toR1s'][R1].append(fid)
        nodes = sorted(lociDiGraph.nodes())
        for i in range(0, len(nodes)-1):
            lociDiGraph.add_edge(nodes[i], nodes[i+1])  

        taskDatar._soft_sortedPoints = nodes            
        taskDatar._soft_lociDiGraph = lociDiGraph

        delattr(chromDatar, 'map') # we do not need it more, to save space

        #---------------------------------------------------
        # calculate the EXTENTION for amplification section
        #---------------------------------------------------
        chromSize = chromDatar.size
        taskDatar._soft_EXTENTION = numpy.ceil(1.0*chromSize/len(nodes))
        return taskDatar         
        
    def getLocation(self, readid):
        m = re.search(self.REGEX, readid)
        if (m is not None):
            return list(map(int,m.groups()))
        else:
            return None
    
    def getLinkedR1Count(self, R2, lociDiGraph):
        toR1sDict = lociDiGraph.nodes[R2]['toR1s']
        localPoints = []
        for R1 in map(abs, toR1sDict.keys()):
            fragments = toR1sDict.get(R1,[])
            fragments.extend(toR1sDict.get(-R1,[]))
            size = len(fragments)
            localPoints.append(size)
        return numpy.median(localPoints)
        
    def buildLocalPoints(self, R2, taskDatar):
        lociDiGraph = taskDatar._soft_lociDiGraph
        pointNode = R2
        extendedR2List = []
        # forward step
        while (True):
            prevs = list(lociDiGraph.predecessors(pointNode))
            if (prevs):
                prev = prevs[0]
                if (R2 - prev < taskDatar._soft_EXTENTION):
                    extendedR2List.append(prev)
                else:
                    break
                pointNode = prev    
            else:
                break 
        extendedR2List.reverse()                
        # backward step
        pointNode = R2 # reset
        while (True):
            nexts = list(lociDiGraph.successors(pointNode))
            if (nexts):
                next = nexts[0]
                if (next - R2 < taskDatar._soft_EXTENTION):
                    extendedR2List.append(next)
                else:
                    break
                pointNode = next
            else:
                break            
        # count R1 points when given R2      
        localPoints = [self.getLinkedR1Count(point, lociDiGraph) for point in extendedR2List]
        return localPoints
        
    def findDuplicates(self, taskDatar, chrid):
        amplifyFragments, opticalFragments = [], []
        lociDiGraph = taskDatar._soft_lociDiGraph
        sortedPoints = taskDatar._soft_sortedPoints
        for R2 in sortedPoints:
            toR1sDict = lociDiGraph.nodes[R2]['toR1s']
            localPoints = None # prepare firstly, filled when necessary
            for R1 in map(abs, toR1sDict.keys()):
                fragments = toR1sDict.get(R1,[])
                fragments.extend(toR1sDict.get(-R1,[]))
                size = len(fragments)
                
                #------------------------------------------------------
                # Amplification section
                #------------------------------------------------------
                # filtering duplicate reads caused by PCR amplification, 
                # based on R2 reads distribution. In ChIPexo, we can't 
                # use R1 because of it is exo-digesting site. And for R2, 
                # we could assume that has an even distribution around. 
                # those amplification reads would introduce an abrupt can 
                # be detected by a statistical scanning.
                #------------------------------------------------------
                if (size >= 5 ): 
                # here we set the baseline (5) due to  
                # 1.0 - stats.skellam.cdf(5, 1, 1) = 0.00025
                # less than 5, the qvalue will become in-significant                
                    if (localPoints is None):
                        localPoints = self.buildLocalPoints(R2, taskDatar)      
                    if (len(localPoints) == 0):
                        #if there are no neighbouring reads, it means this fragment 
                        #itself is too sparse here, and be a rare event
                        #amplifyFragments.extend(fragments)
                        lociLambda = 1.0 # suppose the pseudo count
                    else:                        
                        lociLambda = numpy.mean(localPoints)
                    #go to Skellam distribution routine    
                    k = numpy.ceil(size - lociLambda)
                    if (k <= 0):
                        continue
                    else:
                        prb = 1.0 - stats.skellam.cdf(k, lociLambda, lociLambda) 
                        if (prb <= self.qvalue):
                            ampids = random.sample(fragments, int(k)-1 ) 
                            #amplifyFragments.extend(ampids)
                            for item in ampids:
                                amplifyFragments.append((R2, R1, item, k, size, lociLambda, prb))
                #------------------------------------------------------                
                
                #------------------------------------------------------
                # Optical section
                #------------------------------------------------------
                for i in range(0, size):
                    fieldi= self.getLocation(fragments[i])   #taili,xi,yi
                    if (fieldi is None or len(fieldi)!=3):
                        continue
                    for j in range(i+1,size):
                        fieldj= self.getLocation(fragments[j])  #tailj,xj,yj
                        if (fieldj is None or len(fieldj)!=3):
                            continue
                        if (fieldi[0]!=fieldj[0]): #tile
                            continue
                        if (abs(fieldi[1]-fieldj[1]) > self.optDistance): #x-location
                            continue
                        if (abs(fieldi[2]-fieldj[2]) <= self.optDistance): #y-location
                            opticalFragments.append( (R2, R1, fragments[j], fragments[i]) )
                #------------------------------------------------------               
        taskDatar._soft_amplifyFragments = amplifyFragments
        taskDatar._soft_opticalFragments = opticalFragments
        return taskDatar
    
    def cleaningTask(self, chromNode, taskDatar):
        chromDatar = chromNode.data
        #------------------------------------------
        # change the fragments to count information
        # and remove those duplicated reads
        #------------------------------------------
        R2R1_fragmentsDict = collections.defaultdict(dict)
        for record in taskDatar._soft_amplifyFragments:
            R2, R1, frag, k, size, lociLambda, prb = record
            if (R1 not in R2R1_fragmentsDict[R2]):
                R2R1_fragmentsDict[R2][R1] = set()
            R2R1_fragmentsDict[R2][R1].add(frag)
                        
        for record in taskDatar._soft_opticalFragments:
            R2, R1, optfid, matchfid = record
            if (R1 not in R2R1_fragmentsDict[R2]):
                R2R1_fragmentsDict[R2][R1] = set()
            R2R1_fragmentsDict[R2][R1].add(optfid)            
        
        lociDiGraph = taskDatar._soft_lociDiGraph
        sortedPoints = taskDatar._soft_sortedPoints
        for R2 in sortedPoints:
            dupR1FragDict = R2R1_fragmentsDict.get(R2, None)
            toR1sDict = lociDiGraph.nodes[R2]['toR1s']
            for R1 in toR1sDict.keys():
                fragments = toR1sDict[R1]
                if (dupR1FragDict and abs(R1) in dupR1FragDict):
                    dupFragSet = dupR1FragDict[abs(R1)]
                    fragments = set(fragments).difference(dupFragSet)
                toR1sDict[R1] = len(fragments)
            if (sum(toR1sDict.values()) == 0):
                prevs = list(lociDiGraph.predecessors(R2))
                nexts = list(lociDiGraph.successors(R2))
                if (prevs and nexts):
                    lociDiGraph.add_edge(prevs[0], nexts[0])
                lociDiGraph.remove_node(R2)                    
        sortedPoints = sorted(lociDiGraph.nodes())
        chromDatar.lociDiGraph = lociDiGraph
        chromDatar.sortedPoints = sortedPoints
        del(taskDatar)
        return  
        
    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        #self.log.info("%s: Starting Dedup-task on %s..." % (time.ctime(), chrid) )
        self.log.info("Starting Dedup-task on %s..." % (chrid) )
        taskDatar = self.configureTask(chromNode)

        taskDatar = self.findDuplicates(taskDatar, chrid)
        #---------------------------------------
        # output section
        #---------------------------------------
        outfile = self.pargs.output + '/Dedup/dedup_on_' + chrid + '.txt'
        handle = open(outfile, 'w')
        headline = "#-----------------------------\n"
        headline+= "# Optical Section \n"
        headline+= "#-----------------------------\n"
        headline+= "chrid\tstart\tend\toptfid\tmatchfid\n"
        handle.write(headline)
        for item in taskDatar._soft_opticalFragments:
            R2, R1, optfid, matchfid = item
            start, end = min(R2, R1), max(R2, R1)
            linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t' \
                    + optfid + '\t' + matchfid + '\n'
            handle.write(linestr)
            
        headline = "\n\n#-----------------------------\n"
        headline+= "# Amplification Section \n"
        headline+= "#-----------------------------\n"
        headline+= "chrid\tstart\tend\ampfid\tcount\tmean\tdiff\tprb\n"   
        for item in taskDatar._soft_amplifyFragments:
            R2, R1, frag, k, size, lociLambda, prb = item
            start, end = min(R2, R1), max(R2, R1)
            linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t' \
                    + frag + '\t' + str(size) +'\t' + str(lociLambda) + '\t' \
                    + str(k) + '\t' + str(prb) + '\n'
            handle.write(linestr)
        handle.close() 
        #---------------------------------------
        
        self.cleaningTask(chromNode, taskDatar)
     
        self.log.info("Dedup-task is finished on %s..." % ( chrid) )  
        self.log.info("Please refer to %s!" % (outfile) )		
    
