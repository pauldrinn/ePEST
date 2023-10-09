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

import bll.dedup as dedup
import bll.peak as peak
import bll.border as border
import bll.core as core

def processWorkBee(chromNodes, results, pargs, plog):
    while True:
        next_chrom = chromNodes.get()
        if next_chrom is None:
            #print("%s: Exiting" % proc_i)
            chromNodes.task_done()
            break
        chromId = next_chrom.identifier   
        #time.sleep(10)
        #==============================
        taskSuite = decorateTaskSuite(pargs, plog, next_chrom)
        taskSuite.runSuite()
        #==============================
        chromNodes.task_done()
        results.put(chromId) 

def decorateTaskSuite(pargs, plog, chromNode):
    taskSuite = core.BioTaskSuite(chromNode)
    # at current, we do not deal with the single-end ChIP-exo data
    # pair-mode
    #if(pargs.mode == 'pair'):
    if(True):
        #plog.info("%s: Preparing for PairMode... " % (time.ctime())) 
        #----------------Dedup task section---------------------#
        if(pargs.dedup):
            dedupTask = dedup.DeduplicateWorker(pargs, plog)
            taskSuite.addTask(dedupTask)
        else:
            #prepare the initial data structure...working...
            pass
                
        #----------------Peak calling section-------------------#   
        peakScanner = peak.PeakScanner(pargs, plog) 
        taskSuite.addTask(peakScanner)     
            
        #----------------Border Scanning section----------------#   
        borderScanner = border.BorderScanner(pargs, plog) 
        taskSuite.addTask(borderScanner)      
    # single-mode    
    else: 
        #plog.info("%s: Preparing for SingleMode... " % (time.ctime())) 
        pass
    return taskSuite             
    