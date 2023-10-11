"""Module Description

Copyright (c) 2010,2011 Zhenqing Ye <yez@uthscsa.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Zhenqing Ye
@contact: yez@uthscsa.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import bll.core as core
import bll.util as util
        
class ChromDatar(object):
    '''represent the data-content of each chromosome. Additional features 
       could be attached to or remove from it.
    '''
    def __init__(self, chrid):
        self.tag = chrid
            
class ExoBlocker(core.BioBlocker):
    """The Blocker encapsulate those routines for ChIP-exo analysis."""

    def __init__(self, context):
        super().__init__(context)
        self.log = self.appContext.log
        self.pargs = self.appContext.pargs
        self.configureBlocker()
        
    def configureBlocker(self):
        '''In this stage, we would initialize the data-model and
           set-up the suitable TaskSuite according to the pargs.'''
        self.initDataModel() 
    
    def initDataModel(self):
        #self.log.info("%s: Initializing data for model..." % (time.ctime()))
        self.log.info("Initializing data for model..." )
        input = self.pargs.input
        included_chroms = self.pargs.included_chroms
        chrSizeDict = util.getChromSizeDict(input)
        chrFragmentMapDict = util.loadBAMtoReadMapDict(input, included_chroms)
        baseTree = self.dataModel.baseTree
        for chr in chrFragmentMapDict.keys():
            chromDatar = ChromDatar(chr)
            fragmentMapDict = chrFragmentMapDict[chr]
            chromDatar.map = fragmentMapDict
            chromDatar.size= chrSizeDict[chr]
            baseTree.create_node(chr.upper(), chr.lower(), parent='root', data=chromDatar)  
        return
     
    """     
    def decorateTaskSuite(self, node):
        taskSuite = BioTaskSuite(node)
        # pair-mode
        if(self.pargs.mode == 'pair'):
            self.log.info("%s: Preparing for PairMode... " % (time.ctime())) 
            #----------------Dedup task section---------------------#
            if(self.pargs.dedup):
                dedupTask = dedup.DeduplicateWorker(self.appContext)
                taskSuite.addTask(dedupTask)
            else:
                #prepare the initial data structure...working...
                pass
                
            #----------------Peak calling section-------------------#   
            #peakScanner = peak.PeakScanner(self) 
            #taskSuite.addTask(peakScanner)      
            
            #----------------Border Scanning section----------------#   
            #borderScanner = border.BorderScanner(self) 
            #taskSuite.addTask(borderScanner)      
        # single-mode    
        else: 
            self.log.info("%s: Preparing for SingleMode... " % (time.ctime())) 
            pass
        return taskSuite    


    def processBranch(self, node):
        chrid = node.identifier
        print("branch " + chrid + " start")
        def startBee(node):
            self.taskSuite.runSuite(node)
        processBee = Process(target=startBee, args=(node,), name=chrid)   
        processBee.start()        
        processBee.join()  ##?? how about comments this line???
        print("branch " + chrid + " finished")
        return chrid
        
    def startBlocker(self):    
        self.log.info("%s: TaskSuite has been set-up, about to start... " % (time.ctime())) 
        baseTree = self.dataModel.baseTree
        treeNodes = baseTree.all_nodes()
        branches = [ node for node in treeNodes if baseTree.level(node.identifier)==1 ]      
        
        #Here we combined the threading with multiprocessing for parallelism,
        #within the producer thread, we create several processes
        t = self.pargs.threads # the number of cores
        que = Queue(t)
        def producer(que, branches):
            for treenode in branches:
                que.put(treenode, True)
        finished = []
        def consumer(que, n, t):
            nodeList = []
            print("Start consumer")
            while(len(finished) < n):
                node = que.get(True)
                nodeList.append(node)
                k = len(nodeList)
                if( k > 0 and (k == t or k == n - len(finished) ) ):
                    print("Start ThreadPool")
                    with ThreadPoolExecutor(max_workers=t) as executor:
                        futures = [executor.submit(self.processBranch, bee) for bee in nodeList]    
                        for f in as_completed(futures):
                            finished.append(f.result())
                            print("future  result:" + f.result())                    
                        nodeList=[]
        n = len(branches)               
        prod_thread = threading.Thread(target=producer, args=(que, branches))
        cons_thread = threading.Thread(target=consumer, args=(que, n, t))
        prod_thread.start()
        cons_thread.start()
        prod_thread.join()
        cons_thread.join()    
        print(finished)
    """
    

     
