#!/usr/bin/env python3
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

import time
import signal

from cement.core import backend
from cement.core import foundation
from cement.core import hook
from cement.core import exc

from context.env import env_checking
from context.opt import opt_validating

import bll.core as core
from bll.exo import ExoBlocker
import bll.dedup as dedup
import bll.param as param

import multiprocessing
import bll.workbee as workbee

#=============================================================================
#  Set default config options
#=============================================================================
defaults = backend.defaults("exoApp", "log")
defaults['exoApp']['debug'] = False
timelabel = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
defaults['log']['file'] = 'ePEST_' + timelabel + '.log' ##

#=============================================================================
#  Create an application
#=============================================================================
app = foundation.CementApp("exoApp", config_defaults=defaults)

#=============================================================================
#  Register any framework hook functions after app creation, and before
#  app.setup()
#=============================================================================
def my_pre_setup_hook(app):
    #=========================================================================
	#  In this stage, we would like to check the environment.
    #  note: app.log is still not available at this time.
	#=========================================================================
    env_checking(app)

def my_post_setup_hook(app):
    #=========================================================================
	#  In this stage, we would like to decorate the application logger.
    #  note: the app.Logger is available, but app.pargs is not available
	#=========================================================================
    pass
    
def my_pre_run_hook(app): 
    #=========================================================================
	#  In this stage, app.pargs is still not available due to app.argv has 
    #  not been parsed.
	#========================================================================= 
    pass
  
def my_post_run_hook(app):
    #=========================================================================
	#  In this stage, all cement-framework related configuration has finished,
    #  and the app.pargs also available now. We can do options-validation in 
    #  regarding to our business logics.
	#=========================================================================
    opt_validating(app)
    
hook.register('pre_setup', my_pre_setup_hook)
hook.register('post_setup', my_post_setup_hook)
hook.register('pre_run', my_pre_run_hook)
hook.register('post_run', my_post_run_hook)

try:
    app.setup()

    app.args.description = "ePEST: ChIP-exo paired-end sequencing processing toolkit, for peak-calling and border identify."
    
    app.args.add_argument('input', metavar='--input.bam', help="the bamfile from ChIP-exo. REQUIRED")
    
    app.args.add_argument('-o', '--output', help="the output fold for results. \
                                                  Default is the current directory.")
        
    #app.args.add_argument('-mode', '--mode', choices=['single', 'pair'], default='pair',
    #                      help="The sequencing mode. Default=pair")
    
    app.args.add_argument('-r', '--regex', action='store', metavar='REGEX', \
                           default="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*",
                           help="the read id pattern. The Default is setted with Solexa platform: [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*")
                          
    app.args.add_argument('-p', '--pvalue', action='store', metavar='float', type=float, default=1e-5,
                          help="the p value of statistical significant for trigging peak-calling. Default=1e-5.")
                 
    app.args.add_argument('-D', '--dedup', action='store', choices=['True', 'False'], default=False,
                          help="the flag for removing those duplicated reads; only works for paired sequencing. Default=True.")    
    
    app.args.add_argument('-d', '--dist', action='store', type=int, metavar='int', default=100,
                          help="the cut-off distance between reads touched on the flow chip for optical duplicate detection. Default=100.")
    
    app.args.add_argument('-q', '--qvalue', action='store', type=float, metavar='float', default=0.001,
                          help="the pvalue used for amplification fragments detection. Default=0.001.") 

    app.args.add_argument('-R', '--rscan', action='store', type=int, metavar='int', default=20,
                          help="the r-scan parameter, a minimal number of R. Default=20.")

    app.args.add_argument('-s', '--minstep', action='store', type=int, metavar='int', default=1,
                          help="the minimal number of base pairs for r-window spanning. Default=1.")

    app.args.add_argument('-S', '--maxstep', action='store', type=int, metavar='int', default=8,
                          help="the maximal number of base pairs for r-window spanning. Default=8.")

    app.args.add_argument('-c', '--chernoff', action='store', type=float, metavar='float', default=0.05,
                          help="the probability of chernoff inequality. Default=0.05.")

    app.args.add_argument('-k', '--outlier', action='store', type=float, metavar='float', default=2.0,
                          help="the kth-fold std deviation for outlier detection. Default=2.0.")                          
    app.args.add_argument('-t', '--threads', action='store', metavar='int', type=int, default=2,
                          help="the number of threads enabled for parallel computation.")
                          
                          
    ##                       
    app.run()
    #=========================================================================
    # after this run() called, the opt_validating() will execute before going
    # into detailed business logic layer
    #=========================================================================
    
    #app.log.info("%s: Checked good! Starting..." % (time.ctime()) )
    app.log.info("Checked good! Starting..." )
    
    #=========================================================================
    # now we can go to our specific application functionality
    # prepare for our application business (business logic layer: bll)
    #=========================================================================
    
    #=========================================================================
    # may assign the exoBlocker to a member of app, for extension in future
    #=========================================================================
    exoBlocker = ExoBlocker(app) 
    print(app.pargs)
    #'''
    chromNodes = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    # Start consumers
    t = app.pargs.threads # the number of cores ##
    pargs = app.pargs
    plog = app.log
    consumers = [ multiprocessing.Process(target=workbee.processWorkBee, args=(chromNodes, results, pargs, plog)) for i in range(t) ]
    for w in consumers:
        w.start()
    
    # Enqueue jobs
    baseTree = exoBlocker.dataModel.baseTree
    treeNodes = baseTree.all_nodes()
    branches = [ node for node in treeNodes if baseTree.level(node.identifier)==1 ]    
    for node in branches:
        chromNodes.put(node)
        
    # Add a poison pill for each consumer
    for i in range(0, t):
        chromNodes.put(None)
    # Wait for all of the nodes to finish
    chromNodes.join()
    
    # Start printing results
    num_jobs = len(branches)
    while num_jobs:
        result = results.get()
        print('Result:', result)
        num_jobs -= 1  
    #'''      
except exc.CaughtSignal as e:
    if e.signum == signal.SIGTERM:
        print("Caught signal SIGTERM...")
    elif e.signum == signal.SIGINT:
        print("Caught signal SIGINT...")
finally:
    app.close() 

"""
# for debug and test purpose ##
#python ePEST.py  -D True -p 1e-8 -R 25  -t 12 -c 0.05 -k 2.0 -o testExo5 /data/yezhq/Projects/ChIPexo/data/exo-seq/FoxA1/FoxA1.ChipEXO5.R1R2.Paired.Align.bam 

import sys
import time
import imp  

class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose
    def __enter__(self):
        self.start = time.time()
        return self
    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print ('elapsed time: %f ms' % self.msecs)

sys.argv += ["-oDebugOut"]

sys.argv += ["-DTrue"]

sys.argv += ["-p1e-8"]

sys.argv += ["/data/yezhq/Projects/ChIPexo/data/exo-seq/FoxA1/FoxA1.ChipEXO5.R1R2.Paired.Align.bam"]

with open('ePEST.py') as f:
    code = compile(f.read(), 'ePEST.py','exec')
    exec(code, None, None)

baseTree = exoBlocker.dataModel.baseTree
treeNodes = baseTree.all_nodes()
branches = [node for node in treeNodes if baseTree.level(node.identifier)==1]
import bll.dedup as dedup
pargs, plog = app.pargs, app.log
dedupTask = dedup.DeduplicateWorker(pargs, plog)
chromNode = branches[0]
dedupTask.runTask(chromNode)

import bll.peak as peak
peakTask = peak.PeakScanner(pargs, plog) 
peakTask.runTask(chromNode)

"""    

