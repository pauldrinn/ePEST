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

from cement import App, init_defaults, Controller, ex, CaughtSignal

from bll.exo import ExoBlocker
import multiprocessing
import bll.workbee as workbee

from context.env import env_checking
from context.opt import opt_validating

defaults = init_defaults("ePEST", "log.logging")
defaults['ePEST']['debug'] = False
timelabel = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
defaults['log.logging']['file'] = 'logs/ePEST_' + timelabel + '.log' ##

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

def my_pre_close_hook(app):
    pass

def my_post_close_hook(app):
    pass

class ExoAppBaseController(Controller):
    class Meta:
        label = 'base'
        description = 'ePEST: ChIP-exo paired-end sequencing processing toolkit, for peak-calling and border identify.'
        arguments = [
            (['input'], {'help': 'the bamfile from ChIP-exo. REQUIRED'}),
            (['-o', '--output'], {'help': 'the output folder for results. Default is the current directory.'}),
            (['-C', '--include-chroms'], {'action': 'store', 'metavar': 'chroms', 'type': str, 'nargs': '*', 'default': None, 'dest': 'included_chroms', 'help': 'the chromosomes to be included for analysis. Default=None'}),
            (['-r', '--regex'], {'action': 'store', 'metavar': 'REGEX', 'default': '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*', 'help': 'the read id pattern. The Default is setted with Solexa platform: [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*'}),
            (['-pv', '--pvalue'], {'action': 'store', 'metavar': 'float', 'type': float, 'default': 1e-5, 'help': 'the p value of statistical significant for trigging peak-calling. Default=1e-5.'}),
            (['-D', '--dedup'], {'action': 'store', 'choices': ['True', 'False'], 'default': False, 'help': 'the flag for removing those duplicated reads; only works for paired sequencing. Default=True.'}),
            (['-dt', '--dist'], {'action': 'store', 'type': int, 'metavar': 'int', 'default': 100, 'help': 'the cut-off distance between reads touched on the flow chip for optical duplicate detection. Default=100.'}),
            (['-qv', '--qvalue'], {'action': 'store', 'type': float, 'metavar': 'float', 'default': 0.001, 'help': 'the pvalue used for amplification fragments detection. Default=0.001.'}),
            (['-R', '--rscan'], {'action': 'store', 'type': int, 'metavar': 'int', 'default': 20, 'help': 'the r-scan parameter, a minimal number of R. Default=20.'}),
            (['-s', '--minstep'], {'action': 'store', 'type': int, 'metavar': 'int', 'default': 1, 'help': 'the minimal number of base pairs for r-window spanning. Default=1.'}),
            (['-S', '--maxstep'], {'action': 'store', 'type': int, 'metavar': 'int', 'default': 8, 'help': 'the maximal number of base pairs for r-window spanning. Default=8.'}),
            (['-c', '--chernoff'], {'action': 'store', 'type': float, 'metavar': 'float', 'default': 0.05, 'help': 'the probability of chernoff inequality. Default=0.05.'}),
            (['-k', '--outlier'], {'action': 'store', 'type': float, 'metavar': 'float', 'default': 2.0, 'help': 'the kth-fold std deviation for outlier detection. Default=2.0.'}),
            (['-t', '--threads'], {'action': 'store', 'metavar': 'int', 'type': int, 'default': 2, 'help': 'the number of threads enabled for parallel computation.'})
        ]

    @ex(hide=True)
    def _default(self):
        self.app.log.info("Checked good! Starting...")

        exoBlocker = ExoBlocker(self.app) 
        print(self.app.pargs)

        chromNodes = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()

        t = self.app.pargs.threads
        pargs = self.app.pargs
        plog = self.app.log
        consumers = [multiprocessing.Process(target=workbee.processWorkBee, args=(chromNodes, results, pargs, plog)) for _ in range(t)]
        for w in consumers:
            w.start()

        baseTree = exoBlocker.dataModel.baseTree
        treeNodes = baseTree.all_nodes()
        branches = [node for node in treeNodes if baseTree.level(node.identifier)==1]
        for node in branches:
            chromNodes.put(node)

        for _ in range(0, t):
            chromNodes.put(None)
        chromNodes.join()

        num_jobs = len(branches)
        while num_jobs:
            result = results.get()
            print('Result:', result)
            num_jobs -= 1

class ExoApp(App):
    class Meta:
        label = 'ePEST'
        config_defaults = defaults
        handlers = [
            ExoAppBaseController
        ]
        hooks = [
            ('pre_setup', my_pre_setup_hook),
            ('post_setup', my_post_setup_hook),
            ('pre_run', my_pre_run_hook),
            ('post_run', my_post_run_hook),
            ('pre_close', my_pre_close_hook),
            ('post_close', my_post_close_hook)
        ]

with ExoApp() as app:
    try:
        app.run()
    except CaughtSignal as e:
        if e.signum == signal.SIGTERM:
            print("Caught signal SIGTERM...")
        elif e.signum == signal.SIGINT:
            print("Caught signal SIGINT...")
