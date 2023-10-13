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

from cement import App, init_defaults, Controller, ex
from cement.core import exc

from bll.exo import ExoBlocker
import multiprocessing
import bll.workbee as workbee

defaults = init_defaults("exoApp", "log")
defaults['exoApp']['debug'] = False
timelabel = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
defaults['log']['file'] = 'ePEST_' + timelabel + '.log' ##

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
        consumers = [multiprocessing.Process(target=workbee.processWorkBee, args=(chromNodes, results, pargs, plog)) for i in range(t)]
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

with ExoApp() as app:
    try:
        app.run()
    except exc.CaughtSignal as e:
        if e.signum == signal.SIGTERM:
            print("Caught signal SIGTERM...")
        elif e.signum == signal.SIGINT:
            print("Caught signal SIGINT...")
    finally:
        app.close()
