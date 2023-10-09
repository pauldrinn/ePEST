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

import numpy
import collections
from bll import util

    
class InchWorm:
    STATETABLE = { 'INIT':        0, 
                   'HEADFORWARD': 1, 
                   'TAILFORWARD': 2,  
                   'RESET':       3, 
                   'FINISH':     -1 }
    Peak = collections.namedtuple('Peak', 'start, end, w, r, score, pvalue')               
    # Worm will never go backward, only forward
    
    def __init__(self, hostGraph, poslist, L, N, R=10, s=1, S=5, pvalue=0.0001):
        self.poslist = poslist
        self.depthlist = [ hostGraph.node[node]['depth'] for node in poslist ]
        self.size = len(self.poslist) 
        self.L, self.N, self.R = L, N, R
        self.s, self.S = s, S
        self.pvalue = pvalue
        self.peaks = []
        self.handlers = {}
        self.startState = None
        self.endStates = []
    
    def add_state(self, name, handler, end_state=0):
        self.handlers[name] = handler
        if end_state:
            self.endStates.append(name)
    
    def set_start(self, name):
        self.startState = name
        
    def configure(self):
        self.set_start(InchWorm.STATETABLE['INIT'])
        self.add_state(InchWorm.STATETABLE['INIT'], self.handler_INIT)
        self.add_state(InchWorm.STATETABLE['HEADFORWARD'], self.handler_HEADFORWARD)
        self.add_state(InchWorm.STATETABLE['TAILFORWARD'], self.handler_TAILFORWARD)
        self.add_state(InchWorm.STATETABLE['RESET'], self.handler_RESET)
        self.add_state(InchWorm.STATETABLE['FINISH'], self.handler_FINISH, end_state=1)
    
    def run(self):
        try:
            handler = self.handlers[self.startState]
        except:
            print("must call set_start() before run()")
        if not self.endStates:
            print("at least one state must be an end_state")
        prevState = None
        while True:
            nextState, nowState = handler(prevState)
            prevState = nowState
            if nextState in self.endStates:
                handler = self.handlers[nextState]
                handler(prevState)
                break 
            else:
                handler = self.handlers[nextState] 
    
    def calc_rParam(self, tail, head):
        # head included
        rParam = sum(self.depthlist[tail:head+1]) + 1
        return rParam
        
    def calculateScanScore(self, r, w, L, N):
        logx    = numpy.log(w) - numpy.log(L) + (1.0 + 1.0/r) * numpy.log(N)
        score  = r*logx - util.fastLogFactorial(r)
        lda = numpy.exp(score)
        prob = 1.0 - numpy.exp(-lda)
        return score, prob 
        
    def estimateHead(self, tail):
        head = tail
        depth = self.depthlist[head]
        while(head < self.size-1 and depth+1 < self.R):
            head += 1
            depth += self.depthlist[head]
        if(depth + 1 < self.R):
            head = None # not available of head
        return head    
                    
    def handler_INIT(self, prevState=None):
        self.tail = 0
        self.head = self.estimateHead(self.tail)
        if(not self.head):
            nextState = InchWorm.STATETABLE['FINISH']
            return nextState, InchWorm.STATETABLE['INIT']
        self.r = self.calc_rParam(self.tail, self.head) 
        self.w = self.poslist[self.head] - self.poslist[self.tail] + 1
        score, prob = self.calculateScanScore(self.r, self.w, self.L, self.N)
        if( prob <= self.pvalue):
            start, end = self.poslist[self.tail], self.poslist[self.head]
            nowPeak = InchWorm.Peak(start, end, self.w, self.r, score, prob)
            self.peaks.append(nowPeak)
            nextState = InchWorm.STATETABLE['HEADFORWARD']
        else:
            nextState = InchWorm.STATETABLE['RESET']
        return nextState, InchWorm.STATETABLE['INIT']
    
    def handler_RESET(self, prevState): 
        #jump to next walking journey
        if(prevState == InchWorm.STATETABLE['TAILFORWARD']): 
            tail = min(self.head + 1, self.size - 1)
        else:
            tail = min(self.tail + 1, self.size - 1)
        head = self.estimateHead(tail)
        if(not head):
            nextState = InchWorm.STATETABLE['FINISH']
            return nextState, InchWorm.STATETABLE['RESET']    
        self.head, self.tail = head, tail
        self.r = self.calc_rParam(self.tail, self.head) 
        self.w = self.poslist[head] - self.poslist[tail] + 1 
        score, prob = self.calculateScanScore(self.r, self.w, self.L, self.N)
        if(prob <= self.pvalue):
            start, end = self.poslist[self.tail], self.poslist[self.head]
            nowPeak = InchWorm.Peak(start, end, self.w, self.r, score, prob)
            status = self.updateLastestPeak(nowPeak)
            if(status == "FULLAPPEND" or status == "NULLAPPEND" or status == "REPLACES"):
                nextState = InchWorm.STATETABLE['HEADFORWARD'] 
            else: #(status == "WEAKPEAK"):
                nextState = InchWorm.STATETABLE['RESET']
        else:
            nextState = InchWorm.STATETABLE['RESET']                
        return nextState, InchWorm.STATETABLE['RESET']         
    
    def handler_FINISH(self, prevState=None):
        #for filtering
        badPeaks = []
        for peak in self.peaks:
            if(peak.w <= 5): #minimal size of peak
                badPeaks.append(peak)
        self.peaks = [peak for peak in self.peaks if peak not in badPeaks]
        return       
                
    def handler_HEADFORWARD(self, prevState=None):
        nextState = None    
        for step in range(self.s, self.S+1):  #[s, S], closed region
            head = self.head + step
            if(self.head + step >= self.size):
                nextState = InchWorm.STATETABLE['FINISH']
                break
            r = self.calc_rParam(self.tail, head)
            w = self.poslist[head] - self.poslist[self.tail] + 1      
            score, prob = self.calculateScanScore(r, w, self.L, self.N)
            if( prob <= self.pvalue):
                start, end = self.poslist[self.tail], self.poslist[head]
                nowPeak = InchWorm.Peak(start, end, w, r, score, prob)
                status = self.updateLastestPeak(nowPeak)
                if(status in ["FULLAPPEND", "NULLAPPEND", "REPLACED"]):
                    self.r = r
                    self.head = head
                    self.w = w
                    nextState = InchWorm.STATETABLE['HEADFORWARD']  
                    break
        if(nextState==None):            
            if(prevState in [ InchWorm.STATETABLE['INIT'], InchWorm.STATETABLE['RESET'] ]):
                nextState = InchWorm.STATETABLE['RESET']  
            else: #(prevState == InchWorm.STATETABLE['HEADFORWARD']):
                nextState = InchWorm.STATETABLE['TAILFORWARD']               
        return nextState, InchWorm.STATETABLE['HEADFORWARD']  
  
    def handler_TAILFORWARD(self, prevState=None):
        nextState = None
        for step in range(self.s, self.S+1):
            tail = min(self.tail + step, self.size - 1)
            head = self.estimateHead(tail)
            if(not head or head > self.head):
                break          
            r = self.calc_rParam(tail, self.head)
            w = self.poslist[self.head] - self.poslist[tail] + 1          
            score, prob = self.calculateScanScore(r, w, self.L, self.N)
            if( prob <= self.pvalue):
                start, end = self.poslist[tail], self.poslist[self.head]
                nowPeak = InchWorm.Peak(start, end, w, r, score, prob)
                status = self.updateLastestPeak(nowPeak)
                if(status in ["FULLAPPEND", "NULLAPPEND", "REPLACED"]):
                    self.r = r
                    self.tail = tail
                    self.w = w 
                    nextState = InchWorm.STATETABLE['HEADFORWARD']
                    break 
        if(nextState==None):
            nextState = InchWorm.STATETABLE['RESET']                     
        return nextState,  InchWorm.STATETABLE['TAILFORWARD']                  
        
    def updateLastestPeak(self, comePeak):
        status = None
        if(len(self.peaks) == 0):
            self.peaks.append(comePeak)  
            status = "NULLAPPEND"
        else:
            prePeak = self.peaks[-1]
            if( comePeak.start <= prePeak.end ):   #start position < end position
            # means the two peaks are overlapped 
                if( comePeak.score < prePeak.score ):  #comparing their scores
                    self.peaks[-1] = comePeak
                    status = "REPLACED"
                else:
                    status = "WEAKPEAK"
            else: # not overlap
                self.peaks.append(comePeak)
                status = "FULLAPPEND"                
        return status