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
import abc
from treelib import Tree


class BioBaseError(Exception):
    """
    General environment related errors.
    
    :param msg: The error message.
    
    """
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    
    def __str__(self):
        return self.msg
        
        
class BioDataModel(object):
    """A singleton pattern of BioDataModel, for storing the data model 
       of the specific BioBlocker. Only one instance during the whole
       running-time of the BioBlocker's lifetime. Several intermediate
       results can be attached to or removed from this instance, and 
       shared by the global BioBlocker members. 
       
        >>> o1 = BioDataModel.getInstance()
        >>> o2 = BioDataModel.getInstance()
        >>> o1
        <__main__.BioDataModel object at 0x7f5aa59dfe10>
        >>> o2
        <__main__.BioDataModel object at 0x7f5aa59dfe10>
        >>> o1 == o2
        True
        >>> o1 is o2
        True
    """
    _singleton = None
    
    def __new__(cls, *args, **kwargs):
        if not cls._singleton:
            cls._singleton = super(BioDataModel, cls).__new__(cls, *args, **kwargs)
        return cls._singleton
    
    def __init__(self):
        self.description = "The singleton data-model for BioBlocker"
        self.baseTree  = Tree()  #tree storing input/output data
        self.baseTree.create_node("ROOT", "root")
    
    @staticmethod    
    def getInstance():
        if not BioDataModel._singleton:
            BioDataModel._singleton = BioDataModel()
        return BioDataModel._singleton
    #getInstance = staticmethod(getInstance)
   
   
class BioBlocker(object):
    """The blocker represents one functional module hosted in the CementApp,
       and each blocker may have several BioTasks to accomplish the designated
       goal of the functional module.       
    """
    def __init__(self, context):
        self.appContext = context
        self.dataModel = BioDataModel.getInstance()
        
    def getAppContext(self):
        return self.appContext    
    
    @abc.abstractmethod
    def initDataModel(self):
        '''should be implemented by subclasses'''    
        return
        
    @abc.abstractmethod    
    def decorateTaskSuite(self):
        return   

        
class BioTaskSuite(object):
    def __init__(self, chromNode):
        self.tasks = []
        self.chromNode = chromNode
    
    def addTask(self, task):
        self.tasks.append(task)

    def runSuite(self):
        for task in self.tasks:
            task.runTask(self.chromNode)    
        return         


class TaskDatar(object):
    '''represent the data-content of each task. Additional features 
       could be attached to or remove from it. And the important thing
       is that this datar item just for temporary information manipulation.
       After the task done, this datar should be removed for memory release.
       
       taskDatar._soft_XXX = []
       taskDatar._hard_XXX = {}
       
       Those features with '_soft_' would be removed after task done.
    '''
    def __init__(self):
        pass
        
class BioTasker(object):
    """
    """
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def runTask(self, treeNode):
        '''shoule be implemented by subclasses'''    
        return