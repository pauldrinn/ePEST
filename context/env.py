"""Module Description

Copyright (c) 2010,2011 Zhenqing Ye <Zhenqing.Ye@osumc.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Zhenqing Ye
@contact: Zhenqing.Ye@osumc.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import os

"""Business environment-checking module."""

class ExoEnvError(Exception):
    """
    General environment related errors.
    
    :param msg: The error message.
    
    """
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    
    def __str__(self):
        return self.msg
        
# ------------------------------------
# Misc functions
# ------------------------------------
def env_checking ( app ):
    """Checking the running environment required by the app.
    
    :param app: The underlying application.
    
    :raises: context.env.ExoEnvError
    """
    envFlag = True
        
    #checking for the python version
    if sys.version_info[0] != 3 or sys.version_info[1] < 3:
        raise ExoEnvError(
                          "EnvError: You are using the Python %s.%s, " % (sys.version_info[0], \
                          sys.version_info[1]) + \
                          "But the %s requires Python 3.3 or above!" % app._meta.label
                         )

    #checking for scipy library
    try:
        import scipy
        scipy_version_info = scipy.__version__
        scipy_version_info = scipy_version_info.split('.')
        if int(scipy_version_info[1]) != 13 or int(scipy_version_info[2]) < 2:
            envFlag = False
    except ImportError:
        envFlag = False
    finally:
        if(not envFlag):
            raise ExoEnvError(
                              "EnvError: checking the scipy module, the %s requires " \
                              % app._meta.label + \
                              "scipy-0.13.2 or above. Please see http://www.scipy.org."
                             )

                             
  
    #checking for numpy library
    try:
        import numpy
        numpy_version_info = numpy.__version__
        numpy_version_info = numpy_version_info.split('.')
        if int(numpy_version_info[0]) < 1 or int(numpy_version_info[1]) < 8:
            envFlag = False
    except ImportError:
        envFlag = False
    finally:
        if(not envFlag):
            raise ExoEnvError(
                              "EnvError: checking the numpy module, the %s requires " 
                              % app._meta.label + 
                              "numpy-1.8.0 or above. Please see http://www.numpy.org."
                             )
                             
    #checking for pysam library
    try:
        import pysam
        pysam_version_info = pysam.__version__
        pysam_version_info = pysam_version_info.split('.')
        if int(pysam_version_info[1]) < 7 or int(pysam_version_info[2]) < 7:
            envFlag = False
    except ImportError:
        envFlag = False
    finally:
        if(not envFlag):
            raise ExoEnvError(
                              "EnvError: checking the pysam module, the %s requires " 
                              % app._meta.label + 
                              "pysam-0.7.7 or above. Please see http://code.google.com/p/pysam/"
                             )                         
                             


