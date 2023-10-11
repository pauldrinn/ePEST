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
import os
import pysam
import bll.util as util

"""Business environment-checking module."""

class ExoOptError(Exception):
    """
    General options related errors.

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
def opt_validating(app):
    """
    validating the options.
    :param app: the underlying application.
    :raises: context.env.ExoEnvError
    """

    optFlag = True

    # the format of input file should be BAM
    try:
        bamfile = app.pargs.input
        #app.log.debug("Bamfile:" + bamfile)
        #parsing the header information
        gSizeDict = util.getChromSizeDict(bamfile)
        if (len(gSizeDict)==0):
            optFlag = False
        #testing one record...
        bamhandle = pysam.Samfile(bamfile,'rb')
        record = next(bamhandle)
        bamhandle.close()
    except:
        optFlag = False
    finally:
        if (not optFlag):
            raise ExoOptError(
                              "OptionError: loading the input bamfile, please check the correct BAM format! " +
                              "And be aware the header information that should contain the size of each chromosome."
                             )
	
    # checking the output directory
    try:
        if app.pargs.output:
            output = app.pargs.output
            if(not os.path.exists(output)):
                tmp = os.mkdir(output)
        else:
            output = os.getcwd()
    except:
        optFlag = False
    finally:
        if(not optFlag):
            raise ExoOptError(
                              "OptionError: The directory name for output or volume label syntax is incorrect!"
                             )
