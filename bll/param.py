from bll import util

class ParameterEstimator(object):
    '''calculating some useful parameters for application.'''
    
    def __init__(self, context):
        self._params = {}
        self.context = context
        
    def getParameters(self):
        self._calcFragmentSize()
        self._calcGSizeDict()
        return self._params
        
    def _calcFragmentSize(self):
        chrFragmentMapDict = self.context.getInBaseModel()
        fragMean, fragStd = util.getFragSizeOfPairedLib(chrFragmentMapDict)
        self._params['_fragMean'] = fragMean
        self._params['_fragStd']  = fragStd	    
      
    def _calcGSizeDict(self):          
        bamfile = self.context.pargs.input  
        gsizeDict = util.getGSizeDictfromBAMHEADER(bamfile)   
        self._params['_gsize'] = gsizeDict        