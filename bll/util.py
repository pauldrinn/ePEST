import pysam
import numpy
import networkx
    
def loadBAMtoReadMapDict(bamfile):         
    #strand should be True or False
    samrecords = pysam.Samfile(bamfile,'rb')  #sorted bamfile; if paired, should be properly paired
    chrFragmentMapDict = {}  #organized by chromosome
    #count = 0
    for record in samrecords:
        #if(count >= 8000000):
        #    break
        #count += 1
        chrid = samrecords.getrname(record.tid)
        readid = record.qname.split()[0]
        start = record.pos  #on the reference, 0-based
        #end   = start + record.rlen - 1 # [start,end]
        end   = record.aend - 1 # [start,end]
        #'-' means negative strand, the pos is the 5'-end of reads
        pos   = -end if record.is_reverse else start
        if(chrid not in chrFragmentMapDict):
	        chrFragmentMapDict[chrid] = {}
        if(readid not in chrFragmentMapDict[chrid]):
            chrFragmentMapDict[chrid][readid] = {}
        if(record.is_read1):
            #R1 means exo-side of fragment
            chrFragmentMapDict[chrid][readid]['R1'] = pos
        elif(record.is_read2):
            #R2 means sonication-side of fragment
            chrFragmentMapDict[chrid][readid]['R2'] = pos
        else:
            #Single End
            chrFragmentMapDict[chrid][readid]['SE'] = pos  
    return chrFragmentMapDict
    
 
def getChromSizeDict(bamfile):
    chrSizeDict = {}
    bamhandle = pysam.Samfile(bamfile, 'rb')
    header = bamhandle.header
    if('SQ' in header):
        for item in header['SQ']:
            chrSizeDict[item['SN']] = item['LN']
    bamhandle.close()
    return chrSizeDict 

_LOGFACTORIALDICT = {0:0}
_LOGFACTORIALDICT = { i:sum(map(numpy.log, range(1,i+1))) for i in range(1,201) }
def fastLogFactorial(r):
    if(r <= 200):
        value = _LOGFACTORIALDICT[r]
    else:
        #Stirling's approximation (http://en.wikipedia.org/wiki/Stirling_formula)
        value = 0.5*numpy.log(2.0*numpy.pi*r) + r*(numpy.log(r) - 1) + 0.083333/r
    return value  
    
def buildLociLinkGraph(fragmentMapDict, strand, primer):
    refLociDiGraph = networkx.DiGraph()
    for fid in fragmentMapDict.keys():
        endsite = fragmentMapDict[fid][primer]
        flag = False
        if( strand == '-' and endsite < 0 ): #negative strand
            flag = True
        if( strand == '+' and endsite > 0 ): #positive strand
            flag = True
        if(flag):
            node = abs(endsite)
            if( not refLociDiGraph.has_node( node ) ):
                refLociDiGraph.add_node(node, depth=0)
            refLociDiGraph.node[node]['depth'] += 1		
    nodes = refLociDiGraph.nodes()
    nodes = sorted(nodes)
    for i in range(0, len(nodes)-1):
        refLociDiGraph.add_edge(nodes[i], nodes[i+1])
    return nodes, refLociDiGraph    
    

def getFragSizeOfPairedLib(chrFragmentMapDict):
    #those paired reads has been validated as properly paired
    #so we don't need to valid their position here
    fragmentSizes = []
    for chrid in chrFragmentMapDict.keys():
        fragmentDict = chrFragmentMapDict[chrid]	
        for fid in fragmentDict.keys():
            R1 = fragmentDict[fid]['R1']
            R2 = fragmentDict[fid]['R2']
            fragmentSizes.append(abs(abs(R1)-abs(R2)))	
    fragMean = numpy.mean(fragmentSizes)
    fragStd  = numpy.std(fragmentSizes)
    return fragMean, fragStd 
    
#how we could evaluate the fragment size from single-end library?    
def getFSizefromSingleLib(chrFragmentMapDict):
    fragmentSizes = []
    for chrid in chrFragmentMapDict.keys():
        fragmentDict = chrFragmentMapDict[chrid]
        for fid in fragmentDict.keys():
            pos = fragmentDict[fid]['SE']  
            if( pos < 0 ): # negative strand
                pass
            else:  # positive strand
                pass            
    return None

