#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 08:22:31 2015

@author: acorreia
"""

from ecmwfapi import ECMWFDataServer
  
server = ECMWFDataServer()
  


 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20070101/to/20071231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2007.nc"
    })



 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20080101/to/20081231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2008.nc"
    })
    

 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20090101/to/20091231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2009.nc"
    })
    

 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20100101/to/20101231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2010.nc"
    })
    


 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20110101/to/20111231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2011.nc"
    })
    

 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20120101/to/20121231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2012.nc"
    })
    

 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20130101/to/20131231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2013.nc"
    })
    


 
server.retrieve({
    'repres'  : "SH",
    'dataset' : "interim",
    'expver'  : "0001",
    'stream'  : "OPER",
    'class'   : "EI",
    'type'    : "AN",
    'levtype' : "PL",
    'resol'   : "AV",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20140101/to/20141231",
    'area'    : "-3/-60/-3/-60",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "erainterim-MN-2014.nc"
    })
    









    
