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
    'levelist': "975/1000",
    'param'   : "r/clwc",
    'time'    : "00/06/12/18",    
    'step'    : "0",   
    'date'    : "20040101/to/20040102",
    'area'    : "-11/-62/-11/-62",
    'grid'    : "1/1",
    'format'  : "netcdf",
    'target'  : "teste.nc"
    })
    
    

  
