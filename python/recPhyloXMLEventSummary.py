#!/usr/bin/python
# -*- coding: utf-8 -*-


#########################################
##  Author:         Wandrille Duchemin  
##  Created:        12-Sept-2017         
##  Last modified:  10-Apr-2018        
##
##  
##  This script is used to extract a summary of events count per species
##  from a recPhyloXML file
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             ReconciledTree
##             ReconciledTreeIO
##             
## 
##  developped for python2.7
##
#########################################

from ReconciledTree import ReconciledTreeList
from ReconciledTreeIO import recPhyloXML_parser


import sys
import os



if __name__ == "__main__":

    help =  """
                This script is used to extract a summary of events count per species from a recPhyloXML file.

                usage : python recPhyloXMLEventSummary.py -i inputRecPhyloXML [-o outputFile] [--include.transfer.departure]

                            -i inputRecPhyloXML         : input recPhyloXML file

                            -o outputFile               : (optional) file to write in (by default stdout will be used)

                            --include.transfer.departure: (optional) if true, transfer departure (or branchingOut) events will be counted

               """


    nextKEY = None
    params = {
                "-i" : None, #: input recPhyloXML file
                "-o" : None, #: (optional) file to write in (by default stdout will be used)
                "--include.transfer.departure" : False #(optional) if true, transfer departure (or branchingOut) events will be counted
            }

    flagArgs = ["--include.transfer.departure"]



    for i in range(1,len(sys.argv)):

        if not nextKEY is None:
            params[nextKEY] = sys.argv[i]
            print "argument ",nextKEY,":", sys.argv[i]
            nextKEY = None
            continue

        if sys.argv[i] in params.keys():

            if sys.argv[i] in flagArgs:
                params[sys.argv[i]] = True
                print sys.argv[i],"flag activated"
            else:
                nextKEY = sys.argv[i]
            continue


    #print additionalArguments 

    if params["-i"] is None:
        print "error: input file name not given."
        print help
        exit(1)

    elif not os.path.isfile(params["-i"]):
            print "error: " + params["-i"] + " is not an existing file."
            exit(1)




    ## loading the data

    parser = recPhyloXML_parser()
    RTL = parser.parse(params["-i"])

    if RTL is None:
        print "an error occured while parsing the file" , params["-i"] , "."
        exit(1)


    if not isinstance(RTL,ReconciledTreeList):
        print "The file" , params["-i"] , "should contain a recPhylo object, none was detected."
        exit(1)        

    if not RTL.hasSpTree():
        print "The recPhylo object in file" , params["-i"] , "does not contain a species tree. Please add one (you can use the combineRecPhyloXMLfiles.py script to do so)."
        exit(1)

    ## getting the summary

    summary = RTL.getEventsSummary( includeTransferReception  = True , includeTransferDeparture = params["--include.transfer.departure"]  , indexBySpecies=True)

    ### now ouput
    if len(summary) == 0:
        print "no events to report."
        exit(0)

    OUT = sys.stdout

    if not params["-o"] is None:
        OUT = open( params["-o"] , "w" )

    ## setting parameters

    events = summary.values()[0].keys()
    events.sort()

    species = summary.keys()
    species.sort()

    columnWidths = []

    columnWidths.append( max( [len("species")] + [len(i) for i in species]) + 3 )

    for e in events:
        columnWidths.append( len(e) + 3 )

    def filledUp(s,l,c = " "):
        return s  + max(1,(l - len(s))) * c ##at minimum 1 space

    def lineFilledUp(S,L,c=" "):
        l = ""
        for i,s in enumerate(S):

            l += filledUp(s,L[i],c)
        return l

    ## printing first line
    OUT.write( lineFilledUp( ["species"] + events , columnWidths," ") + "\n" )

    for sp in summary.keys():

        elements = []
        elements.append(sp)
        elements += [str( summary[sp][e] ) for e in events]

        line = lineFilledUp(elements , columnWidths , " ")

        OUT.write( line + "\n" )

    OUT.close()
