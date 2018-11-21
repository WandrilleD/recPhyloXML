#!/usr/bin/python
# -*- coding: utf-8 -*-


#########################################
##  Author:         Wandrille Duchemin  
##  Created:        21-Nov-2018         
##  Last modified:  21-Nov-2018        
##
##  
##  This script is used to extract some trees from a recPhyloXML file
##  
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             ReconciledTree
##             ReconciledTreeIO
##             
## 
##  developped for python2.7
##
#########################################

from ReconciledTree import ReconciledTree, RecEvent, ReconciledTreeList
from ReconciledTreeIO import recPhyloXML_parser


import sys
import os



if __name__ == "__main__":

    help =  """
                This script is used to extract one or several trees from a recPhyloXML file (containing different reconciled gene trees).
                ( NB: positions start at index 0 )

                usage : python RecPhyloXMLtoNewick.py -i inputRecPhyloXML  [-o outputRecPhyloXML]

                            -i inputRecPhyloXML         : input recPhyloXML file
                            -p positionFile             : (optional) a file containing positions (one per line)
               """


    nextKEY = None
    params = {
                "-i" : None, #: input recPhyloXML file
                "-o" : None, #: (optional) file to write in (by default stdout will be used)
            }

    flagArgs = []

    additionalArguments = []

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
        else:
            additionalArguments.append(sys.argv[i])


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

    
    #for RT in RTL:
    #    print RT.getTreeNewick(topoOnly = True)

    ### now ouput
    

    OUT = sys.stdout

    if not params["-o"] is None:
        OUT = open( params["-o"] , "w" )


    for RT in RTL:
        OUT.write( RT.getTreeNewick(topoOnly = True) + "\n" )

    OUT.close()
