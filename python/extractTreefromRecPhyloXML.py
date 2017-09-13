#!/usr/bin/python
# -*- coding: utf-8 -*-


#########################################
##  Author:         Wandrille Duchemin  
##  Created:        11-Sept-2017         
##  Last modified:  11-Sept-2017        
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

                usage : python combineRecPhyloXMLfiles.py -i inputRecPhyloXML  p1 [p2 p3 ...  -p positionFile -o outputRecPhyloXML]

                            -i inputRecPhyloXML         : input recPhyloXML file

                            p1 p2 p3 ...                :  a set of positions in the recPhyloXML object (integers)

                            -o outputRecPhyloXML        : (optional) file to write in (by default stdout will be used)
                            -p positionFile             : (optional) a file containing positions (one per line)
                            --include.species.tree      : (optional) whether the species tree should be included in the output file

                example : python combineRecPhyloXMLfiles.py -i file.xml 1 4 5 
                          will extract trees at index 1, 4 and 5 in the file file.xml and write them in stdout

               """


    nextKEY = None
    params = {
                "-i" : None, #: input recPhyloXML file
                "-o" : None, #: (optional) file to write in (by default stdout will be used)
                "-p" : None, #: (optional) a file containing positions (one per line)
                "--include.species.tree" : False #: (optional) whether the species tree should be included in the output file            }
            }

    flagArgs = ["--include.species.tree"]

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


    ## getting more files

    if not params["-p"] is None:

        if not os.path.isfile(params["-p"]):
            print "error: " + params["-p"] + " is not an existing file."
            exit(1)

        IN = open(params["-p"], "r")
        l = IN.readline()
        while l !="":
            additionalArguments.append(l.strip())
            l = IN.readline()

        IN.close()


    if len(additionalArguments) == 0:
        print 'no position to extract specified. Please specify at least one.'
        print help
        exit(1)

    ## loading the data

    parser = recPhyloXML_parser()
    RTL = parser.parse(params["-i"])

    if RTL is None:
        print "an error occured while parsing the file" , params["-i"] , "."
        exit(1)

    newRTL = ReconciledTreeList()

    if params["--include.species.tree"]:
        newRTL.setSpTree( RTL.spTree ) 


    for strP in additionalArguments:
        intP = None
        try:
            intP =int(strP)

        except:
            print "error : position" , strP, "does not convert into an integer."
            exit(1)

        if intP >= len(RTL):
           print "error : position" , strP, "is too big (There are", len(RTL) ,"reconciled trees in the provided input file)."
           exit(1)

        newRTL.append( RTL[intP] )


    ### now ouput
    
    lines = newRTL.getRecPhyloXMLLines()


    OUT = sys.stdout

    if not params["-o"] is None:
        OUT = open( params["-o"] , "w" )

    for l in lines:
        OUT.write( l + "\n" )

    OUT.close()
