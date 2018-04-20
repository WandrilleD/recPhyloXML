#!/usr/bin/python
# -*- coding: utf-8 -*-


#########################################
##  Author:         Wandrille Duchemin  
##  Created:        20-Sept-2017         
##  Last modified:  10-Apr-2018        
##
##  
##  
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

from ReconciledTree import ReconciledTreeList, ReconciledTree, RecEvent
from ReconciledTreeIO import recPhyloXML_parser


import sys
import os


def MakeLossIndependentNode( node , LossIndex , lostSpecies = "", lostTS = None, lostAdditional = {} , keptChildNameSuffix = ".c"):
    """
    *modifies node*

    Takes:
         - node (ReconciledTree): reconciled node where the *Loss event occurs
         - LossIndex (int): index of the speciationLoss or branchingOutLoss event
         - lostSpecies (str) [default = ""] : species of the loss
         - lostTS (int) [default = None]: timeSlice is the loss
         - lostAdditional [default = {}]: additional information to give to the new loss event
         - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss 
    """



    # 1. create the loss child

    lossNode = ReconciledTree()
    lossNode.addEvent( RecEvent("loss" , lostSpecies, ts= lostTS , additionnalInfo = lostAdditional) )
    lossNode.name="LOSS"

    # 2. create the kept child

    keptNode = ReconciledTree()
    keptNode.name= node.name + keptChildNameSuffix

    while len(node.eventRecs) >  (LossIndex+1):
        #print LossIndex, LossIndex+1 , len(node.eventRecs)
        keptNode.addEvent( node.popEvent(LossIndex+1) )

    # 3. link children to kept child

    while len(node.children) > 0:
        c = node.children[0]
        c.detach()
        keptNode.add_child(c)


    # 4. branching loss and kept to original node

    node.add_child(lossNode)
    node.add_child(keptNode)

    # 5. editing the event

    e = node.eventRecs[LossIndex].eventCode
    node.eventRecs[LossIndex].eventCode = e[:-1]

    return

def ConvertRTtoLossIndepVersion(RT , speciesTree = None, keptChildNameSuffix = ".c"):
    """
    *modifies RT*
    *RECURSIVE*

    Takes:
        - RT (ReconciledTree): reconciled tree or subtree to convert
        - speciesTree (ete3.Tree) [default = None] : species tree
        - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss 
    """

    for i, e in enumerate(RT.eventRecs):

        if len(e.eventCode)>1 and e.eventCode.endswith("L"):

            species = ""

            if not speciesTree is None:
                species = RT.getLostSpecies( i , speciesTree)
                

            lostTS = e.timeSlice
            if not lostTS is None:
                lostTS += 1

            MakeLossIndependentNode( RT , i , lostSpecies = species, lostTS = lostTS, keptChildNameSuffix = keptChildNameSuffix)

    for c in RT.children:
        ConvertRTtoLossIndepVersion(c , speciesTree , keptChildNameSuffix)

    return




if __name__ == "__main__":

    help =  """
                This script is used to extract a summary of events count per species from a recPhyloXML file.

                usage : python .py -i inputRecPhyloXML [-o outputFile]

                            -i inputRecPhyloXML         : input recPhyloXML file
                            -o outputFile               : (optional) file to write in (by default stdout will be used)

               """


    nextKEY = None
    params = {
                "-i" : None, #: input recPhyloXML file
                "-o" : None, #: (optional) file to write in (by default stdout will be used)

            }

    flagArgs = []


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
    RTL = parser.parse(params["-i"], 0)

    if RTL is None:
        print "an error occured while parsing the file" , params["-i"] , "."
        exit(1)


    if not isinstance(RTL,ReconciledTreeList):
        print "The file" , params["-i"] , "should contain a recPhylo object, none was detected."
        exit(1)        

    #if not RTL.hasSpTree():
    #    print "The recPhylo object in file" , params["-i"] , "does not contain a species tree. Please add one (you can use the combineRecPhyloXMLfiles.py script to do so)."
    #    exit(1)


    ## converting


    for RT in RTL.recTrees:
        #print RT
        ConvertRTtoLossIndepVersion(RT , speciesTree = RTL.spTree, keptChildNameSuffix = ".c")
        #print RT

    ## writing

    OUT = sys.stdout

    if not params["-o"] is None:
        OUT = open( params["-o"] , "w" )

    lines = RTL.getRecPhyloXMLLines()
    for l in lines:
        OUT.write(l+"\n")
    OUT.close()