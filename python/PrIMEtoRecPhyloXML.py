#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        23-Mar-2017         
##  Last modified:  10-Apr-2018        
##
##  Decribes functions to transform a reconciled tree in 
##  PrIME format into a tree in the recPhyloXML format
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             ReconciledTree , and ete3-based representation of a reconciled tree
## 
##  developped for python2.7
##
#########################################


import ete3 
from ReconciledTree import RecEvent, ReconciledTree


def setNodeAsXloss(node):
    """ transform a S into a SL event, a bro into a broL, ... """
    node.getEvents()[0].eventCode += "L"

def isXnode(node,X,i=0):
    """ checks if the event at index i (defautl = 0) of the node is a X """

    return node.getEvents()[i].eventCode == X

def isLossNode(node):
    """ checks if the first event of the node is a LOSS """
    return isXnode(node,"L")


def isTansferNode(node):
    """ checks if the first event of the node is a branchingOut """
    return isXnode(node,"bro")


def addTbEvent(node):
    """ adds a transferBack event at the beginning of a nodes chain of event"""
    sp = node.getEvents()[0].species
    tbEvent = RecEvent("Tb",sp)
    node.addEvent(tbEvent,append=False)


def PrIMEAnnotationToRecEvent(annot):
    """ returns a RecEvent OR None in case of problem """

    evtType = None
    name = annot.get("name",None)
    sp = None
    if not name is None:
        sp = name.rpartition("_")[2]


    Vtype = annot.get("VERTEXTYPE",None)
    if not Vtype is None: 
        if Vtype == "Duplication":
            evtType = "D"
        elif Vtype == "Speciation":
            evtType = "S" ## may be a SL. we check this later as this info is in the child
        elif Vtype == "Transfer":
            evtType = "bro" ## we will have to add the transferBack later to the child
        elif Vtype == "Loss":
            evtType = "L"
        elif Vtype == "Leaf":
            evtType = "C"
        else:
            print "UNKNOWN EVENT TYPE", Vtype,"!"
            
    if ( sp is None ) or (evtType is None):
        print "error when trying to assign event. data:", annot
        return None

    addInfo = {}
    if evtType == "C":
        addInfo["geneName"] = name

    return RecEvent(evtType,sp, additionnalInfo = addInfo)


def parsePrIMEAnnotations(info):

    nameBL,j,comment = info.partition("[")

    if comment.endswith("]"):
        comment = comment[:-1]

    name, j, BL = nameBL.partition(":")

    commentList = comment.split(" ")
    
    result = {"name" : name }
    if not BL == "":
        result["BL"] = BL

    for e in commentList:
        if e.startswith("&"):
            continue
        else:
            k,j,v = e.partition("=")
            result[k] = v

    return result


def read_leaf(s):
    annot = parsePrIMEAnnotations(s)

    evt = PrIMEAnnotationToRecEvent(annot)

    leafNode = ReconciledTree()
    leafNode.setName(annot["name"])

    if annot.has_key("BL"):
        leafNode.dist = annot["BL"]

    leafNode.addEvent(evt)

    #print leafNode.getTreeRecPhyloXML()

    return leafNode

def readPrIMEParenthesis(s):
    
    if not s.startswith("("):
        print "error, should begin with a parenthesis"
        return


    #### parsiing parentheses

    i = 1    
    new_n = ""
    info = "" 

    childrenNodes = []


    inComment = False
    beginCommentChar = "["
    endCommentChar = "]"

    while i < len(s):

        if s[i] == beginCommentChar:
            inComment = True
        elif s[i] == endCommentChar:
            inComment = False

        #print i , s[i] , inComment

        if inComment:
            new_n += s[i]
        else:
            if s[i] == "(":##new parenthesis
                offset , child = readPrIMEParenthesis(s[i:])
                i += offset
                childrenNodes.append(child)
            #print "added new internal."

            if s[i] in [")",","]:
                if not new_n == "":
                    #print "added leaf",new_n,"parent: current node"
                    childrenNodes.append( read_leaf(new_n) )
                    #info = new_n
                new_n = ""
                if s[i] == ")":##end of the parenthesis
                    i+=1
                    break
            else:
                new_n += s[i]

        i+=1
        
    ##reading the internal node info

    while i < len(s):
        if s[i] == beginCommentChar:
            inComment = True
        elif s[i] == endCommentChar:
            inComment = False

        #print i , s[i] , inComment

        if not inComment:
            if s[i] in [")",";",","]:
                break
        info += s[i]
        i+=1



    #### building tree

    annot = parsePrIMEAnnotations( info )

    #print annot , [c.name for c in childrenNodes]

    evt = PrIMEAnnotationToRecEvent(annot)

    IntNode = ReconciledTree()
    IntNode.setName(annot["name"])

    if annot.has_key("BL"):
        IntNode.dist = annot["BL"]

    IntNode.addEvent(evt)

    ## checking if we are transfer
    if isTansferNode(IntNode):
        for c in childrenNodes:
            if not c.sameSpeciesAsParent(parent = IntNode): ## we found a children node
                addTbEvent(c) ## we add the transferBack event



    if len(childrenNodes) > 1:
        ## two or more children, add them
        for c in childrenNodes:
            IntNode.add_child(c)
    elif len(childrenNodes) == 1:
        ## only one children. add its events to the current node
        for e in childrenNodes[0].getEvents():
            #print e 
            IntNode.addEvent(e)

        for c in childrenNodes[0].get_children():
            IntNode.add_child(c)


    #print annot["name"],"-->",[n.getTreeNewick() for n in childrenNodes]
    
    return i, IntNode


def pruneReconciledTree(Rnode):
    """
    *RECURSIVE*
    For the special cases of unpruned trees.
    We want to avoid cases of duplicationLoss (absent from the current format)
    We also want to avoid cases of SpeciationOutLoss with only one child is the same species
    In these two cases, the corresponding RecEvent is simply popped

    NB : both these events cannot (in theory) be the last event of the chain
    """

    toPop = []

    evts = Rnode.getEvents()
    for i , e in enumerate(evts):
        if e.eventCode == "DL":
            toPop.append(i)
        elif e.eventCode == "broL":
            if e.species != evts[i+1].species: ## broL with the transfered child deleted
                toPop.append(i)

    for i in toPop[::-1]:
        Rnode.popEvent(i)

    for c in Rnode.get_children():
        pruneReconciledTree(c)

    return 


### tesing lines ###

primeLine = "((G0_15:1.4895169[&&PRIME ID=0 VERTEXTYPE=Leaf],G1_16:1.4895169[&&PRIME ID=1 VERTEXTYPE=Leaf])G2_17:1.3832329[&&PRIME ID=2 VERTEXTYPE=Speciation DISCPT=(1,0)],((G3_2:1.5396005[&&PRIME ID=3 VERTEXTYPE=Leaf],(G4_9:0.34305935[&&PRIME ID=4 VERTEXTYPE=Leaf],G5_15:0.34305935[&&PRIME ID=5 VERTEXTYPE=Leaf])G6_9:1.1965411045745056[&&PRIME ID=6 VERTEXTYPE=Transfer FROMTOLINEAGE=(9,15) DISCPT=(0,1)])G7_9:0.061006042747293265[&&PRIME ID=7 VERTEXTYPE=Transfer FROMTOLINEAGE=(9,2) DISCPT=(1,1)],(((G8_0:0.27117633[&&PRIME ID=8 VERTEXTYPE=Leaf],G9_12:0.27117633[&&PRIME ID=9 VERTEXTYPE=Leaf])G10_0:0.9953092605365245[&&PRIME ID=10 VERTEXTYPE=Transfer FROMTOLINEAGE=(0,12) DISCPT=(0,1)],G11_15:1.2664856[&&PRIME ID=11 VERTEXTYPE=Leaf])G12_15:0.22303131363405873[&&PRIME ID=12 VERTEXTYPE=Transfer FROMTOLINEAGE=(15,0) DISCPT=(0,2)],(((G13_15:0.4357554[&&PRIME ID=13 VERTEXTYPE=Leaf],G14_11:0.4357554[&&PRIME ID=14 VERTEXTYPE=Leaf])G15_15:0.4148762069265683[&&PRIME ID=15 VERTEXTYPE=Transfer FROMTOLINEAGE=(15,11) DISCPT=(0,1)],G16_16:0.85063161[&&PRIME ID=16 VERTEXTYPE=Leaf])G17_16:0.46392360870188326[&&PRIME ID=17 VERTEXTYPE=Transfer FROMTOLINEAGE=(16,15) DISCPT=(0,2)],(G18_5:0.43323511[&&PRIME ID=18 VERTEXTYPE=Leaf],(G19_12:0.037397094[&&PRIME ID=19 VERTEXTYPE=Leaf],G20_16:0.037397094[&&PRIME ID=20 VERTEXTYPE=Leaf])G21_12:0.39583801225313864[&&PRIME ID=21 VERTEXTYPE=Transfer FROMTOLINEAGE=(12,16) DISCPT=(0,1)])G22_5:0.8813201085697117[&&PRIME ID=22 VERTEXTYPE=Transfer FROMTOLINEAGE=(5,12) DISCPT=(0,1)])G23_16:0.17496168551161417[&&PRIME ID=23 VERTEXTYPE=Transfer FROMTOLINEAGE=(16,5) DISCPT=(0,2)])G24_17:0.1110896[&&PRIME ID=24 VERTEXTYPE=Speciation DISCPT=(1,0)])G25_17:1.272143277023072[&&PRIME ID=25 VERTEXTYPE=Transfer FROMTOLINEAGE=(17,9) DISCPT=(1,1)])G26_17:4.9420893[&&PRIME ID=26 VERTEXTYPE=Duplication DISCPT=(4,2)][&&PRIME NAME=PrunedTree];"

primeLine2 = "(G5_0:0.023926212087861037[&&PRIME ID=5 VERTEXTYPE=Loss],((G0_1:0.050143376[&&PRIME ID=0 VERTEXTYPE=Leaf],G1_2:0.050143376[&&PRIME ID=1 VERTEXTYPE=Leaf])G2_3:0.32098622[&&PRIME ID=2 VERTEXTYPE=Speciation DISCPT=(1,0)],G3_4:0.3711296[&&PRIME ID=3 VERTEXTYPE=Leaf])G4_5:0.78187421[&&PRIME ID=4 VERTEXTYPE=Speciation DISCPT=(2,0)])G6_6:0.0[&&PRIME ID=6 VERTEXTYPE=Speciation DISCPT=(3,0)][&&PRIME NAME=UnprunedTree];"

#i,T = readPrIMEParenthesis(primeLine)
#print T.getTreeRecPhyloXML()
#print "***"
#i,T = readPrIMEParenthesis(primeLine2)
#print T.getTreeRecPhyloXML()


import sys

if __name__ == "__main__":

    help =  """
                Given a file containing a reconciled tree in PrIME format, 
                this script writes the tree in recPhyloXML format.

                usage : python PrIMEtoRecPhyloXML.py -g geneFileIn [-o fileOut ]
                            -g geneFileIn       : name of the file containing PrIME reconciliations
                            -o fileOut          : (optional) name of the output file (default is geneFileIn + ".xml" )
               """


    OK = True

    nextKEY = None
    params = {
                            "-g"    : None ,#name of the file containing NHX reconciliations                        
                            "-o"    : None ,#(optional) name of the output file (default is geneFileIn + ".xml" )

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
        else:
            print "unknown argument", sys.argv[i]


    if params["-g"] is None:
        OK = False
        print "error: gene input file not given."



    if not OK:
        print help
        exit(1)



    defaultOutputSuffix = ".xml"
    if params["-o"] is None:
        params["-o"] = params["-g"] + defaultOutputSuffix


    #print "reading input gene tree."

    IN = open(params["-g"],"r")
    l = IN.readline() ## this very crude converter doesn't care about the other lines
    IN.close()

    l = l.strip()

    i,T = readPrIMEParenthesis(l) ## this function read the notung like and build a ReconciledTree object from it
    pruneReconciledTree(T)

    OUT = open(params["-o"],"w")
    OUT.write(T.getTreeRecPhyloXML()) ## getting the recPhyloXML of a Reconciled tree is easy
    OUT.close()

    print "reconciled tree converted and written."