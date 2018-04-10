#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        22-Mar-2017         
##  Last modified:  10-Apr-2018        
##
##  Decribes functions to transform a reconciled tree in 
##  NOTUNG format into a tree in the recPhyloXML format
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             ReconciledTree , and ete3-based representation of a reconciled tree
## 
##  developped for python2.7
##
#########################################


import ete3 
from ReconciledTree import RecEvent, ReconciledTree

def myBasicTreeXMLLinesAux(tree):
    """
    Takes:
        - tree (ete3.TreeNode)

    Returns:
        (list): list of xml lines
    """

    indentChar = "  "

    lines = ["<clade>"]

    lines.append( indentChar + "<name>" + tree.name + "</name>" )

    for c in tree.children:
        tmp = myBasicTreeXMLLinesAux(c)
        for l in tmp:
            lines.append( indentChar + l )

    lines.append("</clade>")

    return lines

def myBasicTreeXMLLines(tree):
    """
    Takes:
        - tree (ete3.TreeNode)

    Returns:
        (list): list of xml lines
    """
    lines = ["<phylogeny>"]
    indentChar = "  "
    tmp = myBasicTreeXMLLinesAux(tree)
    for l in tmp:
            lines.append( indentChar + l )

    lines.append("</phylogeny>")

    return lines


def AddLeafEvent(node):
    """ adds a Leaf event to the node """
    sp = node.getEvents()[-1].species
    node.addEvent( RecEvent("C",sp) )
    node.getEvents()[-1].additionnalInfo["geneName"] = node.name

def setNodeAsLeaf(node):
    """ transform a X into a Leaf event """
    node.getEvents()[0].eventCode = "C"
    node.getEvents()[0].additionnalInfo["geneName"] = node.name


def setNodeAsTransfer(node):
    """ transform a S into a bro event """
    node.getEvents()[0].eventCode = "bro"

def setNodeAsXloss(node):
    """ transform a S into a SL event, a bro into a broL, ... """
    node.getEvents()[0].eventCode += "L"


def isXnode(node,X,i=0):
    """ checks if the event at index i (defautl = 0) of the node is a X """

    return node.getEvents()[i].eventCode == X

def isLossNode(node):
    """ checks if the first event of the node is a LOSS """
    return isXnode(node,"L")

def isSpeciationNode(node):
    """ checks if the first event of the node is a speciation """
    return isXnode(node,"S")

def isLeafNode(node):
    """ checks if the LAST event of the node is a leaf event """
    return isXnode(node,"C",-1)


def isTransferBackNode(node):
    """ checks if the first event of the node is a transferBack """
    return isXnode(node,"Tb")


def NotungAnnotationToRecEvent(annot):
    """ returns a RecEvent OR None in case of problem """

    evtType = None
    sp = annot.get("S",None)

    if annot.has_key("D"): # duplication
        if annot["D"] == "Y":
            evtType = "D"

    if annot.has_key("H"):# transfer reception
        if annot["H"].startswith("Y"):
            evtType = "Tb"
            ## the destination species is already in the S field. The sending species
            #sHinfo = annot["H"].split("@")
            ## typical line : 'Y@MOUSE@GORILLA' for a transfer from mouse to gorilla


    if evtType is None:
        # speciation, branchingOut or Loss or leaf!
        if annot["name"].endswith("*LOST"):
            evtType = "L"
        else:
            evtType = "S" ## between leaf, S and bro, always put S and we'll complement later

    if ( sp is None ) or (evtType is None):
        print "error when trying to assign event. data:", annot
        return None
    return RecEvent(evtType,sp)


def parseNotungAnnotations(info):

    nameBL,j,comment = info.partition("[")

    if comment.endswith("]"):
        comment = comment[:-1]

    name, j, BL = nameBL.partition(":")

    commentList = comment.split(":")
    
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
    annot = parseNotungAnnotations(s)

    evt = NotungAnnotationToRecEvent(annot)

    leafNode = ReconciledTree()
    leafNode.setName(annot["name"])

    if annot.has_key("BL"):
        leafNode.dist = annot["BL"]

    leafNode.addEvent(evt)

    if isSpeciationNode(leafNode):
        setNodeAsLeaf(leafNode)

    if ( not isLeafNode(leafNode) ) and (not isLossNode(leafNode)):
        AddLeafEvent(leafNode)

    return leafNode

def readNotungParenthesis(s):
    
    if not s.startswith("("):
        print "error, should begin with a parenthesis"
        return

    i = 1    
    new_n = ""
    info = "" 

    childrenNodes = []

    while i < len(s):
        if s[i] == "(":##new parenthesis
            offset , child = readNotungParenthesis(s[i:])
            i += offset
            childrenNodes.append(child)
#            print "added new internal."

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

    while s[i] not in [")",";",","] and i < len(s):
        info += s[i]
        i+=1


    annot = parseNotungAnnotations( info )

    evt = NotungAnnotationToRecEvent(annot)

    IntNode = ReconciledTree()
    IntNode.setName(annot["name"])

    if annot.has_key("BL"):
        IntNode.dist = annot["BL"]

    IntNode.addEvent(evt)

    ## checking if : one of the child is a loss and/or a child is a transfer reception
    hasLoss = False
    hasTransfer = False
    j = 0
    while j < len(childrenNodes):
        c  = childrenNodes[j]
        if isLossNode(c):
            hasLoss = True
            #childrenNodes.pop(j)
            #continue
        if isTransferBackNode(c):
            hasTransfer = True
        j +=1


    if hasTransfer:
        setNodeAsTransfer(IntNode)

    #if hasLoss:
    #    setNodeAsXloss(IntNode)

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



    #print IntNode.getTreeRecPhyloXML()

    #print annot["name"],"-->",[n.getTreeNewick() for n in childrenNodes]
    
    return i, IntNode


## test chain for testing
TESTnotung = '((((gB_human[&&NHX:S=HUMAN:Nset=<HUMAN>],gA_human[&&NHX:S=HUMAN:Nset=<HUMAN>])n2[&&NHX:S=HUMAN:Nset=<HUMAN>:D=Y],GORILLA*LOST[&&NHX:S=GORILLA:Nset=<GORILLA>])r21[&&NHX:S=n28:Nset=<HUMAN@GORILLA>],((gA_mouse[&&NHX:S=MOUSE:Nset=<MOUSE>],g_gorilla[&&NHX:S=GORILLA:Nset=<GORILLA>:H=Y@MOUSE@GORILLA])n5:73.0[&&NHX:S=MOUSE:Nset=<MOUSE>:B=73.0],gB_mouse[&&NHX:S=MOUSE:Nset=<MOUSE>])n7[&&NHX:S=MOUSE:Nset=<MOUSE>:D=Y])n8:100.0[&&NHX:S=n30:Nset=<n28@MOUSE>:B=100.0],(gY_cow[&&NHX:S=COW:Nset=<COW>],gX_cow[&&NHX:S=COW:Nset=<COW>])n11:100.0[&&NHX:S=COW:Nset=<COW>:D=Y:B=100.0])n12[&&NHX:S=n32:Nset=<COW@n30>];'



import sys

if __name__ == "__main__":

    help =  """
                Given a file containing a reconciled tree in NOTUNG format, 
                this script writes the tree in recPhyloXML format.

                usage : python NOTUNGtoRecPhyloXML.py -g geneFileIn [-o fileOut ]
                            -g geneFileIn       : name of the file containing NOTUNG reconciliations
                            -o fileOut          : (optional) name of the output file (default is geneFileIn + ".xml" )

                            --include.species   : (optional) whether the species tree should be included in the XML file (using the <spTree> tag)
               """


    OK = True

    nextKEY = None
    params = {
                            "-g"    : None ,#name of the file containing NHX reconciliations                        
                            "-o"    : None ,#(optional) name of the output file (default is geneFileIn + ".xml" )
                            "--include.species" : False #(optional) whether the species tree should be included in the XML file (using the <spTree> tag)
            }

    flagArgs = ["--include.species"]

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

    geneLine = IN.readline()

    speciesLine = IN.readline()

    IN.close()

    geneLine = geneLine.strip()

    i,RT = readNotungParenthesis(geneLine) ## this function read the notung like and build a ReconciledTree object from it

    speciesLine = "(" + speciesLine.strip().partition("(")[2][:-1] + ";" ## [&&NOTUNG-SPECIES-TREE(((GORILLA,HUMAN)n28,MOUSE)n30,COW)n32] to the newick alone
    speciesTree = ete3.Tree(speciesLine , format = 1)


    OUT = open(params["-o"],"w")

    indentLevel = 0
    indentChar = "  "

    if params["--include.species"]:

        OUT.write( "<recPhylo>" + "\n" )
    
        indentLevel += 1

        OUT.write( indentLevel * indentChar + "<spTree>" + "\n")

        indentLevel += 1
        
        lines = myBasicTreeXMLLines(speciesTree)

        for xmlline in lines:
                OUT.write( indentLevel * indentChar + xmlline + "\n" )             

        indentLevel -= 1
        OUT.write( indentLevel * indentChar + "</spTree>" + "\n")


    XMLlines = RT.getTreeRecPhyloXMLLines()
    for xmlline in XMLlines:
        OUT.write( indentLevel * indentChar + xmlline + "\n" ) 





    if params["--include.species"]:
        OUT.write( "</recPhylo>" + "\n" )

    OUT.close()

    print "reconciled tree converted and written."