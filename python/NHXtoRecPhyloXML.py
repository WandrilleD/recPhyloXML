#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        24-Feb-2017         
##  Last modified:  20-Sept-2017        
##
##  Decribes functions to transform a reconciled tree in 
##  NHX format into a tree in the recPhyloXML format
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             ReconciledTree , and ete3-based representation of a reconciled tree
## 
##  developped for python2.7
##
#########################################


from ete3 import Tree, TreeNode
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

def NHXtreeToBasicRecTree(nhxTree , spTree = None):
    """ 
    *RECURSIVE*

    From a tree read in a NHX format file to a ReconciledTree without any intermediary events (SpeciationLoss events)

    Takes:
        - nhxTree (ete3.TreeNode) : tree read from a NHX formatted line
        - spTree (ete3.Tree or None) [default = None] : if different from None, 
                                            internal node's events associated to species whose name is in the species tree will be kept as such
                                                        if equal to None,
                                            only leaves get to keep their associated species (and species of other events will have to be re-associated later)

    Returns:
        (ReconciledTree)

    """

    RT = ReconciledTree()

    eventCode = None
    species = None

    ## only terminal events in a DL context are considered here : leaf, speciation or duplication

    if nhxTree.is_leaf():
        eventCode = "C"
        species = nhxTree.S 
        ## we only get the species for the leaves 
        ##( as they are the only one where we   are sure the species is one that is present in the species tree)

    elif nhxTree.D == "Y":
        eventCode = "D"
    else:
        eventCode = "S"

    if not spTree is None:
        if len(spTree.search_nodes(name = nhxTree.S)) == 1: ## not at all efficient ...
            species = nhxTree.S


    ##additional info:
    for f in nhxTree.features:
        RT.add_feature( f , nhxTree.__getattribute__(f) ) ## stinks that we call thi __method__ ...

    evt = RecEvent(eventCode , species)
    RT.addEvent(evt)

    for c in nhxTree.children:
        RT.add_child( NHXtreeToBasicRecTree(c) )
    return RT


def completeTreeNames(tree):
    """
    Takes:
        - tree (ete3.Tree)

    Returns:
        (ete3.Tree) : the tree, but where the nodes without a name now have one that correspond
                    to their post-order

    """

    for i,n in enumerate(tree.traverse('postorder')):
        if n.name == "":
            n.name = str(i)
    return tree

def annotateIncompleteRecRTree(recTree , spTree):
    """
    Takes:
        - recTree (ReconciledTree)
        - spTree (ete3.Tree)

    Returns:

    """
    ### post order traverse
    ## -> setup species using LCA of children
    ## -> add SL using dist from parent to CH

    ## for quicker access
    dspNameToNode = {n.name : n for n in spTree.traverse()}

    NODETODO = [n for n in recTree.traverse("postorder")]
    for n in NODETODO:

        if n.is_leaf():
            continue

        childrenSpecies = []
        for c in n.children:
            spNode = dspNameToNode.get( c.getEvents()[-1].species , None)

            if spNode is None: ##overkill
                print "ERROR : species",c.getEvents()[-1].species,"unknown..."
                exit(1)
            childrenSpecies.append( spNode )

        if n.getEvents()[-1].species is None:
            ## we have to assign a species
            ## we assign the LCA (this is a little 'by default' and presume parsimony)

            n.getEvents()[-1].species = spTree.get_common_ancestor( childrenSpecies ).name


        currentSpNode = dspNameToNode.get( n.getEvents()[-1].species , None)

        if currentSpNode is None:
            print "ERROR : species",n.getEvents()[-1].species,"unknown..."
            exit(1)


        ChSpMustBeEqual = ( n.getEvents()[-1].eventCode != "S" ) ##unless this is a speciation, the next event shall have the same species


        ## now we want to add SL events
        for i,c in enumerate(n.children):

            childrenSpNode = childrenSpecies[i]
            
            if childrenSpNode == currentSpNode:
                continue ## already ok


            while childrenSpNode.up != currentSpNode:
                #print "adding SL in ", childrenSpNode.up.name , "from" , childrenSpNode.name, 'loss in', getSister(childrenSpNode).name
                #print "in node ", c.name, "father is ", n.name
                addSpeciationAndLoss(c , childrenSpNode)
                #e = RecEvent( "SL" , childrenSpNode.up.name )
                #c.addEvent( e,append = False) ##will insert the evt in first position

                childrenSpNode = childrenSpNode.up

            ## last SL in case the parent is not a speciation
            if ChSpMustBeEqual:
                #print "adding SL in ", childrenSpNode.up.name , "from" , childrenSpNode.name, 'loss in', getSister(childrenSpNode).name
                #print "in node ", c.name, "father is ", n.name
                addSpeciationAndLoss(c , childrenSpNode)
                #e = RecEvent( "SL" , childrenSpNode.name )
                #c.addEvent( e,append = False) ##will insert the evt in first position                

    return recTree


def addSpeciationAndLoss(node , keptSpeciesNode):
    """
    *modifies node in place*

    Takes:
        - node (ReconciledTree): node where a SpeciationLoss must take place
        - keptSpeciesNode (ete3.TreeNode): node of the species tree where the lineage survived (ie. the sister species of the one where the loss occured)
    """
    parentSpeciesNode = keptSpeciesNode.up
    lossSpeciesNode = getSister(keptSpeciesNode)

    lossNode = ReconciledTree()
    lossNode.addEvent( RecEvent("loss" , lossSpeciesNode.name) )
    lossNode.name="LOSS"

    # 2. create the kept child

    keptNode = ReconciledTree()

    ##transfering the events of node to keptNode
    while len(node.eventRecs) >  0:
        keptNode.addEvent( node.popEvent(0), append=False )

    # 3. link children to kept child
    while len(node.children) > 0:
        c = node.children[0]
        c.detach()
        keptNode.add_child(c)


    # 4. branching loss and kept to original node

    node.add_child(lossNode)
    node.add_child(keptNode)

    # 5. editing the event
    e = RecEvent( "S" , keptSpeciesNode.up.name )
    node.addEvent( e,append = False) ##will insert the evt in first position
    
    return





def getSister(node):
    """ presumes a binary tree """
    if node.is_root():
        return None

    C = node.up.children
    for c in C:
        if c != node:
            return c
    raise Exception("tree format error: expected a binary tree")

import sys

if __name__ == "__main__":

    help =  """
                Given a file containing reconciled trees in NHX format (and containing no transfers), 
                and their species tree in newick format
                this script writes the trees in recPhyloXML format.

                usage : python NHXtoRecPhyloXML.py -g geneFileIn -s speciesFileIn [-o fileOut --include.species]
                            -g geneFileIn       : name of the file containing NHX reconciliations
                            -s speciesFileIn    : name of the species tree file
                            -o fileOut          : (optional) name of the output file (default is geneFileIn + ".xml" )
                            --include.species   : (optional) whether the species tree should be included in the XML file (using the <spTree> tag)
               """


    OK = True

    nextKEY = None
    params = {
                            "-g"    : None ,#name of the file containing NHX reconciliations
                            "-s"    : None ,#name of the species tree file
                            "-o"    : None ,#(optional) name of the output file (default is geneFileIn + ".xml" )
                            "--include.species"   : False #(optional) whether the species tree should be included in the XML file (using the <spTree> tag)
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
    if params["-s"] is None:
        OK = False
        print "error: species input file not given."

    if OK:
        
        ## treating positive float options
        for pname in []:
            try:
                params[pname] = float(params[pname])
                if params[pname] < 0:
                    print "error: ",pname ,"must be a positive number."
                    OK = False
            except:
                print "error:",pname,"must be a positive number."
                OK = False

        ## treating positive int options
        for pname in []:
            try:
                params[pname] = int(params[pname])
                if params[pname] < 1:
                    print "error: ",pname ,"must be a positive integer."
                    OK = False
            except:
                print "error:",pname,"must be a positive number."
                OK = False

    if not OK:
        print help
        exit(1)



    defaultOutputSuffix = ".xml"
    if params["-o"] is None:
        params["-o"] = params["-g"] + defaultOutputSuffix


    print "reading input species tree."

    spTree = Tree( params["-s"] )

    spTree = completeTreeNames( spTree )

    
    OUT = open(params["-o"],"w")

    OUT.write( "<recPhylo>" + "\n" )

    indentLevel = 1
    indentChar = "  "

    if params["--include.species"]:

        OUT.write( indentLevel * indentChar + "<spTree>" + "\n")

        indentLevel += 1
        lines = myBasicTreeXMLLines(spTree)
        for xmlline in lines:
                OUT.write( indentLevel * indentChar + xmlline + "\n" )             

        indentLevel -= 1
        OUT.write( indentLevel * indentChar + "</spTree>" + "\n")

    print "reading input reconciled trees."

    IN = open(params["-g"],"r")

    l = IN.readline()

    while l != "":

        if l != "\n":##special ignore white lines

            NHXTree = Tree(l)

            RT = NHXtreeToBasicRecTree(NHXTree , spTree)

            RT = annotateIncompleteRecRTree(RT,spTree)

            RT = completeTreeNames(RT)

            XMLlines = RT.getTreeRecPhyloXMLLines()

            for xmlline in XMLlines:
                OUT.write( indentLevel * indentChar + xmlline + "\n" ) 


        l = IN.readline()

    IN.close()


    OUT.write( "</recPhylo>" + "\n" )

    OUT.close()