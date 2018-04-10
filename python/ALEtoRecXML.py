#!/usr/bin/python

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        20-June-2016        
##  Last modified:  10-Apr-2018
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
from ReconciledTree import RecEvent, ReconciledTree, myBasicTreeXMLLines

def completeTreeNames(tree, useBS = False ) :
    """
    Takes:
        - tree (ete3.Tree)
        - useBS (bool) [default = False] : uses bootstrap to name nodes

    Returns:
        (ete3.Tree) : the tree, but where the nodes without a name now have one that correspond
                    to their post-order OR their bootstrap 

    """

    for i,n in enumerate(tree.traverse('postorder')):
        if n.name == "":
            print n.support
            if useBS:
                n.name = str(int(n.support))
            else:
                n.name = str(i)
    return tree



def split_node_annotation(node_annotation, isUndated = False):
    """
    Takes:
        - node_annotation (str): reconciliation annotation on a node
    
    Returns:
        (list): list of of the string representation of the different events
    """


    if len(node_annotation) != 0:
    
        if node_annotation.startswith("."):
            node_annotation = node_annotation[1:]


        s_ann = []

        cur_annot = ""
        arobaseCount = 0
        symbolsToCheck = None
        for i in range(len(node_annotation)):

            if symbolsToCheck is None:
                if node_annotation[i:].startswith("."):

                    s_ann.append(cur_annot)##anotation finished
                    

                    cur_annot = ""
                    arobaseCount = 0

                    continue

                if node_annotation[i] == "@":
                    arobaseCount += 1

                    if arobaseCount > 1:
                        ## in case of transfer and loss like: T@27|SYNP2@26|ACAM1
                        s_ann.append(cur_annot)##anotation finished
                        ### treat this annotation!
                        cur_annot = ""
                        arobaseCount = 0

                    ## looking for a symbol that must happen before there is any new event
                    if isUndated:
                        if cur_annot[-1] != "D": ## if not a duplication
                            symbolsToCheck = "->"
                    else:
                        symbolsToCheck = "|"




            elif node_annotation[i:].startswith(symbolsToCheck):
                symbolsToCheck = None


            cur_annot += node_annotation[i]


        if len(cur_annot)> 0:
            s_ann.append(cur_annot)##anotation finished


    # when dated, events are from oldest to newest
    # when undated they are from newest to oldest
    if isUndated:
        s_ann = s_ann[::-1]

    return s_ann
            
#split_node_annotation("T@27|SYNP2@26|ACAM1",  isUndated = False)
#split_node_annotation(".4.2.T@1|a",  isUndated = False)
#split_node_annotation("@0|d",  isUndated = False)
#split_node_annotation("Tb@3|3" , isUndated = False)
#split_node_annotation(".T@4|3@1|g" , isUndated = False)
#split_node_annotation("D@5|5" , isUndated = False)
#split_node_annotation(".5" , isUndated = False)

#split_node_annotation("T@GorillaSTL->39", isLeaf = False, isDead = False, isUndated = True)
#split_node_annotation("39.T@CRBM.05->39", isLeaf = False, isDead = False, isUndated = True)
#split_node_annotation(".D@D" , isUndated = True)
#split_node_annotation(".T@D->A" , isUndated = True)
#split_node_annotation(".D@D.4" , isUndated = True) ###### ambiguity of the format? 

#split_node_annotation("T@Saki->Chimp1.T@38->Saki" , isUndated = True)
#split_node_annotation("D@43.46.49" , isUndated = True)




def parse_node_annotation(node_annotation, isLeaf = False, isDead = False, isUndated = False):
    """
    Takes:
        - node_annotation (str): reconciliation annotation on a node
    
    Returns:
        (list): list of dicts coding a particular event
    """
    #print "annot : ",node_annotation

    l_events = []

    if len(node_annotation) != 0:
    
        if node_annotation.startswith("."):
            node_annotation = node_annotation[1:]

        s_ann = split_node_annotation(node_annotation, isUndated )

        #print node_annotation, "->", s_ann
        
        for ann in s_ann:
            
            if ann[0].isdigit(): ##starts with a number spe,dup or spe+loss
                if ann.isdigit(): ##only numbers: spe or spe+loss
                    target = ann
                    
                    ts = int(target)

                    if isUndated:
                        ts = None

                    l_events.append( RecEvent( "S" , target , ts ) )
                    continue
            
            if ann.startswith("T@"): ##Transfer out
                
                ## time slice of the source
                source_ts = None
                source_sp = None
                if isUndated:
                    ## of the shape : "T@D->A"
                    source_sp = ann[2:].partition("->")[0]
                else:
                    source_ts,junk,source_sp = ann[2:].partition("|")## partitionning something like T@22|22
                    source_ts = int(source_ts)

                ##adding the event
                l_events.append( RecEvent( "bro" , source_sp , source_ts ) )
                                
            
            if ann.startswith("@"): # or ann.startswith("Tb@"):##transfer in or back
            
                pre = 3##cropping size
                if ann.startswith("@"):
                    pre = 1
    
                target_ts,junk,target_sp = ann[pre:].partition("|")## partitionning something like @22|22 to get the time slice and specie
    
                ##adding the event
                l_events.append( RecEvent( "Tb" , target_sp, int(target_ts) ) )


    
            if ann.startswith("Tb@"):
                l_events.append( RecEvent("Bo", "-1" ) )
                                
            if ann.startswith("D@"):##Duplication

                ts = None
                sp = None
                if isUndated:
                    sp = ann[2:]
                else:    
                    ts,junk,sp = ann[2:].partition("|")## partitionning something like D@23|22
                    ts = int(ts)
                
                l_events.append( RecEvent( "D" , sp, ts ) )

            
    if isLeaf and ( len(l_events)==0 or l_events[-1].eventCode !="C" ) :
        ts = 0
        if isUndated:
            ts = None
        l_events.append( RecEvent("C",None,ts) ) ##temporary placeholder for the leaf species

    if isDead: ## we start in the dead so the first event MUST be Bout or Tloss

        if not l_events[0].eventCode in ["Tb", "Bo"]:

            target_ts = l_events[0].timeSlice
            target_sp = l_events[0].species

            l_events.insert(0, RecEvent( "Tb" , target_sp , target_ts ) )


    ##adding loss labels
    for i in range(len(l_events)-1):##all events but the last one
        if l_events[i].eventCode in ["bro","S"]:
            l_events[i].eventCode += "L"


    return l_events




def separateLeafNameFromLeafAnnotation( leafName , sepSp = "_" , sepAnnot = (".","@") ):
    """
    Takes:
        - leafName (str) : name of a leaf, potentially containing reconciliation information (exemple: "g_g3.T@4|3@1|g" )
        - sepSp (str) [default = "_" ] : separator between species name and gene name
        - sepAnnot (tuple) [default = (".","@") ] : possible separators between gene name and annotations

    Returns:
        (tuple)
            (str) : gene name
            (str) : reconciliation annotation  (empty string if there is none)

    """

    #print leafName,

    x = 0
    AnnotFound = False

    while (not AnnotFound) and (x < len(leafName)):
       
        if leafName[x] in sepAnnot:
            AnnotFound=True
            break

        x+=1
    #print "->", leafName[:x] , leafName[x:]
    return leafName[:x] , leafName[x:]


def getLeafSpeciesFromLeafName(leafName, sepSp = "_"):
    """
    Takes:
         - leafName (str) : name of a leaf, in the format: species separator gene
         - sepSp (str) [default = "_" ] : separator between species name and gene name
    
    Returns:
        (str) : species name
    """
    for i,c in enumerate(leafName):
        if c in [sepSp , "." , "@"]:
            return leafName[:i]

    return leafName.partition(sepSp)[0]


def ALEtreeToReconciledTree(ALEtree, isDead = False, isUndated = False):
    """
    Recursively builds the reconciled tree 

    Takes:
        - ALEtree (ete3.TreeNode) : a tree read from a reconciled ttree in the ALE format (ie. reconciliation annotations are in the .name field)
        - isDead (bool) [default = False] : indicates whether or not this lineage starts in a dead/unsampled species

    Returns:
        (ReconciledTree)
    """

    isLeaf  = ALEtree.is_leaf()

    annotation = None
    name = ALEtree.name
    if isLeaf:
        name , annotation = separateLeafNameFromLeafAnnotation(ALEtree.name)
    else:
        annotation = ALEtree.name

    #print "name : ",ALEtree.name

    events = parse_node_annotation(annotation, isLeaf, isDead = isDead, isUndated = isUndated)



    if isLeaf:
        ## we specify the species of the leaf event
        events[-1].species = getLeafSpeciesFromLeafName( ALEtree.name )

    #print [str(e) for e in events]

    
    if events[-1].eventCode == "Bo":
        isDead = True ## means that both children must begin by an event in the dead


    RT = ReconciledTree()
    RT.setName(name)

    current = RT

    for e in events:
#        if e.eventCode.endswith("L"):
#            #print "plep"
#            current = treatLossEvent(current, e , ".c")
#        else:
#            current.addEvent(e)
#
        current.addEvent(e)

    for c in ALEtree.children: ##recursion on successors
        current.add_child( ALEtreeToReconciledTree(c, isDead , isUndated) )

    return RT


def refineReconciledTreeWithTransferBack(RT):
    """
    adds transferBack events where they were omitted

    Takes:
        - RT (ReconciledTree) : a reconciled tree obtained from an ale string

    """

    for n in RT.traverse():

        if n.is_root():
            continue

        lastEventParent = n.up.getEvent(-1)

        print [x.eventCode+","+x.species for x in n.eventRecs]

        if lastEventParent.eventCode in ["branchingOut", "bro","branchingOutLoss", "broL"]:
            ## parent is an outgoing transfer

            firstEvent = n.getEvent(0)

            print lastEventParent.eventCode , lastEventParent.species , "->",
            print firstEvent.eventCode , firstEvent.species

            if firstEvent.species == lastEventParent.species:
                continue ## this is the "kept" child --> continue

            if firstEvent.eventCode in ["Tb", "transferBack","Bo","bifurcationOut"]:
                continue ## already the correct annotation

            print "  creating event"

            TbEvent = RecEvent("Tb" ,  species = firstEvent.species)

            n.addEvent( TbEvent , append = False)

def refineReconciledTreeLosses(RT, spTree):
    """ 
    adds species to the losses 

    Takes:
        - RT (ReconciledTree) : a reconciled tree obtained from an ale string
        - spTree (ete3.Tree) : a species tree

    """

    for n in RT.traverse():
        if n.is_root():
            continue

        firstEvent = n.getEvent(0)

        if first.eventCode in ["L", "loss"]:
            ## loss!

            lastEventParent = n.up.getEvent(-1)

            if firstEvent.species == lastEventParent.species:
                continue ## this is the "kept" child --> continue

            if firstEvent.eventCode in ["Tb", "transferBack","Bo","bifurcationOut"]:
                continue ## already the correct annotation

            TbEvent = RecEvent("Tb" ,  species = firstEvent.species)

            n.addEvent( TbEvent , append = False)


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



TESTING_SPECIES_TREE_NEWICK = "(((a,b)2,(c,d)1)4,((e,f)3,g)5)6;"
TESTING_GENE_TREE_NEWICK = "((a_a1,b_b1@0|b).4.2.T@1|a,((((d_d1@0|d,e_e2.3)Tb@3|3,g_g2)T@4|g,g_g3.T@4|3@1|g).5,(e_e1.3,g_g1).5)D@5|5).6;"


#undated ex with a duplication : ((A:0.1,B:0.2).5:0.3,(C:0.3,(E:0.4,(D_1:0.5,D_2:0.5).D@D:0.5).6:0.05).7:0.3).8:0;
#undated transfer ((B:0.2,A_1:0.1).5:0.3,(C:0.3,(E:0.4,(D:0.5,A_2:0.5).T@D->A:0.5).6:0.05).7:0.3).8:0;

if __name__ == "__main__":

    import sys


    ##############################################
    help =  """
                Given a file containing reconciled trees in ALE reconciled tree format, 
                this script writes the trees in recPhyloXML format.

                usage : python ALEtoRecPhyloXML.py -g geneFileIn [-o fileOut]
                            -g geneFileIn       : name of the file containing NHX reconciliations
                            -o fileOut          : (optional) name of the output file (default is geneFileIn + ".xml" )

               """
#                            (TODO:)
#                usage : python ALEtoRecPhyloXML.py -g geneFileIn -s speciesFileIn [-o fileOut --include.species]
#                            (-s speciesFileIn    : (optional) name of the species tree file
#                            (--include.species   : (optional) whether the species tree should be included in the XML file (using the <spTree> tag)


    OK = True

    nextKEY = None
    params = {
                            "-g"    : None ,#name of the file containing NHX reconciliations
                            "-o"    : None #(optional) name of the output file (default is geneFileIn + ".xml" )
#                            "-s"    : None ,#name of the species tree file
#                            "--include.species"   : False #(optional) whether the species tree should be included in the XML file (using the <spTree> tag)
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



    
    OUT = open(params["-o"],"w")

    OUT.write( "<recPhylo>" + "\n" )

    indentLevel = 1
    indentChar = "  "



    print "reading input reconciled trees."

    spTree = None
    isUndated = False

    IN = open(params["-g"],"r")

    l = IN.readline()

    while l != "":

        if l != "\n":
            
            if l.startswith("("):##special ignore white lines

                ALEtree = Tree( l, format = 1 )
    

                RT = ALEtreeToReconciledTree(ALEtree, isUndated = isUndated)
    
                ConvertRTtoLossIndepVersion(RT , speciesTree = spTree, keptChildNameSuffix = ".c")

                if isUndated:
                    refineReconciledTreeWithTransferBack(RT)

    
    
                XMLlines = RT.getTreeRecPhyloXMLLines()
    
                for xmlline in XMLlines:
                    OUT.write( indentLevel * indentChar + xmlline + "\n" ) 
                #break


            elif l.startswith("S:"):
                ## found a species tree!
                treeLine = l[2:].strip()
                print treeLine

                spTree = Tree( treeLine , format = 1 )
        
                spTree = completeTreeNames( spTree , True)
        
                OUT.write( indentLevel * indentChar + "<spTree>" + "\n")
        
                indentLevel += 1
                lines = myBasicTreeXMLLines(spTree)
                for xmlline in lines:
                        OUT.write( indentLevel * indentChar + xmlline + "\n" )             
        
                indentLevel -= 1
                OUT.write( indentLevel * indentChar + "</spTree>" + "\n")


            elif l.startswith("#ALE"): ## trying to recognise "#ALE****_undated " prefix
                if l.partition("_")[2].startswith("undated "):
                    isUndated = True

        l = IN.readline()

    IN.close()


    OUT.write( "</recPhylo>" + "\n" )

    OUT.close()

    print "reconciled tree converted and written."