#!/usr/bin/python

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        20-June-2016        
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

def treatLossEvent(node, LossEvent , keptChildNameSuffix = ".c"):
    """
    Takes:
        - node (ReconciledTree)
        - LossEvent (recEvent)
        - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss 
    """


    # 1. create the loss child

    #determining ts if there is one
    lostTS = LossEvent.timeSlice
    if not lostTS is None:
        lostTS -= 1 # one TS more recent than the parent

    lossNode = ReconciledTree()
    lossNode.addEvent( RecEvent("loss" , "", ts= lostTS ) )
    lossNode.name="LOSS"

    # 2. create the kept child

    keptNode = ReconciledTree()
    keptNode.name = node.name + keptChildNameSuffix

    # 4. branching loss and kept to original node

    node.add_child(lossNode)
    node.add_child(keptNode)

    # 5. editing the event and adding to node

    e = LossEvent.eventCode
    LossEvent.eventCode = e.rpartition("L")[0]

    node.addEvent(LossEvent)

    return keptNode


def parse_node_annotation(node_annotation, isLeaf = False, isDead = False):
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
        tmp_ann = node_annotation.split(".")
    
        ##further splitting multiple transfer
        s_ann = []
        for ann in tmp_ann:
            if ann.count("@") < 1:
                s_ann.append(ann)
                continue
            ## in case of transfer and loss like: T@27|SYNP2@26|ACAM1
            new_anns = ann.split("@")
            
            s_ann.append("@".join(new_anns[0:2]))##first tranfer, a transfer out
            
            for a in new_anns[2:]:##for each transfer after that (should be only one)
                s_ann.append("@" + a)
            
        
        
        for ann in s_ann:
            
            if ann[0].isdigit(): ##starts with a number spe,dup or spe+loss
                if ann.isdigit(): ##only numbers: spe or spe+loss
                    target = ann
                    
                    l_events.append( RecEvent( "S" , target , int(target) ) )
                    continue
            
            if ann.startswith("T@"): ##Transfer out
                
                ## time slice of the source
                source_ts,junk,source_sp = ann[2:].partition("|")## partitionning something like T@22|22
                
                ##adding the event
                l_events.append( RecEvent( "So" , source_sp , int(source_ts) ) )
                                
            
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
                ts,junk,sp = ann[2:].partition("|")## partitionning something like D@23|22
                
                l_events.append( RecEvent( "D" , sp, int(ts) ) )

            
    if isLeaf and ( len(l_events)==0 or l_events[-1].eventCode !="C" ) :
        l_events.append( RecEvent("C",None,0) ) ##temporary placeholder for the leaf species

    if isDead: ## we start in the dead so the first event MUST be Bout or Tloss

        if not l_events[0].eventCode in ["Tb", "Bo"]:

            target_ts = l_events[0].timeSlice
            target_sp = l_events[0].species

            l_events.insert(0, RecEvent( "Tb" , target_sp , target_ts ) )


    ##adding loss labels
    for i in range(len(l_events)-1):##all events but the last one
        if l_events[i].eventCode in ["So","S"]:
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


    x = 0
    SpFound = False
    AnnotFound = False

    while (not AnnotFound) and (x < len(leafName)):

        if not SpFound:
            if leafName[x] == sepSp:
                SpFound=True
        elif leafName[x] in sepAnnot:
            AnnotFound=True
            break

        x+=1

    return leafName[:x] , leafName[x:]


def getLeafSpeciesFromLeafName(leafName, sepSp = "_"):
    """
    Takes:
         - leafName (str) : name of a leaf, in the format: species separator gene
         - sepSp (str) [default = "_" ] : separator between species name and gene name
    
    Returns:
        (str) : species name
    """
    return leafName.partition(sepSp)[0]


def ALEtreeToReconciledTree(ALEtree, isDead = False):
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

    events = parse_node_annotation(annotation, isLeaf, isDead = isDead)



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
        if e.eventCode.endswith("L"):
            print "plep"
            current = treatLossEvent(current, e , ".c")
        else:
            current.addEvent(e)

    for c in ALEtree.children: ##recursion on successors
        current.add_child( ALEtreeToReconciledTree(c, isDead ) )

    return RT






TESTING_SPECIES_TREE_NEWICK = "(((a,b)2,(c,d)1)4,((e,f)3,g)5)6;"
TESTING_GENE_TREE_NEWICK = "((a_a1,b_b1@0|b).4.2.T@1|a,((((d_d1@0|d,e_e2.3)Tb@3|3,g_g2)T@4|g,g_g3.T@4|3@1|g).5,(e_e1.3,g_g1).5)D@5|5).6;"

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


#    spTree = None
#
#    if not params["-s"] is None:
#
#        print "reading input species tree."
#
#        spTree = Tree( params["-s"] )
#
#        spTree = completeTreeNames( spTree )
#
    
    OUT = open(params["-o"],"w")

    OUT.write( "<recPhylo>" + "\n" )

    indentLevel = 1
    indentChar = "  "

#    if params["--include.species"]:
#
#        if spTree is None:
#            print "activated --include.species without providing a species tree (-s option). Abort."
#            exit(1)
#
#        OUT.write( indentLevel * indentChar + "<spTree>" + "\n")
#
#        indentLevel += 1
#        lines = myBasicTreeXMLLines(spTree)
#        for xmlline in lines:
#                OUT.write( indentLevel * indentChar + xmlline + "\n" )             
#
#        indentLevel -= 1
#        OUT.write( indentLevel * indentChar + "</spTree>" + "\n")


    print "reading input reconciled trees."

    IN = open(params["-g"],"r")

    l = IN.readline()

    while l != "":

        if l != "\n":##special ignore white lines

            ALEtree = Tree( l, format = 1 )

            RT = ALEtreeToReconciledTree(ALEtree)


            XMLlines = RT.getTreeRecPhyloXMLLines()

            for xmlline in XMLlines:
                OUT.write( indentLevel * indentChar + xmlline + "\n" ) 


        l = IN.readline()

    IN.close()


    OUT.write( "</recPhylo>" + "\n" )

    OUT.close()

    print "reconciled tree converted and written."