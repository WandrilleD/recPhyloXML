#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        08-Sept-2017         
##  Last modified:  08-Sept-2017        
##
##  
##  This script is used to combine different recPhyloXML 
##    (containing different reconciled gene trees) files into one.
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             
## 
##  developped for python2.7
##
#########################################



import ete3

import sys
import os


RECPHYLOTAG = "recPhylo"
RECTREETAG = "recGeneTree"
SPTREETAG = "spTree"


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

def readSpeciesTree( filename , isNewick ):
    """

    Takes:
        - filename (str)
        - isNewick (str)
    Returns:
        (tuple)
            (int) : 0 -> no problem
                    1 -> file does not exists
                    2 -> file is not in recPhyloXML (no <recGeneTree> or <recPhylo> tag)

            (list) : species tree lines

    """


    if not os.path.isfile(filename):
        return 1, []

    spLines = ["<" + SPTREETAG + ">\n"]


    returnCode = 0

    if isNewick:
        t = ete3.Tree(filename, format = 1)

        XMLlines = myBasicTreeXMLLines(t)

        for l in XMLlines:
            spLines.append( "  " + l + "\n" ) 
        
    else:
        TAG = "phylogeny"
        foundTAG = False


        offset = 0

        IN = open(filename, "r")
        l = IN.readline()


        while l !="":

            stripped = l.strip()
            if foundTAG:

                spLines.append( "  " + l[offset:] )
                
                if stripped.startswith("</" + TAG) :
                    break
            
            elif stripped.startswith("<" + TAG) :
                    foundTAG = True

                    offset = l.index("<")
                    spLines.append( "  " + l[offset:] )

            l = IN.readline()

        if not foundTAG:
            returnCode = 2

        IN.close()

    spLines.append( "</" + SPTREETAG + ">\n" )

    return returnCode, spLines


def readOneTreeFile( filename ):
    """
    Takes:
        filename (str): name of a file containing a recPhyloXML tree

    Returns:
        (tuple)
            (int) : 0 -> no problem
                    1 -> file does not exists
                    2 -> file is not in recPhyloXML (no <recGeneTree> or <recPhylo> tag)

            (list) : species tree lines
            (list) : recGeneTree lines (potentially several)
    """



    if not os.path.isfile(filename):
        return 1, [],[]

    InvalidTree = True
    
    spLines = []
    recLines = []


    IN = open(filename, "r")
    l = IN.readline()


    readingSpecies = False
    readingRec = False

    offset = 0

    while l !="":

        stripped = l.strip()

        if InvalidTree:
            if stripped.startswith("<" + RECPHYLOTAG) or stripped.startswith("<" + RECTREETAG):
                InvalidTree = False
        
        if readingSpecies:

            ##adding this line to the species lines

            spLines.append( l[offset:] ) ## offsetting to be sure have a nice indentation in the output file

            if stripped.startswith("</" + SPTREETAG): ## end of the species tree
                readingSpecies = False

        elif readingRec:

            ##adding this line to the rec lines

            recLines.append( l[offset:] ) ## offsetting to be sure have a nice indentation in the output file

            if stripped.startswith("</" + RECTREETAG): ## end of the species tree
                readingRec = False

        else:
            if stripped.startswith("<"  + SPTREETAG):

                readingSpecies = True
                offset = l.index("<")
                spLines.append( l[offset:] ) ## offsetting to be sure have a nice indentation in the output file

            elif stripped.startswith("<"  + RECTREETAG):

                readingRec = True
                offset = l.index("<")
                recLines.append( l[offset:] ) ## offsetting to be sure have a nice indentation in the output file

        l = IN.readline()

    IN.close()

    returnCode = 0
    if InvalidTree:
        returnCode = 2

    return returnCode, spLines , recLines








if __name__ == "__main__":

    help =  """
                This script is used to combine different recPhyloXML (containing different reconciled gene trees) files into one.

                NB: if a species tree is present in at least one of the recPhyloXML file, it will be included in the output file (if several files contain a species tree, then only the species tree of the last file will remain).
                    If the -s option is used, then the species tree that is included will be the one in the provided file.


                usage : python combineRecPhyloXMLfiles.py -o fileOut recPhyloXMLfile1 recPhyloXMLfile2 recPhyloXMLfile3 ... [-f XML files] [-s speciesFile]

                            -o fileOut              : name of the output file

                            -f XMLfiles             : (optional) name of the file containing the names of recPhyloXML files to combine
                            -s speciesFile          : (optional) name of the file containing a species tree 
                            --species.tree.newick   : (optional) whether the provided species tree is is newick format (otherwise it is assumed to be in the phyloXML format)
               """


    nextKEY = None
    params = {
                 "-o" : None, #name of the output file
                 "-f" : None, # (optional) name of the file containing the names of recPhyloXML files to combine
                 "-s": None, # (optional) name of the file containing a species tree
                 "--species.tree.newick" : False # (optional) whether the species tree is is newick format (otherwise it is assumed to be in the phyloXML format
            }

    flagArgs = ["--species.tree.newick"]

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

    if params["-o"] is None:
        print "error: output file name not given."

        print help
        exit(1)

    ## getting more files

    if not params["-f"] is None:

        if not os.path.isfile(params["-f"]):
            print "error: " + params["-f"] + " is not an existing file."
            exit(1)

        IN = open(params["-f"], "r")
        l = IN.readline()
        while l !="":
            additionalArguments.append(l.strip())
            l = IN.readline()

        IN.close()



    spLines, recLines = [],[]

    for f in additionalArguments:
        returnCode, spLines, recLinesTMP = readOneTreeFile( f )

        if returnCode != 0:
            ## tere was a problem
            if returnCode == 1:
                print "error :" , f , "is not an existing file."
            elif returnCode == 2:
                print "error :" , f , "is not a valid recPhyloXMLfile (missing a <recGeneTree> or <recPhylo> tag?)."
            else:
                print "error : unknown error no.",returnCode
        
            exit(returnCode)

        recLines += recLinesTMP


    ### include here treatment on the -s option

    if not params["-s"] is None:

        returnCode, spLines = readSpeciesTree( params["-s"] , params["--species.tree.newick"] )

        if returnCode != 0:
            if returnCode == 1:
                print "error :" , params["-s"] , "is not an existing file."
            elif returnCode == 2:
                print "error :" , params["-s"] , "is not a valid phyloXML file (missing <phylogeny> tag?)."
            else:
                print "error : unknown error no.", returnCode
            exit(returnCode)


    ### writing

    OUT = open(params["-o"],"w")

    OUT.write("<" + RECPHYLOTAG + ">" + "\n" )

    offset = 1
    offsetChar = "  "

    for l in spLines:
        OUT.write( offsetChar*offset + l )

    for l in recLines:
        OUT.write( offsetChar*offset + l )


    OUT.write("</" + RECPHYLOTAG + ">" + "\n" )
    OUT.close()