#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  
##  Created:        19-Jul-2017         
##  Last modified:  19-Jul-2017         
## 
##  Decribes one classes : ...
##
##  requires : ReconciledTree.py
##             ete3 ( http://etetoolkit.org/ )
##             xml ( in standard library )
## 
##  developped for python2.7
##
#########################################

import ete3 
import xml.etree.ElementTree as ET

from ReconciledTree import ReconciledTree, RecEvent, ReconciledTreeList




class recPhyloXML_parser:
    def __init__(self):
        pass


    def parse(self, fileName):
        """
        Takes:
            -   fileName (str) : name of a recPhyloXML file containing a single reconciled gene tree

        Returns:
            None : error
                or
            (ReconciledTree) : the reconciled tree
                or
            (ReconciledTreeList) : a set of reconciled trees
        """

        tree = ET.parse(fileName)
        
        root = tree.getroot()

        TAGtoFUNCTION = { "recPhylo" : self.parse_recPhylo,
                          "recGeneTree" : self.parse_recGeneTree }
    
        parseFunction  = TAGtoFUNCTION.get( self.tagCorrection(root.tag) , None)


        if parseFunction is None:
            raise Exception("recPhyloXML exception. Problem while parsing the xml file : no recPhylo or recgeneTree tag found at the root of the file.")
            return None

        else:
            obj = parseFunction(root)

        if obj is None:
            raise Exception("recPhyloXML exception. Problem while parsing the xml file : no phylogeny or clade found?")

        return obj


    def tagCorrection(self, tag):
        """
        Takes:
            - tag (str) : tag with or without the "{***}" prefix

        Returns:
            (str) : the tag without this prefix
        """
        return tag.rpartition("}")[2]

    def isOfTag(self, element, tag):
        """
        Takes:
            - element (Element) : an element from xml.etree.ElementTree
            - tag (str) : a tag to check
    
        Returns:
            (bool) : True if the element has the desired tag, False otherwise
        """
        if self.tagCorrection(element.tag) != tag:
            return False
        return True
    
    
    def parseSimpletextElement(self, element):
        """
        Takes:
            - element (Element) : element containing some text only
    
        Returns:
            (str) : text contained in the element
        """
        return element.text


    def parse_recPhylo(self, element):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "recPhylo" tag
    
        Returns:
            None : error
                or
            (ReconciledTreeList) : a representation of the different recGeneTrees in the file, with their species tree if it is present
        """
        TAG = "recPhylo"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        RTL = ReconciledTreeList()

        for ch in children:
            if self.isOfTag(ch,  "recGeneTree" ) :
                RT = self.parse_recGeneTree(ch)
                RTL.append(RT)

            elif self.isOfTag(ch,  "spTree" ):
                #parsing a simple tree 
                RTL.setSpTree( self.parse_SpTree(ch) ) 



        return RTL


    def parse_SpTree(self, element):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "spTree" tag
    
        Returns:
            None : error
                or
            (ete3.Tree) : the species tree
        """
        TAG = "spTree"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        node = None

        for ch in children:
            if self.isOfTag(ch,  "phylogeny" ) :
                node = self.parse_phylogeny(ch, reconciled = False)
                break
    
        return node 

    
    def parse_recGeneTree(self, element):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "recGeneTree" tag
    
        Returns:
            None : error
                or
            (ReconciledTree) : the reconciled tree
        """
        TAG = "recGeneTree"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        node = None

        for ch in children:
            if self.isOfTag(ch,  "phylogeny" ) :
                node = self.parse_phylogeny(ch, reconciled = True)
                break
    
        return node 
    
    
    def parse_phylogeny(self, element , reconciled = True):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "phylogeny" tag
            - reconciled (bool) [default = True] : whether the element passed should be considered a ReconciledTree or not
    
        Returns:
            None : error
                or
            (ReconciledTree) : the reconciled tree
                or
            (ete3.Tree) : the tree (if reconciled is True)
        """
        TAG = "phylogeny"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        node = None
    
        additionnalInfo = {}
    
        for ch in children:
            if self.isOfTag(ch ,  "clade"):
                if node is None:
                    node  = self.parse_clade(ch, reconciled)
                else:
                    raise Exception("BadTagException. A " + TAG + " element has more than one clade children (only one is expected).")                
            else:
                ### treatment for other children
                additionnalInfo[ch.tag] = ch
    
    
    
        if node is None:
            raise Exception("BadTagException. A " + TAG + " element has no clade children (one is expected).")
        
    
        ### treatment for keys
        for k,v in element.items():
            additionnalInfo[k] = v
    
        if len(additionnalInfo) > 0:
            node.add_features( **additionnalInfo )
    
    
        return node
    
    def parse_clade(self, element, reconciled = True):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "clade" tag
            - reconciled (bool) [default = True] : whether the element passed should be considered a ReconciledTree or not
    
        Returns:
            None : error
                or
            (ReconciledTree) : the reconciled tree
        """
        TAG = "clade"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        name = None
        childrenNodes = []
        events = []
    
        additionnalInfo = {}
    
        for ch in children:
            if self.isOfTag(ch ,  "clade" ):
                childrenNodes.append( self.parse_clade(ch , reconciled) )
    
            elif self.isOfTag(ch ,  "name" ):
                name = self.parseSimpletextElement(ch)            
    
            elif self.isOfTag(ch ,  "eventsRec" ):
                events = self.parse_eventsRec(ch)
    
            else:
                ### treatment for other children
                additionnalInfo[ self.tagCorrection( ch.tag ) ] = ch
    
    
        ### treatment for keys
        for k,v in element.items():
            if k != "rooted":
                additionnalInfo[k] = v

        node = None

        if reconciled:
            node = ReconciledTree()
        else:
            node = ete3.Tree()

        node.name = name

        if reconciled:
            for e in events:
                node.addEvent(e)
    
        for ch in childrenNodes:
            node.add_child( ch )
    
        if len(additionnalInfo) > 0:
            node.add_features( **additionnalInfo )
    
    
        return node
    
    def parse_eventsRec(self, element):
        """
        *recursive funtion*
    
        Takes:
            - element (Element) : element with the "eventsRec" tag
    
        Returns:
            None : error
                or
            (ReconciledTree) : the reconciled tree
        """
        TAG = "eventsRec"
    
        if not self.isOfTag(element, TAG):
            raise Exception('BadTagException. The element is of tag ' + element.tag + " instead of " + TAG + "." )
    
        children = element.getchildren()
    
        events = []
    
        for ch in children:
    
            evtCode = self.tagCorrection( ch.tag )
    
            species = None
            speciesTAGs  = ["destinationSpecies" , "speciesLocation"]
            ts = None
            tsTAG = "ts"
            additionnalInfo = {}
    
            it = ch.items()
            for k,v in it:
                if k in speciesTAGs:
                    species = v
                elif k == tsTAG:
                    ts = int(v)
                else:
                    additionnalInfo[k] = v
    
    
            evt = RecEvent(evtCode , species, ts, additionnalInfo)
    
            events.append(evt)
    
    
        return events
    

if __name__ == "__main__":

    parser = recPhyloXML_parser()
    
    fileName = '../testFiles/geneFamily0.phyloxml'
    print "Test 1: parses the file" , fileName , "and prints its structure to the screen"


    print parser.parse(fileName)


    fileName = '../testFiles/SeveralRecTrees.xml'
    print "Test 2: parses the file" , fileName , "and prints the number of reconciled trees it contains on the screen"

    RTL = parser.parse(fileName)
    print len(RTL)


    fileName = "../testFiles/genetree_SMALL.reconciled.0.ntg.xml"
    print "Test 3: parses the file" , fileName , "and prints the species tree. it contains on the screen"


    RTL = parser.parse(fileName)
    print RTL.spTree

