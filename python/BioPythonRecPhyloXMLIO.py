#!/usr/bin/python

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        25-August-2016      #
##  Last modified:  10-Apr-2018         #
#########################################

## todo
#      adapt structure to the ones in BiopythonRecPhyloXML 
#      better doc

##requires biopython
from Bio import Phylo
import Bio.Phylo.PhyloXMLIO as PhyloXMLIO
from Bio.Phylo import PhyloXML as PX

import BioPythonRecPhyloXML as BPrecPhyloXML

from xml.etree import ElementTree as ElementTree

EVENTSRECTAG = 'eventsRec'

RECGENETREETAG = 'recGeneTree'
RECPHYLOTAG = 'recPhylo'
SPPHYLOTAG = 'spTree'




class recGeneTreeXMLParser (PhyloXMLIO.Parser):
    """  A parser of reconciled gene tree in the recPhyloXML format """
    def __init__(self, file):
        PhyloXMLIO.Parser.__init__(self,file)

    _clade_recPhyloXML_list_type = { EVENTSRECTAG : 'eventsRec'}

    _clade_complex_types = ['color', 'events', 'binary_characters', 'date']
    _clade_list_types = {
        'confidence': 'confidences',
        'distribution': 'distributions',
        'reference': 'references',
        'property': 'properties',
    }

    _clade_tracked_tags = set(_clade_complex_types).union(_clade_list_types.keys()).union(
        ['branch_length', 'name', 'node_id', 'width']).union( _clade_recPhyloXML_list_type.keys() )#.union( _clade_recPhyloXML_tags )



    def _parse_clade(self, parent):
        """Parse a Clade node and its children, recursively."""

        clade = BPrecPhyloXML.Clade(**parent.attrib)
        if clade.branch_length is not None:
            clade.branch_length = float(clade.branch_length)
        # NB: Only evaluate nodes at the current level
        tag_stack = []
        for event, elem in self.context:
            namespace, tag = PhyloXMLIO._split_namespace(elem.tag)
            #print event, namespace, tag
            if event == 'start':
                if tag == 'clade':
                    clade.clades.append(self._parse_clade(elem))
                    continue
                if tag == 'taxonomy':
                    clade.taxonomies.append(self._parse_taxonomy(elem))
                    continue
                if tag == 'sequence':
                    clade.sequences.append(self._parse_sequence(elem))
                    continue
                if tag == EVENTSRECTAG: ## list of reconciliation events
                    clade.eventsRec = self._parse_eventsRec(elem)
                    continue
                if tag in self._clade_tracked_tags:
                    tag_stack.append(tag)
            if event == 'end':
                if tag == 'clade':
                    elem.clear()
                    break
                if tag != tag_stack[-1]:
                    continue
                tag_stack.pop()
                # Handle the other non-recursive children
                if tag in self._clade_list_types:
                    getattr(clade, self._clade_list_types[tag]).append(
                        getattr(self, tag)(elem))
                elif tag in self._clade_complex_types:
                    setattr(clade, tag, getattr(self, tag)(elem))
                elif tag == 'branch_length':
                    # NB: possible collision with the attribute
                    if clade.branch_length is not None:
                        raise PhyloXMLIO.PhyloXMLError(
                            'Attribute branch_length was already set '
                            'for this Clade.')
                    clade.branch_length = PhyloXMLIO._float(elem.text)
                elif tag == 'width':
                    clade.width = PhyloXMLIO._float(elem.text)
                elif tag == 'name':
                    clade.name = PhyloXMLIO._collapse_wspace(elem.text)
                elif tag == 'node_id':
                    clade.node_id = PX.Id(elem.text.strip(),
                                          elem.attrib.get('provider'))
                elif namespace != PhyloXMLIO.NAMESPACES['phy']:
                    clade.other.append(self.other(elem, namespace, tag))
                    elem.clear()
                elif tag in self._clade_recPhyloXML_list_type:
                    #clade.eventsRec = self.other(elem, namespace, tag)
                    continue
                    #getattr(clade, self._clade_recPhyloXML_list_type[tag]).append(
                    #    getattr(self, tag)(elem))
                else:
                    raise PhyloXMLIO.PhyloXMLError('Misidentified tag: ' + tag)
        return clade



    _events_recPhyloXML_tags = ['loss',
                                'transferBack',
                                'speciation',
                                'branchingOut',
                                'bifurcationOut',
                                'duplication',
                                'leaf']

    def _parse_eventsRec(self, parent):
        """ parsing reconciliation events """
        eventsRec = BPrecPhyloXML.eventsRec( None )
        for event, elem in self.context:
            namespace, tag = PhyloXMLIO._split_namespace( elem.tag )
            if event == 'end':
                if tag == EVENTSRECTAG:
                    parent.clear()
                    break
                if tag in self._events_recPhyloXML_tags:
                    eventsRec.events.append( BPrecPhyloXML.RecEvent( tag , elem.attrib ) )
                else:
                    print 'non valid tag event ' + tag
                    raise Exception('non valid tag event ' + tag)

        return eventsRec



#    def parse(self):
#        """Parse the phyloXML file incrementally and return each phylogeny."""
#        phytag = PhyloXMLIO._ns('phylogeny')
#        print "plep" , phytag
#        for event, elem in self.context:
#            if event == 'start' and elem.tag == phytag:
#                yield self._parse_phylogeny(elem)




class recPhyloXMLParser (recGeneTreeXMLParser):
    """ A parser for a set of reconciled gene tree in the recPhyloXML format """
    def __init__(self, file):
        # Get an iterable context for XML parsing events
        #function actually copied from the biopython.phyloxml parser class
        recGeneTreeXMLParser.__init__(self,file)
        #context = iter(ElementTree.iterparse(file, events=('start', 'end')))
        #event, root = next(context)
        #self.root = root
        #self.context = context
        
        if len(PhyloXMLIO._split_namespace(self.root.tag)) == 1:
            self.splittedRootTag = self.root.tag
        else:
            self.splittedRootTag = PhyloXMLIO._split_namespace(self.root.tag)[1]

    def parse(self):
        """ Parse the recPhyloXML file incrementally and return the next recGeneTree or recPhylo object """
        #print self.splittedRootTag
        if self.splittedRootTag == RECGENETREETAG :
            ## case where the root object is a recGeneTree -> 
            yield next( recGeneTreeXMLParser.parse(self) )
        elif self.splittedRootTag == RECPHYLOTAG :

            recs = BPrecPhyloXML.recPhylo()

            for event, elem in self.context:
                #print event, elem.tag, PhyloXMLIO._split_namespace(elem.tag)[1]
                tmp = PhyloXMLIO._split_namespace(elem.tag)
                SRT = elem.tag
                if len(tmp)> 1:
                    SRT = tmp[1]
                if event == 'start' and SRT == RECGENETREETAG:
                    ## reading a rec phylogeny and adding it

                    for r in recGeneTreeXMLParser.parse(self):
                        recs.addRecGeneTree( r )

                elif event == 'start' and SRT == SPPHYLOTAG:
                    ## reading a species tree
                    if recs.hasSpeciesTree(): ## there already is a species tree -> replace it but give out some error?
                        print "Warning: found several", SPPHYLOTAG, "in a", RECPHYLOTAG, ". There should be only one."
                        continue

                    stree = [x for x in  PhyloXMLIO.Parser.parse(self)]
                    print stree
                    if len(stree) != 1:
                        raise Exception("error while parsing the species tree!")

                    if not stree is None:
                        recs.setSpTree(stree[0])
            
            yield recs


    def read(self):
        """ read the recPhyloXML file return the recGeneTree or recPhylo objects it contains"""
        p =self.parse()
        X = [i for i in p]
        self.root.clear()
        self.splittedRootTag = ""
        return X


def _handle_complex(tag, attribs, subnodes, has_text=False):

    def wrapped(self, obj):
        #print tag, obj , attribs, subnodes
        elem = ElementTree.Element(tag, PhyloXMLIO._clean_attrib(obj, attribs))
        for subn in subnodes:
            if isinstance(subn, basestring):
                # singular object: method and attribute names are the same
                if getattr(obj, subn) is not None:
                    elem.append(getattr(self, subn)(getattr(obj, subn)))
            else:
                # list: singular method, pluralized attribute name
                method, plural = subn
                for item in getattr(obj, plural):
                    elem.append(getattr(self, method)(item))
        if has_text:
            elem.text = PhyloXMLIO._serialize(obj.value)
        #print "exit _handle_complex",tag, obj, type(elem)
        return elem
    wrapped.__doc__ = "Serialize a %s and its subnodes, in order." % tag
    return wrapped





class recGeneTreeXMLWriter(PhyloXMLIO.Writer):
    """Methods for serializing a recPhyloXML object to recXML."""

    def __init__(self, phylo):
        """Build an ElementTree from a PhyloXML object."""
        PhyloXMLIO.Writer.__init__(self, phylo.to_phyloxml_container())


    def phyloxml(self, obj):
        elem = ElementTree.Element(RECGENETREETAG, obj.attributes)  # Namespaces
        for tree in obj.phylogenies:
            T = self.phylogeny(tree)
            #print type(T)
            elem.append(T)
        for otr in obj.other:
            elem.append(self.other(otr))
        return elem

    def other(self, obj):
        #print "other",obj, obj.attributes
        elem = ElementTree.Element(PhyloXMLIO._ns(obj.tag, obj.namespace), obj.attributes)
        elem.text = obj.value
        #print "other done , doing children"
        for child in obj.children:
            elem.append(self.other(child))
        #print "exit other ", obj, obj.attributes
        return elem

    def RecEvent(self, obj):
    	#print "RecEvent",obj.tag, obj.attributes
        elem = ElementTree.Element(obj.tag, obj.attributes)
        elem.text = obj.value
        #print "RecEvent done", type(elem)
        return elem


    phylogeny = _handle_complex('phylogeny',
                                ('rooted', 'rerootable',
                                 'branch_length_unit', 'type'),
                                ('name',
                                 'id',
                                 'description',
                                 'date',
                                 ('confidence', 'confidences'),
                                 'clade',
                                 ('clade_relation', 'clade_relations'),
                                 ('sequence_relation',
                                  'sequence_relations'),
                                 ('property', 'properties'),
                                 ('other', 'other'),
                                 ))




    clade = _handle_complex('clade', ('id_source',),
                                ('name',
                                 'eventsRec', ## contains reconciliation events list
                                 'branch_length',
                                 ('confidence', 'confidences'),
                                 'width',
                                 'color',
                                 'node_id',
                                 ('taxonomy', 'taxonomies'),
                                 ('sequence', 'sequences'),
                                 'events',
                                 'binary_characters',
                                 ('distribution', 'distributions'),
                                 'date',
                                 ('reference', 'references'),
                                 ('property', 'properties'),
                                 ('clade', 'clades'),
                                 ('other', 'other'),
                                 ))

    eventsRec = _handle_complex('eventsRec', (), ( ('RecEvent','events'), ) )
                                


class recPhyloXMLWriter( recGeneTreeXMLWriter ):
    """
    Methods for serializing a recPhylo object into an XML format
    """
    def __init__(self,recs):
        """ builds an element tree from a recPhylo object """
        if isinstance(recs, BPrecPhyloXML.recPhylo ) :
            self._tree = ElementTree.ElementTree(self.recPhylo(recs))
        elif isinstance(recs, Phylo.PhyloXML.Phylogeny ) :
            recGeneTreeXMLWriter.__init__(self,recs)
        else:
            print "not a recPhylo nor a reconciled gene tree."

    def recPhylo(self, obj):
        elem = ElementTree.Element( RECPHYLOTAG )
        if obj.hasSpeciesTree():
            elem.append( self.namedPhyloContainer( obj.spTree , SPPHYLOTAG ) )
        for rectree in obj:
            T = self.namedPhyloContainer( rectree )
            #print type(T)
            elem.append(T)
        return elem

    def namedPhyloContainer(self, obj, name = RECGENETREETAG):
        elem = ElementTree.Element(name)  # Namespaces
        T = self.phylogeny(obj)
        elem.append(T)
        return elem



def read(file):
    """
    Read a recPhyloXML file and build 
    a biopython tree
    or a recPhylo object (that contains a certain number of reconciled gene tree as biopython tree objects)

    Takes:
        - file (str) : the name of a file containing recPhyloXML data 

    Returns:
        ( Bio.phylo.PhyloXML.Phylogeny or BioPythonRecPhyloXML.recPhylo )

    """
    return recPhyloXMLParser(file).read()

def parse(file):
    """
    Parse a recPhyloXML file and build 
    a biopython tree
    or a recPhylo object (that contains a certain number of reconciled gene tree as biopython tree objects)

    Takes:
        - file (str or input handle) : the name of a file containing recPhyloXML data 

    Returns:
        (generator) : a generator of Bio.phylo.PhyloXML.Phylogeny or BioPythonRecPhyloXML.recPhylo objects

    """
    return recPhyloXMLParser(file).parse()

def write(obj,file):
    """
    write a reconciled gene tree (or an ensemble of them) in the recPhyloXML format

    Takes:
        - obj (Bio.phylo.PhyloXML.Phylogeny or BioPythonRecPhyloXML.recPhylo) : object to write 
        - file (str or output handle) : filename or open output handle where the xml data will be written
    """
    return recPhyloXMLWriter(obj).write(file)



if __name__ == "__main__":
#    filename = '../testFiles/geneFamily0.phyloxml'
#
#    p = recPhyloXMLParser(filename).read()[0]
#    print type(p)
#    
#
#    leaf = p.get_terminals()[0]
#    print leaf.eventsRec ,type(leaf.eventsRec)
#
#    for x in  leaf.eventsRec.events:
#        print type(x), x.tag, x.attributes

    

    print "testing Biopython recPhyloXML API"

    filename2 = '../testFiles/lossSeparatedtestFile'
    print "parsing reconciled gene tree from",filename2

    X = recPhyloXMLParser(filename2)
    p = X.read()[0]

    print "generated object:", [type(i) for i in p]
    print "number of gene trees in the file",len(p)

    print p.hasSpeciesTree()
    print p.spTree

    #leaf = p[0].get_terminals()[0]
    #print leaf.eventsRec
    #for x in  leaf.eventsRec.events:
    #    print x.tag, x.attributes


    write(p,'../testFiles/generatedRecPhylo.xml')
    print "written generated recPhylo object in", '../testFiles/generatedRecPhylo.xml'
    write(p[0],'../testFiles/generated.xml')
    print "written generated recGeneTree object in", '../testFiles/generated.xml'
