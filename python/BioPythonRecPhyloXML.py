#!/usr/bin/python

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        25-August-2016      #
##  Last modified:  10-Apr-2018         #
#########################################


#todo:
#     make specific grec gene tree class with util fc
#     make recphylo object taking a (facultative) sp tree and a list of gene trees


import Bio.Phylo.PhyloXML as PhyloXML


#class recPhylogeny(PhyloXML.Phylogeny):
#
#    def __init__(self, root=None, rooted=True,
#                 rerootable=None, branch_length_unit=None, type=None,
#                 # Child nodes
#                 name=None, id=None, description=None, date=None,
#                 # Collections
#                 confidences=None, clade_relations=None, sequence_relations=None,
#                 properties=None, other=None,
#                 ):
#
#                PhyloXML.Phylogeny.__init__(self, root, rooted, 
#                                    rerootable, branch_length_unit, type, 
#                                    name, id, description, date,
#                                    confidences, clade_relations, sequence_relations,
#                                    properties, other
#                                    )

class recPhylo:
    """
    a structure for reconciled phylogeny data as described in the recPhyloXML format (url coming)

    Attributes:
        - spTree : Bio.Phylo.PhyloXML.Phylogeny or None
                            Represents the species tree the gene trees are reconciled with. 
                            This may be left to None to indicated that the species tree data is absent
        - recGeneTrees : list of Bio.Phylo.PhyloXML.Phylogeny instances
                            The clades of these instances have an eventsRec (a class described below) that contains reconciliation information
    """
    def __init__(self, spTree = None , recGeneTrees = []):
        """
        Takes:
            - spTree [default = None] (Bio.Phylo.PhyloXML.Phylogeny or None) :
                            Represents the species tree the gene trees are reconciled with. 
                            This may be left to None to indicated that the species tree data is absent
            - recGeneTrees [ default = [] ] (list) : list of Bio.Phylo.PhyloXML.Phylogeny instances
                            The clades of these instances have an eventsRec (a class described below) that contains reconciliation information
        """
        self.spTree = spTree
        self.recGeneTrees = recGeneTrees

    def setSpTree(self,newspTree):
        """
        set the species tree to a new one

        Takes:
            - newspTree (Bio.Phylo.PhyloXML.Phylogeny) : the new species tree that the gene trees are reconciled on
        """
        self.spTree = newspTree

    def addRecGeneTree(self,tree):
        """
        adds a reconciled gene tree 

        Takes:
            - tree (Bio.Phylo.PhyloXML.Phylogeny): reconciled gene tree to add
        """
        self.recGeneTrees.append(tree)

    def __len__(self):
        """ returns the number of reconciled gene tree """
        return len(self.recGeneTrees)

    def hasSpeciesTree(self):
        """ True iif the species tree is set """
        return not self.spTree is None

    def __getitem__(self,i):
        """ returns the reconciled gene tree at the corresponding index (if it exists) """
        return self.recGeneTrees[i]

    def __iter__(self):
        for x in self.recGeneTrees:
            yield x

class Clade(PhyloXML.Clade):
    """ a Clade class for Biopython.Phylo trees with a list of reconciliation events"""
    def __init__(self,
                 # Attributes
                 branch_length=None, id_source=None,
                 # Child nodes
                 name=None, width=None, color=None, node_id=None, events=None,
                 binary_characters=None, date=None,
                 # Collections
                 confidences=None, taxonomies=None, sequences=None,
                 distributions=None, references=None, properties=None, clades=None,
                 other=None,
                 eventsRec = None ## list of reconciliation events
                 ):
        PhyloXML.Clade.__init__(self,
                 # Attributes
                 branch_length, id_source,
                 # Child nodes
                 name, width, color, node_id, events,
                 binary_characters, date,
                 # Collections
                 confidences, taxonomies, sequences,
                 distributions, references, properties, clades,
                 other,
                 )

        self.eventsRec = eventsRec or [] ## list of reconciliation events



class eventsRec (PhyloXML.PhyloElement):
    """ a object containing a list of reconciliation events """
    def __init__(self, events = None):
        self.tag = 'eventsRec'
        self.events = events or []

    def __str__(self):
        return ",".join([e.__str__() for e in self.events])


class RecEvent(PhyloXML.PhyloElement):

    ok_tag = ['loss',
              'transferBack',
              'speciation',
              'branchingOut',
              'bifurcationOut',
              'duplication',
              'leaf']

    def __init__(self,tag , attr):
        PhyloXML._check_str(tag, self.ok_tag.__contains__)
        self.tag = tag
        self.attributes = attr or []
        self.value=""
