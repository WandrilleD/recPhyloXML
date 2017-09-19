#ifndef RECPHYLOXMLIO_H_
#define RECPHYLOXMLIO_H_


/*

This file contains helper functions to read / write recPhyloXML files

Created the: 18-09-2017
by: Wandrille Duchemin

Last modified the: 18-09-2017
by: Wandrille Duchemin

*/

#include <string>
#include <vector>
#include <map>

#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Exceptions.h>
#include <fstream>


#include "XMLUtils.h"
#include "ReconciledTree.h"

using namespace bpp;
using namespace std;

const string SPTREETAG = "spTree";
const string RECGENETREETAG = "recGeneTree";
const string CLADETAG = "clade";

// reading

void ReadRecPhyLoXMLFile(  vector <ReconciledTree>  &recTreeList, string fileName , TreeTemplate<Node> & Stree, bool replaceSpTree, unsigned verboseLevel);
void readSpTreeFromPhyloXMLFile( TreeTemplate<Node> & Stree, ifstream& fileIN, unsigned verboseLevel);
void addXMLCladeToSpTree(ifstream& fileIN, Node * node, int &nextId, bool VERBOSE);

// writing

void beginLine(ostream& OUT, int indent_level);
void WritePhyloXMLRecEvent(ostream& OUT, ReconciledTree * Rtree, Node * node, int indent_level, bool hasLoss);
void WritePhyloXMLRecTreeAux(ostream& OUT, ReconciledTree * Rtree, Node * node, int indent_level);

void WritePhyloXMLSpeTreeAux(ostream& OUT, TreeTemplate< Node > * Stree, Node * node, int indent_level);
void WritePhyloXMLSpeTree(ostream& OUT,  TreeTemplate< Node > * Stree);

void WritePhyloXML(ostream& OUT, vector <ReconciledTree> &recTreeList , TreeTemplate<Node> * Stree);
void WritePhyloXML(string fileName, vector <ReconciledTree> &recTreeList , TreeTemplate<Node> * Stree);

#endif