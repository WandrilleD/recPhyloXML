#ifndef RECONCILED_TREE_H_
#define RECONCILED_TREE_H_


/*

This file contains a class for reconciled trees

Created the: 15-09-2017
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

using namespace bpp;

//Variables used to get and set some node properties
const string spe = "S";  //name of the node property containing their species
const string ev = "Ev";  //name of the node property containing their event
const string clnum = "Cnum";//name of the node property containing their clade id
const string porder = "Po";//name of the node property containing their post order
//const string alpha = "Al";//name of the node property containing their alpha status (whether they are in a dead lineage or no)
const string ts = "TS";//name of the node property containing their time slice
const string uts = "UTS";//upper time slice boundary
const string lts = "LTS";//lower time slice boundary


//Code of the reconciliation events

const int C     = 0; //Current (leaf)
const int S     = 1; //Speciation (in the species tree)
const int L     = 2; //Loss (in the species tree)
const int D     = 3; //Duplication (in the species tree)
const int Sout  = 4; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
const int R     = 5; //Transfer reception
const int N     = 6; //no event (to account for time slices)
const int Bout  = 7; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)


class ReconciledTree: public TreeTemplate<Node>
{
protected:

    map<int,vector<int> > CladeIdToNodeIds;
    //map<int, Node * > NodeIdToNodeP;


    int nbNodes;
    int MyNextId()
    {
        return nbNodes;
    }

    int TimeSliceStatus; // 0 if no time slices (NTS) , 1 if precise time slices (TS)

    int countAnyEvent( int evttype);


    /*
    Takes:
        - Nodeid (int): a node id
        - evtLine (string): recPhyloXML line describing a reconciliation event
    */
    void setNodeDetailsFromXMLLine(int Nodeid, string evtLine, bool VERBOSE);

    //same as previous; uses the interpreted xml line instead
    void setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value,  bool VERBOSE);

    /*
    read and add a clade from a recPhyloXML stream
    
    Takes:
        - fileIN (ifstream): streaam to a recPhyloXML file
        - VERBOSE (bool) (default = false)

    */
    void addXMLClade(ifstream& fileIN, int Nodeid, bool VERBOSE = false); // will need a species tree


    void readPhyloXMLFile(ifstream& fileIN, TreeTemplate<Node> * Stree, bool VERBOSE );

    void removeNoEventNodes();


    string NodeString(int nodeid, bool hideLosses = false); //pseudo newick representation of the node

public:
    ReconciledTree() {} // empty tree

    ReconciledTree(Node * root, int TSstatus);


    //Constructor that uses a stream from a recPhyloXML file
    ReconciledTree(ifstream& fileIN, TreeTemplate<Node> * Stree, bool VERBOSE = false);

    //Constructor that uses a recPhyloXML file name
    ReconciledTree(string phyloxmlFileName, TreeTemplate<Node> * Stree, bool VERBOSE = false);


    ~ReconciledTree(){}//fairly simple destruction



    //access to NodeId with CladeId
    vector<int> getNodeIdsFromCladeId(int CladeId);
    void addNodeIdToCladeIdMap(int CladeId, int NodeId);


    int getTimeSliceStatus();

    //getter of node properties
    string getNodeSpecies(int nodeid); // retrieve the specified node species. 
    int getNodePostOrder(int nodeid);
    int getNodeTimeSlice(int nodeid);
    int getNodeEvent(int nodeid); // retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
    int getNodeCladeNum(int nodeid); // retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance

    string getNodeSpecies(Node * NP); // retrieve the specified node species. 
    int getNodePostOrder(Node * NP);
    int getNodeTimeSlice(Node * NP);
    int getNodeEvent(Node * NP); // retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
    int getNodeCladeNum(Node * NP); // retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance


    //setter of node properties
    void setNodeSpecies(int nodeid, string speciesid); // set the specified node species.
    void setNodePostOrder(int nodeid, int postorder);
    void setNodeTimeSlice(int nodeid, int timeslice);
    //void setNodeUpperBoundaryTS(int nodeid, int timeslice);
    //void setNodeLowerBoundaryTS(int nodeid, int timeslice);
    void setNodeEvent(int nodeid, int eventid); // assigns an event to the node according to the eventid
    void setNodeCladeNum(int nodeid, int cladenum); // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance

    void setNodeSpecies(Node * NP, string speciesid); // set the specified node species.
    void setNodePostOrder(Node * NP, int postorder);
    void setNodeTimeSlice(Node * NP, int timeslice);
    //void setNodeUpperBoundaryTS(Node * NP, int timeslice);
    //void setNodeLowerBoundaryTS(Node * NP, int timeslice);
    void setNodeEvent(Node * NP, int eventid); // assigns an event to the node according to the eventid
    void setNodeCladeNum(Node * NP, int cladenum); // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance



    void resetNodeCladeNum(int nodeid);//unset the clade id of the node. Always use this function because it also updates the CladeIdToNodeIds map

    void printNode(int nodeid);//prints node information
    void printMe();//prints all nodes
    string NewickString(bool hideLosses = false);


    //void computeBoundedTimeSlices(Tree * SpeciesTree);
    void setToNonTimeSliced();
    void SubdivideTree();//expects that node already have an exact set time slice

    int getNumberOfDuplication();
    int getNumberOfLoss();
    int getNumberOfTransfer();


    bool haveSameSpecies(Node  * n1, Node * n2);// just check species
    bool areTSCompatible(Node  * n1, Node * n2);// just check TS
    bool areCompatible(Node  * n1, Node * n2); //check if the node are compatible (ie. same species and comparable timeslice). Presumes that both nodes have the same timeslice status as the tree
    bool areCompatible(int id1, int id2); //check if the node are compatible (ie. same species and comparable timeslice). Presumes that both nodes have the same timeslice status as the tree


    bool isAncestor(int AncId, int SonId);//Checks if AncId is the ancestor of SonId
    vector <Node *> pathFromAncestorToSon(int AncId, int SonId); // Get the path from ancestor to son.
    vector <int> idPathFromAncestorToSon(int AncId, int SonId); // Get the path from ancestor to son.

    int getIdWithName(const string &name);


    bool isRealLeaf(int id);
    bool isExtant(int evtcode);
    bool isSpeciation(int evtcode);
    bool isLoss(int evtcode);
    bool isDup(int evtcode);
    bool isSout(int evtcode);
    bool isRec(int evtcode);
    bool isNull(int evtcode);
    bool isBout(int evtcode);

    ReconciledTree * cloneSubtree(int newRootId);

    vector <string> getRealLeavesNames();

};



#endif /*RECONCILED_TREE_H_*/