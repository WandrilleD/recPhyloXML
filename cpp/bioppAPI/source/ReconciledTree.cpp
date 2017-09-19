
/*

This file contains a class for reconciled trees

Created the: 15-09-2017
by: Wandrille Duchemin

Last modified the: 15-09-2017
by: Wandrille Duchemin

*/



#include "ReconciledTree.h"




/*
recursively searches a node with the correct name.

Takes:
    - Node * n : current node of the recursion
    - string name : desired name

Returns:
    (Node * ) pointer to desireed node or Null if such a node is absent.

*/
Node * getNodeFromName(Node * n , string name)
{
    if( n->getName() == name )
        return n;

    Node * res = NULL;
    vector <Node * > children =   n->getSons();
    for(unsigned i = 0 ; i < children.size() ; i++)
    {
        res = getNodeFromName(children[i] , name);
        if(res != NULL)
            return res;
    }


    return res;
}



//Variables used to get and set some node properties
/*const string spe = "S";  //name of the node property containing their species
const string ev = "Ev";  //name of the node property containing their event
const string clnum = "Cnum";//name of the node property containing their clade id
const string porder = "Po";//name of the node property containing their post order
//const string alpha = "Al";//name of the node property containing their alpha status (whether they are in a dead lineage or no)
const string ts = "TS";//name of the node property containing their time slice
*/



int ReconciledTree::countAnyEvent( int evttype)
{
    int counter = 0;
    vector<int> NodesIds  = getNodesId();
    int NodeId;

    for (unsigned j=0; j<NodesIds.size() ; j++)
    {
        NodeId = NodesIds[j];
        if(getNodeEvent(NodeId) == evttype)
            counter++;
    }
    return counter;
}



/*
Takes:
    - Nodeid (int): a node id
    - evtLine (string): recPhyloXML line describing a reconciliation event
*/
void ReconciledTree::setNodeDetailsFromXMLLine(int Nodeid, string evtLine, bool VERBOSE)
{

    map <string, string> * properties = new map <string,string>;
    string * value = new string(""); 
    string Tname = InterpretLineXML(evtLine,properties,value);

    setNodeDetailsFromXMLLine( Nodeid,  Tname,  properties,  value, VERBOSE);
}


/*
Takes:
    - Nodeid (int): a node id
    - Tname (string) : name of the tag (here, event name)
    - properties (map <string, string> *): container for informations like species, timeSlice, ...
    - value ( string * ): should be empty, technically not used here

*/
void ReconciledTree::setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value, bool VERBOSE)
{
    char * pEnd; // for str -> int conversion 
    
    
    for (map<string,string>::iterator it=properties->begin(); it!=properties->end(); ++it)
    {
        if(it->first.compare("timeSlice") == 0)
        {
            int TS =  (int) strtol(it->second.c_str(), &pEnd, 10) ;
            setNodeTimeSlice(Nodeid,TS);
        }
        else if (it->first.compare("speciesLocation") == 0)
            setNodeSpecies(Nodeid, it->second);
        else if (it->first.compare("destinationSpecies") == 0)
            setNodeSpecies(Nodeid, it->second);
        else if (it->first.compare("geneName") == 0)
            setNodeName(Nodeid,it->second);
    }


    /*
    C       --> Current (leaf)
    S       --> Speciation (in the species tree)
    L       --> Loss (in the species tree)
    D       --> Duplication (in the species tree)
    Sout    --> Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
    R       --> Transfer reception
    N       --> no event (to account for time slices)
    Bout    --> Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
    */

    int evtCode = -1;
    //3. finding which event is implied
    if(Tname.compare("leaf") == 0)
    {
        setNodeTimeSlice(Nodeid,0);

        evtCode = C;
    }
    else if(Tname.compare("speciation") == 0)
    {
        evtCode = S;
    }
    else if(Tname.compare("speciationOut") == 0)
    {
        evtCode = Sout;
    }
    else if(Tname.compare("bifurcationOut") == 0)
    {
        evtCode = Bout;
        setNodeSpecies(Nodeid, "-1" );
    }
    else if(Tname.compare("duplication") == 0)
    {
        evtCode = D;
    }
    else if(Tname.compare("transferBack") == 0)
    {
        evtCode = R;
    }
    else if(Tname.compare("speciationOutLoss") == 0)
    {
        evtCode = Sout;
        //create some Loss son in the same species, with the same time slice

        int newId = getNextId();
        Node * newLossNode  = new Node(newId);
        Node * currentNode  = getNode(Nodeid);
        currentNode->addSon(newLossNode);

        //setting new loss node properties
        setNodeSpecies(newId,getNodeSpecies(Nodeid));
        setNodeEvent(newId, L); // setting as a LOSS
        setNodeCladeNum(newId, -1); // the clade num of the loss node is -1 -> the empty clade

    }
    else if(Tname.compare("speciationLoss") == 0)
    {
        evtCode = S;
        //create some Loss son, but in a unknown species... with time slice = current -1
        // The loss node will be added in a later post-treatment phase as it requires a species tree to determine in which species the loss is. 
        // This won't get forgotten as these are speciations nodes with only one son

    }

    if(evtCode == -1)
        throw Exception("ReconciledTree::setNodeDetailsFromXMLLine : unknown event : " + Tname );
    else
        setNodeEvent(Nodeid,evtCode);

    if(VERBOSE)
    {
        cout << "set details of node " << Nodeid << ";" << evtCode << ";" << getNodeSpecies(Nodeid) << ";" ;
        if(hasNodeProperty(Nodeid,ts))
            cout <<getNodeTimeSlice(Nodeid) << endl;
        else
            cout << "NOTS" << endl;
    }
        

}

/*
read and add a clade from a recPhyloXML stream

!!! This function is recursive

Takes:
    - fileIN (ifstream): streaam to a recPhyloXML file
    - Nodeid (int) : node to put the clade IN
    - VERBOSE (bool) (default = false)
*/
void ReconciledTree::addXMLClade(ifstream& fileIN,int Nodeid, bool VERBOSE) // will need a species tree
{

    int Cnum = CladeIdToNodeIds.size(); //the clade num will be some equivalent to a BFO index

    setNodeCladeNum(Nodeid , Cnum);

    if(VERBOSE)
        cout << "parsing clade " << Cnum << " id " << Nodeid<< endl;

    vector <string> evtLines;

    //2. parsing of the clade
    string line;
    getline( fileIN, line );

    map <string, string> * properties = new map <string,string>;
    string * value = new string(""); 
    string Tname = InterpretLineXML(line,properties,value);

    //cout << "value:" << *value << endl;

    string name = "";

    while(Tname.compare("/clade") != 0) //--> line signing the end of the clade
    {
        if(Tname.compare("name") == 0) // found name
        {
            if(VERBOSE)
                cout << "found name : " << *value <<  endl;
            name = *value;

        }
        else if(Tname.compare("eventsRec") == 0) // reading reconciliation events
        {
            if(VERBOSE)
                cout <<  "parsing events" << endl;
            
            *value = "";
            delete properties;

            getline( fileIN, line );
            properties = new map <string,string>;
            Tname = InterpretLineXML(line,properties,value);

            while(Tname.compare("/eventsRec") != 0) // here, add some nodes ...
            {
                evtLines.push_back(line); // maybe give in the already interpreted line rather than re-interpret it inside constructor...
                
                //looping
                *value = "";
                delete properties;
                getline( fileIN, line );
                properties = new map <string,string>;
                Tname = InterpretLineXML(line,properties,value);

                if(Tname.compare("") == 0)
                    throw Exception("ReconciledTree::addXMLClade : error while parsing...");
            }
        }
        else if(Tname.compare("clade") == 0) // found a child --> recursion
        {
            //create child clade

            //creating the new node
            int newId = getNextId();

            Node * newNode  = new Node(newId);

            //branching to the existing structure
            Node * currentNode  = getNode(Nodeid);
            currentNode->addSon(newNode);

            //recursing
            addXMLClade(fileIN, newId,  VERBOSE );


        }

        //reading new line
        *value = "";
        delete properties;
        getline( fileIN, line );
        properties = new map <string,string>;
        Tname = InterpretLineXML(line,properties,value);
        
    }



    //3. Event parsing
    if(evtLines.size() == 0)
    {
        throw Exception("ReconciledTree::addXMLClade : no reconciliation events found");
    }
    else if(evtLines.size()>1)
    {
        for(unsigned i = 0; i < evtLines.size() - 1; i++ ) // all events but the last
        {
            //creating parent node 
            //cout << "creating evt" << evtLines[i] << endl;

            int newId = getNextId();

            Node * newNode  = new Node(newId);

            //branching to the existing structure
            Node * currentNode  = getNode(Nodeid);

            if(!isRoot(Nodeid))
            {
                Node * FNode = currentNode->getFather();
                FNode->removeSon(currentNode);
                FNode->addSon(newNode);
            }
            newNode->addSon(currentNode);

            if((i==0)&&(name!=""))
                newNode->setName(name);

            if(isRoot(Nodeid))
                rootAt(newNode); // rerooting

            //new putting correct attributes to in the new node
            setNodeCladeNum(newId , Cnum);

            setNodeDetailsFromXMLLine(newId,evtLines[i], VERBOSE);


        }
    }
    else if(name!="")
    {
        setNodeName(Nodeid,name);
    }


    setNodeDetailsFromXMLLine(Nodeid,evtLines.back(), VERBOSE); // setting details on the last node == the bifurcation or leaf


}

/*
Constructor that uses a stream from a recPhyloXML file
Takes
 - fileIN (ifstream): stream rom a recPhyloXML file
 - Stree (TreeTemplate<Node> *): a pointer to the species tree
 - VERBOSE (bool) [fdefault : false]
*/
void ReconciledTree::readPhyloXMLFile(ifstream& fileIN, TreeTemplate<Node> * Stree, bool VERBOSE )
{
    // presuming that the stream will yield a <clade> line froming the root of the tree
    //first find the root
    if(VERBOSE)
        cout << "Creating Reconciled_Tree." << endl;

    int rootid = 0;

    Node * newnode  = new Node(rootid);

    setRootNode(newnode);//setting as the new root


    addXMLClade( fileIN,rootid,  VERBOSE ); 

    ///post-treatment requiring species tree
    vector < Node * > NodesP  = getNodes();


    bool NoSpeciesTree = false;
    if(Stree == NULL)
        NoSpeciesTree=true;
    else if(Stree->getRootNode() == NULL)
        NoSpeciesTree=true;

    for(unsigned i = 0 ; i < NodesP.size(); i++)
    {
        int evt = getNodeEvent(NodesP[i]);
        vector < Node * > SonsP = NodesP[i]->getSons();

        if(( isSpeciation( evt ) ) && (SonsP.size() == 1)) //speciation with only one child --> speciation loss
        {
            string lostSpecies;
            if(NoSpeciesTree)
            { // failsafe for the case where there is no valid species tree
                lostSpecies = "";
            }
            else
            { 
                string spF = getNodeSpecies(NodesP[i]); //species of the speciation
                string spC = getNodeSpecies(SonsP[0]); //sister species of the loss.
    
                Node * ParentSpNode = Stree->getNode(spF);
    
                if( ParentSpNode == NULL)
                    throw Exception("ReconciledTree::readPhyloXMLFile : could not find species with name "+spF+"."); 
    
    
                vector <Node * > spSons = ParentSpNode->getSons();
    
                //there should be 2 sons in spSons
                if(spSons.size() != 2)
                    throw Exception("ReconciledTree::readPhyloXMLFile : species without son as a speciation."); 
    
                int LostIndex = 0;
                if(spC == spSons[LostIndex]->getName())
                    LostIndex++;

                lostSpecies = spSons[LostIndex]->getName();
            }


            if(VERBOSE)
                cout << "Adding Loss node in species " << lostSpecies << endl;


            //creating new loss node
            int newId = getNextId();
            Node * newLossNode  = new Node(newId);
            Node * currentNode  = NodesP[i];
            currentNode->addSon(newLossNode);
    
            //setting new loss node properties
            setNodeSpecies(newLossNode, lostSpecies );
            setNodeEvent(newLossNode, L); // setting as a LOSS
            setNodeCladeNum(newLossNode, -1); // the clade num of the loss node is -1 -> the empty clade

        }
        else if((evt == S) && (SonsP.size() == 0))
            throw Exception("ReconciledTree::readPhyloXMLFile : speciation without any son.");
    }


    /// determining timeSlice Status
    if(hasNodeProperty(getRootId(),ts))
        TimeSliceStatus = 1;
    else
        TimeSliceStatus = 0;

}

/* 
* *RECURSIVE*
* @arg int nodeid : id of the node whose subtree string is built
* @arg bool hideLosses (default: false): if true, losses and the branches leading to them wil be removed from the newick string
*
* @return string : newick representation of the subtree rooted at nodeid 
*/
string ReconciledTree::NodeString(int nodeid, bool hideLosses)
{
    //to rewrite using pointers instead of ids...
    string Str = "";

    vector<int> sons = getSonsId(nodeid);
    if(sons.size()> 0)
    {

        vector <string> sonStr;

        for(unsigned i = 0; i < sons.size(); i++)
        {
            string tmp =  NodeString(sons[i], hideLosses); // recursion
            if(tmp.size() >0) // don't add it the the representation is an empty str
                sonStr.push_back(tmp);
            
        }

        
        if(sonStr.size() == 0) // no son with a non-empty representation
            return Str; // return empty string

        //building Str
        Str += "(";
        for(unsigned i = 0; i < sonStr.size(); i++)
        {   
            if(i != 0)
                Str += ",";
            Str += sonStr[i];
        }
        Str += ")";
    }

    if(hasNodeName(nodeid))
        Str += getNodeName(nodeid);
    else
    {
        Str += static_cast<ostringstream*>( &(ostringstream() << nodeid) )->str() ;
    }
        

    Str += "|";

    if(!hasNodeProperty(nodeid,ev))
        Str += "NOEVENT";
    else
    {

        int evt = getNodeEvent(nodeid);
        switch(evt)
        { 
            case C :
                Str += "Extant"; //Current (leaf)
                break;
            case S :
                Str += "Spe"; //Speciation (in the species tree)
                break;
            case L :
                if(hideLosses)// soecial option so that losses and their branch oare not in the final newick tree
                {
                    return ""; // return empty string
                }
                Str += "Loss"; //Loss (in the species tree)
                break;
            case D :
                Str += "Dup"; //Duplication (in the species tree)
                break;
            case Sout :
                Str += "SpeOut"; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
                break;
            case R :
                Str += "Reception"; //Transfer reception
                break;
            case N :
                Str += "Null"; //no event (to account for time slices)
                break;
            case Bout :
                Str += "BifOut"; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
                break;
            default:
                Str += "NOEVENT";
        }
    }

    Str += "|";

    if(!hasNodeProperty(nodeid,spe))
        Str += "no species";
    else
    {
        Str += getNodeSpecies(nodeid);
    }

    Str += "|";

    if(hasNodeProperty(nodeid,ts))
    {
        int TS = getNodeTimeSlice(nodeid);
        Str += static_cast<ostringstream*>( &(ostringstream() << TS) )->str() ;

    }
    else
        Str += "NOTS";

    return Str;
}

/*
Constructor from a root and a time slice status.

Intended to be used by CloneSubtree
*/
ReconciledTree::ReconciledTree(Node * root, int TSstatus): TreeTemplate<Node>(root)
{
    TimeSliceStatus = TSstatus;
}



/*
Constructor that uses a stream from a recPhyloXML file
Takes
 - fileIN (ifstream): stream rom a recPhyloXML file
 - Stree (TreeTemplate<Node> *): a pointer to the species tree, where the names of the nodes will correspond to the speciesLocation property of the events in the XML file
 - VERBOSE (bool) [fdefault : false]
*/
ReconciledTree::ReconciledTree(ifstream& fileIN, TreeTemplate<Node> * Stree, bool VERBOSE): TreeTemplate<Node>()
{

    readPhyloXMLFile( fileIN, Stree, VERBOSE );
}

/*
Constructor that uses a recPhyloXML file name

Note that this constructor is VERY crude and will try to build the reconciled tree
out of the first <clade> tag it sees, independently of the presence of other reconciled gene tree or species tree.

One should only use this constructor on files where one is sure that there is a single reconciled gene tree without any species tree.


Takes
 - string phyloxmlFileName:  recPhyloXML file name
 - TreeTemplate<Node> (Tree *): a pointer to the species tree, where the names of the nodes will correspond to the speciesLocation property of the events in the XML file
 - VERBOSE (bool) [fdefault : false]
*/
ReconciledTree::ReconciledTree(string phyloxmlFileName, TreeTemplate<Node> * Stree, bool VERBOSE )
{
    ifstream fileStream(phyloxmlFileName.c_str());
    if( !fileStream.is_open() ) 
    {
        throw Exception("ReconciledTree::ReconciledTree : could not open reconciled tree file : " + phyloxmlFileName);
        exit(1);
    }

    while( !fileStream.eof() ) 
    {

        string line;
        getline( fileStream, line );
        
        map <string, string> * properties = new map <string,string>;
        string * value; 
        string Tname = InterpretLineXML(line,properties,value);

        if(Tname.compare("clade") == 0)
        {
            readPhyloXMLFile(fileStream,Stree, VERBOSE);
            break;
        }

    }


}


/*
Returns 0 if there is no time slices (NTS), 1 if there are time slices (TS) or 2 if there are bounded time slices (BTS)
*/
int ReconciledTree::getTimeSliceStatus()
{
    return TimeSliceStatus;
}

/*
Access to NodeId with CladeId
Takes:
 - CladeId (int): id of a clade

Returns:
    vector<int>: NodeIds of the nodes with that clade OR empty vector if the clade is absent

*/
vector<int> ReconciledTree::getNodeIdsFromCladeId(int CladeId)
{
    if(CladeIdToNodeIds.count(CladeId) == 0)
        return vector<int>();
    return CladeIdToNodeIds[CladeId];
}



/*
Adds NodeId to CladeId
Takes:
 - CladeId (int): id of a clade
 - NodeId (int):  if of the node with that clade

*/
void ReconciledTree::addNodeIdToCladeIdMap(int CladeId, int NodeId)
{
    if(CladeIdToNodeIds.count(CladeId) == 0)//key is absent from the map -> add it
        CladeIdToNodeIds[CladeId] = vector<int>();
    CladeIdToNodeIds[CladeId].push_back(NodeId);
}


//getter of node properties

/*
retrieve the specified node species. 

Takes:
 - NodeId (int): id of a node in the tree

Returns:
    (string): speciesid of the node
*/
string ReconciledTree::getNodeSpecies(int nodeid)
{

    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(spe) )->getValue();
    }
    //else we update the map

    NodeIdToNodeP[nodeid] = getNode(nodeid);*/


    return dynamic_cast<BppString *> ( getNodeProperty(nodeid, spe))->toSTL();
}

int ReconciledTree::getNodePostOrder(int nodeid)
{
    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(porder) )->getValue();
    }
    //else we update the map
    NodeIdToNodeP[nodeid] = getNode(nodeid);*/

    return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, porder))->getValue();
}


/*
Will return -1 of there are no time slices
*/
int ReconciledTree::getNodeTimeSlice(int nodeid)
{
    if(TimeSliceStatus == 1)
        return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, ts))->getValue();
    
    return -1;
}



int ReconciledTree::getNodeEvent(int nodeid) // retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
{

    /*map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(ev) )->getValue();
    }
    //else we update the map
    NodeIdToNodeP[nodeid] = getNode(nodeid);*/

    return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, ev))->getValue();
}

int ReconciledTree::getNodeCladeNum(int nodeid) // retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{

    /*map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(clnum) )->getValue();
    }
    //else we update the map
    NodeIdToNodeP[nodeid] = getNode(nodeid);*/

    return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, clnum))->getValue();
}




//retrieve the specified node species. 
string ReconciledTree::getNodeSpecies(Node * NP)
{
    return dynamic_cast<BppString *> (NP->getNodeProperty(spe))->toSTL();
}

//;
int ReconciledTree::getNodePostOrder(Node * NP)
{
    return dynamic_cast<BppInteger *> (NP->getNodeProperty(porder))->getValue();
}

/*
Will return -1 of there are no time slices
*/
int ReconciledTree::getNodeTimeSlice(Node * NP)
{
    if(TimeSliceStatus == 1)    
        return dynamic_cast<BppInteger *> (NP->getNodeProperty(ts))->getValue();
    
    return -1;
}

// retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
int ReconciledTree::getNodeEvent(Node * NP)
{
    return dynamic_cast<BppInteger *> (NP->getNodeProperty(ev))->getValue();
}

//retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
int ReconciledTree::getNodeCladeNum(Node * NP)
{
    return dynamic_cast<BppInteger *> (NP->getNodeProperty(clnum))->getValue();
}











//setter of node properties
void ReconciledTree::setNodeSpecies(int nodeid, string speciesid) // set the specified node species.
{
    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        it->second->setNodeProperty(spe, BppInteger(speciesid) );
    }
    else
    {
        //else we update the map
        Node * NP = getNode(nodeid);
        NodeIdToNodeP[nodeid] = NP;

        NP->setNodeProperty(spe, BppInteger(speciesid));
    }*/
    Node * NP = getNode(nodeid);
    NP->setNodeProperty(spe, BppString(speciesid));
}

void ReconciledTree::setNodeSpecies(Node * NP, string speciesid) // set the specified node species.
{
    NP->setNodeProperty(spe, BppString(speciesid));    
}

void ReconciledTree::setNodePostOrder(int nodeid, int postorder)
{
    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        it->second->setNodeProperty(porder, BppInteger(postorder) );
    }
    else
    {
        //else we update the map
        Node * NP = getNode(nodeid);
        NodeIdToNodeP[nodeid] = NP;

        NP->setNodeProperty(porder, BppInteger(postorder));
    }*/

    Node * NP = getNode(nodeid);
    NP->setNodeProperty(porder, BppInteger(postorder));

}

void ReconciledTree::setNodePostOrder(Node * NP, int postorder)
{
    NP->setNodeProperty(porder, BppInteger(postorder)); 
}

void ReconciledTree::setNodeTimeSlice(int nodeid, int timeslice)
{
    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        it->second->setNodeProperty(ts, BppInteger(timeslice) );
    }
    else
    {
        //else we update the map
        Node * NP = getNode(nodeid);
        NodeIdToNodeP[nodeid] = NP;

        NP->setNodeProperty(ts, BppInteger(timeslice));
    }*/
    Node * NP = getNode(nodeid);

    NP->setNodeProperty(ts, BppInteger(timeslice));
}

void ReconciledTree::setNodeTimeSlice(Node * NP, int timeslice)
{
    NP->setNodeProperty(ts, BppInteger(timeslice)); 
}


void ReconciledTree::setNodeEvent(int nodeid, int eventid) // assigns an event to the node according to the eventid
{
    /*
    map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
    if(it != NodeIdToNodeP.end())
    {
        it->second->setNodeProperty(ev, BppInteger(eventid) );
    }
    else
    {
        //else we update the map
        Node * NP = getNode(nodeid);
        NodeIdToNodeP[nodeid] = NP;

        NP->setNodeProperty(ev, BppInteger(eventid));
    }*/
    Node * NP = getNode(nodeid);
    NP->setNodeProperty(ev, BppInteger(eventid));
}

void ReconciledTree::setNodeEvent(Node * NP, int eventid)
{
    NP->setNodeProperty(ev, BppInteger(eventid));   
}


void ReconciledTree::setNodeCladeNum(int nodeid, int cladenum) // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{
    if(hasNodeProperty(nodeid, clnum))//property already set
    {
        if(getNodeCladeNum(nodeid) == cladenum)//already the correct property -> don't change anything
            return;
        else
            resetNodeCladeNum(nodeid); // reset the property and update the map cladeid -> nodeid
    }

    addNodeIdToCladeIdMap(cladenum, nodeid);//updating the map cladeid -> nodeid
    setNodeProperty(nodeid, clnum, BppInteger(cladenum));
}

void ReconciledTree::setNodeCladeNum(Node * NP, int cladenum) // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{
    int nodeid = NP->getId();

    if(NP->hasNodeProperty(clnum))
    {
        if(getNodeCladeNum(nodeid) == cladenum)
            return;
        else
            resetNodeCladeNum(nodeid);
    }


    addNodeIdToCladeIdMap(cladenum, nodeid);//updating the map cladeid -> nodeid
    NP->setNodeProperty(clnum, BppInteger(cladenum));

}



//unset the clade id of the node. Always use this function because it also updates the CladeIdToNodeIds map
void ReconciledTree::resetNodeCladeNum(int nodeid)
{
    int CladeId = getNodeCladeNum(nodeid);

    getNode(nodeid)->deleteNodeProperty(clnum);

    ///updating the map
    if(CladeIdToNodeIds.count(CladeId) != 0)//key is not absent from the map
    {
        vector<int>::iterator it;
        for (it = CladeIdToNodeIds[CladeId].begin() ; it != CladeIdToNodeIds[CladeId].end(); ++it)
            if(*it == nodeid)
                break;

        if(it != CladeIdToNodeIds[CladeId].end())
            CladeIdToNodeIds[CladeId].erase(it);//erasing
                

        if(CladeIdToNodeIds[CladeId].size() == 0)
            CladeIdToNodeIds.erase(CladeId);
    }
    //the else case should never happen if the cladenum was set properly
}


/*
Prints node information

Takes:
 - NodeId (int): id of a node in the tree

 */
void ReconciledTree::printNode(int nodeid)
{
    //TODO pointer version

    cout << "Node : " << nodeid ;

    if(hasNodeName(nodeid))
        cout << " " <<getNodeName(nodeid);

    if(hasFather(nodeid))
        cout << " ; father : " << getFatherId(nodeid);
    else
        cout << " ; root ";

    if(hasNodeProperty(nodeid, spe))
        cout << " ; Species " << getNodeSpecies(nodeid);
    else
        cout << " ; no Species";

    if(hasNodeProperty(nodeid, ev))
        cout << " ; Event " << getNodeEvent(nodeid);
    else
        cout << " ; no Event";

    if(hasNodeProperty(nodeid,clnum))
        cout << " ; CladeNum " << getNodeCladeNum(nodeid);
    else
        cout << " ; no CladeNum";

    if(hasNodeProperty(nodeid,uts))
    {
        int TS = getNodeTimeSlice(nodeid);
        cout << " ; Time Slice " << TS;

    }
    else
        cout << " ; no TS";

//  if(!isLeaf(nodeid))
//  {
//      cout << " ; sons :" ;
//      vector <int> SonIds = getSonsId(nodeid);
//      for(unsigned i = 0; i < SonIds.size(); i++)
//          cout << SonIds[i] << " ";
//  }


    cout << endl;
}

void ReconciledTree::printMe()//prints all nodes
{
    //TODO pointer version
    vector <int> NodesIds = getNodesId();
    for (unsigned j=0; j<NodesIds.size() ; j++)
        printNode(NodesIds[j]);

}

string ReconciledTree::NewickString(bool hideLosses )
{
    return NodeString(getRootId(), hideLosses);
}

/*
    Removes all nodes whose event is No Event: 6
*/
void ReconciledTree::removeNoEventNodes()
{
    //TODO pointer version
    vector<int> NodesIds  = getNodesId();
    int NodeId;

    for (unsigned j=0; j<NodesIds.size() ; j++)
    {
        NodeId= NodesIds[j];

        if(getNodeEvent(NodeId) == N)//NoEvent node
        {
            resetNodeCladeNum(NodeId);//forgetting that the node has this cladenum

            Node * NoEventNode = getNode(NodeId);

            Node * fatherNode = NoEventNode->getFather();

            Node * sonNode = NoEventNode->getSon(0);//we want the 1st son (there should be 1 and only 1 son to a NoEvent node)

            fatherNode->removeSon(NoEventNode);
            fatherNode->addSon(sonNode);// also tells the son that its father is the fatherNode

            delete NoEventNode;
        }
    }

}


/*
Removes null event nodes and set to not time sliced
*/
void ReconciledTree::setToNonTimeSliced()
{
    if( getTimeSliceStatus() == 0)//first check if the time slices exists
        return

    //1. delete NoEvent nodes
    removeNoEventNodes();

    //2. setting the TimeSliceStatus to notify of BTS
    TimeSliceStatus = 0;

}

void ReconciledTree::SubdivideTree()
{
    //expects that node already have an exact set time slice
    if( getTimeSliceStatus() == 0)//first check if the time slices exists
        throw Exception("ReconciledTree::SubdivideTree : No time slices to draw from!");


    vector<int> NodesIds  = getNodesId(); //TODO ; use pointer rather than ids -> faster.
    int NodeId;

    for (unsigned j=0; j<NodesIds.size() ; j++)
    {
        NodeId= NodesIds[j];


        Node * currentNode = getNode(NodeId);
    
        int currentTS;

      //cout << NodeId << " " << currentNode->hasNodeProperty(uts) << "->" ;

        if(! currentNode->hasNodeProperty(ts) ) // special loss cases
        {
            Node * fatherNode = currentNode->getFather();

            currentTS = getNodeTimeSlice(fatherNode) - 1;
            if(currentTS < 0)
                currentTS = 0;

            //correcting the tree while we're at it
            setNodeTimeSlice(currentNode, currentTS);//set them at the current time slice
        }
        else
            currentTS = getNodeTimeSlice(currentNode);

        if(!currentNode->hasFather())//this is the root node -> do nothing else
            continue;



        Node * fatherNode = currentNode->getFather();

        int fatherTS = getNodeTimeSlice(fatherNode);

        int nbMissingNodes = fatherTS - currentTS - 1;


        if(nbMissingNodes > 0 )
        {
            string sp = getNodeSpecies(currentNode);
            if( isRec( getNodeEvent(currentNode) ) )
                sp = -1; //<- because the Null parent of a reception is in the dead

            int evt = N; //Null event
            int cladenum = getNodeCladeNum(NodeId);

            Node * newFatherNode = fatherNode;



            for(int i = 0; i< nbMissingNodes;i++)//creating the missing nodes 
            {


                int newId = getNextId();

                Node * newNode  = new Node(newId);

                //branching to the father
                newFatherNode->addSon(newNode);

                //setting the new node properties
                setNodeSpecies(newNode,sp);
                setNodeEvent(newNode,evt);
                setNodeCladeNum(newNode,cladenum);

                setNodeTimeSlice(newNode, fatherTS - i - 1);


                //incrementing
                newFatherNode = newNode;
            }

            fatherNode->removeSon(currentNode);
            newFatherNode->addSon(currentNode);//branching current node to the Null node before it
        }

    }

    //set time slice status to notify of TS
    TimeSliceStatus = 1;    
}


int ReconciledTree::getNumberOfDuplication()
{
    return countAnyEvent(D);
}

int ReconciledTree::getNumberOfLoss()
{
    return countAnyEvent(L);
}

int ReconciledTree::getNumberOfTransfer()
{
    return countAnyEvent(R);
}


/*
Checks if the node have the same species

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node * ): a node with at least the species property set

Returns:
    (bool) : true if the two nodes have the same species; false otherwise

*/
bool ReconciledTree::haveSameSpecies(Node  * n1, Node * n2)
{
    //basic checks
    if((!n1->hasNodeProperty(spe)) || (!n2->hasNodeProperty(spe)))
        throw Exception("ReconciledTree::haveSameSpecies : given nodes don't have a species (S property)");

    string n1s = getNodeSpecies(n1);
    string n2s = getNodeSpecies(n2); 

    if(n1s != n2s)
        return false;

    return true;
}

/*
Checks if the node have compatible time slices

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node * ): a node with at least the species property set

Returns:
    (bool) : true if the two nodes have compatible time slices; false otherwise

*/
bool ReconciledTree::areTSCompatible(Node  * n1, Node * n2)
{

    if(TimeSliceStatus == 0) // no TS -> same species is sufficient
        return true;

    //if(TimeSliceStatus == 1)
    //{
    // else, we can presume ts is 1
    int n1ts = getNodeTimeSlice(n1);
    int n2ts = getNodeTimeSlice(n2);
    if(n1ts == n2ts)
        return true;
    return false;


}


/*
Checks if the node are compatible (ie. same species and comparable timeslaice).
Presumes that both nodes have the same timeslice status as the tree

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node *): a node with at least the species property set

Returns:
    (bool) : true if the two nodes are compatible; false otherwise

*/
bool ReconciledTree::areCompatible(Node  * n1, Node * n2)
{
    //first species
    if(!haveSameSpecies( n1,  n2))
        return false;
    return areTSCompatible(n1, n2);
}

/*
Checks if the node are compatible (ie. same species and comparable timeslice).
Presumes that both nodes have the same timeslice status as the tree

Takes:
 - id1 (int): a node id with at least the species property set
 - id2 (int): a node id with at least the species property set

Returns:
    (bool) : true if the two nodes are compatible; false otherwise

*/
bool ReconciledTree::areCompatible(int id1, int id2)
{
    Node * n1 = getNode(id1);
    Node * n2 = getNode(id2);
    return areCompatible(n1, n2);
}


/*
Checks if AncId is the ancestor of SonId

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
    (bool): true if AncId is an ancestor of SonId
*/
bool ReconciledTree::isAncestor(int AncId, int SonId)
{
    //solving the easiest cases
    if(AncId  == SonId)
        return false;

    int RootId = getRootId();

    if(AncId == RootId) // the potential ancestor is the root so it is ancestor to SonId
        return true;

    if(SonId == RootId) // The descendant is the root -> impossible
        return false;

    Node * SonNode = getNode(SonId);

    Node * FatherNode = SonNode->getFather();

    int FatherId = FatherNode->getId();

    bool ok= false; //storing result

    while(FatherId != RootId)
    {
        if(FatherId == AncId) // found the Ancestor
        {
            ok = true;
            break;
        }
        Node * SonNode = FatherNode;
        FatherNode = SonNode->getFather();      
        FatherId = FatherNode->getId();
    }
    return ok;
}


/*
Get the path from ancestor to son. (If AncId is not an ancestor of SonId, will return the path to the root)

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
    (vector <Node *>): path including start and stop
*/
vector <Node *> ReconciledTree::pathFromAncestorToSon(int AncId, int SonId)
{


    vector <Node *> path;
    path.push_back(getNode(SonId));

    int RootId = getRootId();
    int currentId = SonId;

    while((currentId != RootId)&&(currentId != AncId))
    {
        path.insert(path.begin(),path.front()->getFather());
        currentId = path.front()->getId();
    }

    return path;
}

/*
Get the path from ancestor to son. (If AncId is not an ancestor of SonId, will return the path to the root)

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
    (vector <int >): path including start and stop
*/
vector <int> ReconciledTree::idPathFromAncestorToSon(int AncId, int SonId)
{
    vector <Node *> path = pathFromAncestorToSon(AncId, SonId);

    vector <int> idpath;

    for(unsigned i = 0; i < path.size(); i++)
        idpath.push_back(path[i]->getId());

    return idpath;
}


int ReconciledTree::getIdWithName(const string &name)
{
    //go through the whole tree. Does not throw error if a node has no name
    int IdToReturn = -1;

    vector <Node *> NodeToDo;
    Node * currentNode;

    NodeToDo.push_back(getRootNode());

    while(NodeToDo.size() > 0)
    {
        currentNode = NodeToDo[0];
        if(currentNode->hasName())
        {

            if(currentNode->getName() == name)
            {
                IdToReturn = currentNode->getId();
                break;
            }
        }

        for(unsigned i = 0 ; i < currentNode->getNumberOfSons() ; i++)
        {
            NodeToDo.push_back(currentNode->getSon(i));
        }
        NodeToDo.erase(NodeToDo.begin());
    }

    if(IdToReturn == -1)
        throw Exception("ReconciledTree::getIdWithName: no nodes with name \"" + name  + "\"");

    return IdToReturn;

}


bool ReconciledTree::isRealLeaf(int id)
{
    if(getNodeEvent(id) != C)
        return false;
    return true;
}

bool ReconciledTree::isExtant(int evtcode)
{
    if(evtcode != C)
        return false;
    return true;
}

bool ReconciledTree::isSpeciation(int evtcode)
{
    if(evtcode != S)
        return false;
    return true;
}

bool ReconciledTree::isLoss(int evtcode)
{
    if(evtcode != L)
        return false;
    return true;
}

bool ReconciledTree::isDup(int evtcode)
{
    if(evtcode != D)
        return false;
    return true;
}

bool ReconciledTree::isSout(int evtcode)
{
    if(evtcode != Sout)
        return false;
    return true;
}

bool ReconciledTree::isRec(int evtcode)
{
    if(evtcode != R)
        return false;
    return true;
}

bool ReconciledTree::isNull(int evtcode)
{
    if(evtcode != N)
        return false;
    return true;
}

bool ReconciledTree::isBout(int evtcode)
{
    if(evtcode != Bout)
        return false;
    return true;
}


ReconciledTree * ReconciledTree::cloneSubtree(int newRootId)
{
    Node * newRoot = TreeTemplateTools::cloneSubtree<Node>(*this, newRootId);
    return new ReconciledTree(newRoot, getTimeSliceStatus());
}

vector <string> ReconciledTree::getRealLeavesNames()
{
    vector <string>  names;

    vector <int> LeavesId = getLeavesId();

    for(size_t i = 0; i < LeavesId.size(); i++)
    {
        if( isRealLeaf(LeavesId[i]) )
            names.push_back(getNodeName(LeavesId[i]));

    }
    return names;
}
