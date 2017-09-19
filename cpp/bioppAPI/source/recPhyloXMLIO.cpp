/*

This file contains helper functions to read / write recPhyloXML files

Created the: 18-09-2017
by: Wandrille Duchemin

Last modified the: 18-09-2017
by: Wandrille Duchemin

*/

#include "recPhyloXMLIO.h"


/*
read the reconcilied tree in a recphyloXML file and create a GeneFamily for each of them

Takes:
    - vector <ReconciledTree>  &recTreeList : list of reconciled tree where the read reconciled tree will be appended
    - string fileName : name of the recPhyloXML file
    - TreeTemplate<Node> * Stree : either already read species tree, or NULL pointer to receive the species tree that will be read 
    - bool replaceSpTree : whether a species tree will be read in this file or not.
    - unsigned verboseLevel: 0 -> quiet
                             1 -> just little messages whenever trees are read
                             2 -> clade by clade detailed info on progression

*/
void ReadRecPhyLoXMLFile(  vector <ReconciledTree> &recTreeList, string fileName , TreeTemplate<Node> &Stree, bool replaceSpTree, unsigned verboseLevel)
{
    ifstream fileStream(fileName.c_str());
    if( !fileStream.is_open() ) 
    {
        throw Exception("ReadRecPhyLoXMLFile : could not open reconciled tree file : " + fileName);
        exit(1);
    }

    vector <string> TagToSearch;
    TagToSearch.push_back(RECGENETREETAG);


    if(replaceSpTree)
        TagToSearch.push_back(SPTREETAG);
    
    string FoundTag = goToNextOf(fileStream,TagToSearch);


    while( FoundTag != "" )  //goes to the next RecGeneTreeTag -> next reconciled gene tree
    {

        if( FoundTag == SPTREETAG )  //goes to the next SPTREETAG -> SP TREE
        {
            if( goToNextTag(fileStream,CLADETAG) ) // go to clade -> the root of the reconciled gene tree
            {
                //if(Stree != NULL) // if there is a species tree already, we want to delete it properly before we replace it
                //    delete Stree; // note that, if possible, one should avoid this as it creates overhead computations

               if(verboseLevel>0)
                    cout << "Creating species tree." << endl;


                readSpTreeFromPhyloXMLFile(Stree, fileStream, verboseLevel);

            }
        }
        else if(FoundTag == RECGENETREETAG )
        {
            if( goToNextTag(fileStream,CLADETAG) ) // go to clade -> the root of the reconciled gene tree
            {
                //if(verboseLevel>0)
                //    cout << "Creating reconciled tree." << endl;
                TreeTemplate<Node> * TMP = &Stree;
                ReconciledTree Rtree(fileStream, TMP, verboseLevel>1);
                recTreeList.push_back( Rtree );
            }
            else
            { // no clade ? pb but ignore...
                if(verboseLevel>0)
                    cout << "found a "<< RECGENETREETAG << " without any " << CLADETAG << " in the file " << fileName << " -> ignoring that tree."<< endl;
            }
        }

        FoundTag = goToNextOf(fileStream,TagToSearch);
    }

}

/*

Takes:
    - TreeTemplate<Node> * Stree : pointer to fill
    - ifstream& fileIN : input stream. I presume that the stream has yielded a <clade> line forming the root of the tree
    - unsigned verboseLevel: 0 -> quiet
                             1 -> just a little message
                             2 -> clade by clade info on progression


*/
void readSpTreeFromPhyloXMLFile( TreeTemplate<Node> & Stree, ifstream& fileIN, unsigned verboseLevel)
{
     
    int Id = 0;

    Node * newnode  = new Node(Id);

    Id++;

    Stree.setRootNode(newnode);//setting as the new root


    addXMLCladeToSpTree( fileIN,newnode, Id, (verboseLevel>1) ); 

    return;
}


/*
read and add a clade from a recPhyloXML stream
to a simple tree that is used to represent the species tree

!!! This function is recursive

Takes:
    - fileIN (ifstream): streaam to a recPhyloXML file
    - Node * node : node to put the clade IN
    - int &nextId : gives the id to give to the next node (better than the getNexId of bio++ trees which go over the whole tree to determine that number.)
    - VERBOSE (bool) (default = false)
*/
void addXMLCladeToSpTree(ifstream& fileIN, Node * node, int &nextId, bool VERBOSE) 
{
    /// I presume that the "<clade>" opening line of the current clade has already been read.


    if(VERBOSE)
        cout << "parsing Node " << node->getId() << endl;

    //2. parsing of the clade
    string line;
    getline( fileIN, line );

    map <string, string> * properties = new map <string,string>;
    string * value = new string(""); 
    string Tname = InterpretLineXML(line,properties,value);

    //cout << "value:" << *value << endl;

    while(Tname.compare("/clade") != 0) //--> line signing the end of the clade
    {
        if(Tname.compare("name") == 0) // found name
        {
            if(VERBOSE)
                cout << "found name : " << *value <<  endl;
            node->setName( *value);

        }
        else if(Tname.compare("clade") == 0) // found a child --> recursion
        {
            //create child clade

            //creating the new node
            Node * newNode  = new Node(nextId);

            nextId++;

            //branching to the existing structure
            node->addSon(newNode);

            //recursing
            addXMLCladeToSpTree(fileIN, newNode, nextId,  VERBOSE );


        }

        //reading new line
        *value = "";
        delete properties;
        getline( fileIN, line );
        properties = new map <string,string>;
        Tname = InterpretLineXML(line,properties,value);
        
    }

}








/*
Takes:
    - OUT (ofstream) : an opened stream
    - indent_level (int): number of double-spaces at the beginning of a line
*/
void beginLine(ostream& OUT, int indent_level)
{
    for(unsigned i = 0; i< indent_level ; i++)
        OUT << "  ";
}
    
/*
Writes the event of the reconciled node

Takes:
    - OUT (ofstream) : an opened stream
    - Rtree (ReconciledTree *) : a pointer to a reconciled tree
    - node (Node * ): pointer to the current node
    - indent_level (int): number of double-spaces at the beginning of a line
    - hasLoss (bool): true if one of nodeids child is a loss; false otherwise
*/
void WritePhyloXMLRecEvent(ostream& OUT, ReconciledTree * Rtree, Node * node, int indent_level, bool hasLoss)
{
    int event = Rtree->getNodeEvent(node);

    if(Rtree->isNull(event)) // if the event is a Null event it won't be written
        return;

   //cout << nodeid << " "<< indent_level << endl;

    beginLine(OUT,indent_level);
    OUT << "<";

    string eventString = "";
    string SpeciesString = "speciesLocation=";

    if(Rtree->isExtant(event))
    {
        eventString = "leaf";
    }
    else if(Rtree->isSpeciation(event))
    {
        eventString = "speciation";
    }
    else if(Rtree->isLoss(event))
    {
        eventString = "loss";
    }
    else if(Rtree->isDup(event))
    {
        eventString = "duplication";
    }
    else if(Rtree->isSout(event))
    {
        eventString = "speciationOut";
    }
    else if(Rtree->isRec(event))
    {
        eventString = "transferBack";
        SpeciesString = "destinationSpecies=";      
    }
    else if(Rtree->isBout(event))
    {
        eventString = "bifurcationOut";
    }
    else
    {
        eventString = "unknownEvent";
    }

    OUT << eventString;

    if(hasLoss)
        OUT << "Loss";

    if(!Rtree->isBout(event))
    {
        OUT << " ";
        OUT << SpeciesString;
        OUT << "\"";
        OUT << Rtree->getNodeSpecies(node);
        OUT << "\"";
    }

    if(Rtree->isExtant(event))
        if(node->hasName())
            OUT << " geneName=\"" << node->getName() << "\"" ;

    //... time slice if any.
    if((Rtree->getTimeSliceStatus() != 0) && (eventString != "leaf"))
    {
        OUT << " timeSlice=\"";
        OUT <<Rtree->getNodeTimeSlice(node);
        OUT << "\"";
    }

    OUT << "></";
    OUT << eventString;
    if(hasLoss)
        OUT << "Loss";

    OUT << ">\n";

    return;
}


/*
Writes the reconciled node along with its children (recursively)

Takes:
    - OUT (ofstream) : an opened stream
    - Rtree (ReconciledTree *) : a pointer to a reconciled tree
    - node (Node * ): pointer to the current node
    - indent_level (int): number of double-spaces at the beginning of a line
*/
void WritePhyloXMLRecTreeAux(ostream& OUT, ReconciledTree * Rtree, Node * node, int indent_level)
{
    
    int event = Rtree->getNodeEvent(node);
    if(Rtree->isNull(event))// ignoring any Null event.
    {
        vector < Node * > children = node->getSons();
        WritePhyloXMLRecTreeAux(OUT, Rtree, children[0], indent_level);//recursion on the first (and only) son of the null node
    }
    else
    {

        beginLine(OUT,indent_level);
        OUT << "<clade>\n";
        indent_level++;
    
        //writing the different informations:
    
        //1. name:
        beginLine(OUT,indent_level);
        OUT << "<name>";
        if(node->hasName())
        {
            OUT << node->getName();
        }
        else
        {
            OUT << node->getId();
        }
        OUT << "</name>\n";
    
        //2. the reconciliation event(s)
        beginLine(OUT,indent_level);
        OUT << "<eventsRec>\n";
        indent_level++;
    
        Node * currentNode = node;
        Node * nextNode = currentNode;
        bool hasLoss = false;
        bool donext = true;
    
        while( donext )
        {
            currentNode = nextNode; //iteration
            hasLoss = false;
    
            vector < Node * > children = currentNode->getSons();
    
            if(children.size() == 0)
                donext = false;
    
            for(unsigned i = 0 ; i < children.size(); i++)
            {
                int event = Rtree->getNodeEvent(children[i]);
                if(Rtree->isLoss(event))
                {
                    hasLoss = true;
                }
                else
                {
                    if(nextNode  == currentNode)//this is the first non-loss child of current id -> we accept it as nextid
                        nextNode = children[i];
                    else //this is not the first non-loss child of current id -> we get out of the loop
                    {
                        donext = false;
                        break;
                    }
                }
            }
            //we write the event of the current node
            WritePhyloXMLRecEvent(OUT,  Rtree, currentNode, indent_level,  hasLoss);
        }
    
        indent_level--;
        beginLine(OUT,indent_level);
        OUT << "</eventsRec>\n";
    
    
        //3. Sons, if any
        vector < Node * > children = currentNode->getSons();
            
        for(unsigned i = 0 ; i < children.size(); i++)
        {
            WritePhyloXMLRecTreeAux( OUT,  Rtree, children[i], indent_level);
        }
    
        indent_level--;
        beginLine(OUT,indent_level);
        OUT << "</clade>\n";
    }
}


/*
Writes the sp node along with its children (recursively)

Takes:
    - OUT (ofstream) : an opened stream
    - Stree (MySpeciesTree *) : a pointer to a species tree
    - node (Node * ): current node
    - indent_level (int): number of double-spaces at the beginning of a line
*/
void WritePhyloXMLSpeTreeAux(ostream& OUT, TreeTemplate< Node > * Stree, Node * node, int indent_level)
{

    beginLine(OUT,indent_level);
    OUT << "<clade>\n";
    indent_level++;

    double distance = 0;

    Node * current = node;
    vector < Node * > children = node->getSons();
  
    while(children.size() ==1)
    {
        current = children[0];
        //cout  << "->" << currentId;
        children = current->getSons();

        if(current->hasDistanceToFather())
            distance += current->getDistanceToFather();
    }

    //writing the different informations:

    //1. name:
    beginLine(OUT,indent_level);
    OUT << "<name>";
    if( node->hasName() ) //internal nodes actually have names that do not correspond to their RPO but the TS of the speciation
    {
        OUT << node->getName();
    }
    else
    {
        OUT << node->getId();
    }

    OUT << "</name>\n";



    //2. distance to father
    if(distance != 0)
    {
        beginLine(OUT,indent_level);
        OUT << "<branch_length>" << distance << "</branch_length>\n" ;
    }

    for(unsigned i = 0 ; i < children.size(); i++)
    {
        WritePhyloXMLSpeTreeAux( OUT,  Stree, children[i], indent_level);
    }   

    indent_level--;
    beginLine(OUT,indent_level);
    OUT << "</clade>\n";
    

}


/*
Writes the whole reconciled tree in the of stream

Takes:
    - OUT (ofstream) : an opened stream
    - Rtree (ReconciledTree *) : a pointer to a reconciled tree
    - int index (default: -1): index of the reconciled tree to write (or -1 to omit it)
*/
void WritePhyloXMLRecTree(ostream& OUT, ReconciledTree * Rtree, int index)
{
    //OUT << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    //OUT <<"<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\" xsi:schemaLocation=\"../recxml.xsd\">\n";

    int indent_level = 1;
    beginLine(OUT,indent_level);
    OUT << "<recGeneTree>\n";

    indent_level++;

    beginLine(OUT,indent_level);
    OUT << "<phylogeny rooted=\"true\">\n";

    indent_level++;

    if( index != -1)
    {
        beginLine(OUT,indent_level);
        OUT << "<id>"<< index <<"</id>"<<endl;
    }

    WritePhyloXMLRecTreeAux( OUT,  Rtree, Rtree->getRootNode(), indent_level);

    indent_level--;
    beginLine(OUT,indent_level);
    OUT << "</phylogeny>\n";

    indent_level--;
    beginLine(OUT,indent_level);

    OUT << "</recGeneTree>\n";

}


/*
Writes the whole species tree in the of stream

Takes:
    - OUT (ofstream) : an opened stream
    - Stree (TreeTemplate<Node> *) : a pointer to a species tree

*/
void WritePhyloXMLSpeTree(ostream& OUT, TreeTemplate<Node> * Stree)
{
    int indent_level = 1;
    beginLine(OUT,indent_level);
    OUT << "<spTree>\n";

    indent_level++;

    beginLine(OUT,indent_level);
    OUT << "<phylogeny rooted=\"true\">\n";

    indent_level++;

    WritePhyloXMLSpeTreeAux( OUT,  Stree, Stree->getRootNode(), indent_level);

    indent_level--;
    beginLine(OUT,indent_level);
    OUT << "</phylogeny>\n";

    indent_level--;
    beginLine(OUT,indent_level);
    OUT << "</spTree>\n";
}



/*
Takes:
    - string filename : name of the file to open

*/
void WritePhyloXML(ostream& OUT, vector <ReconciledTree> &recTreeList , TreeTemplate<Node> * Stree)
{

    OUT << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" <<endl;
    OUT << "<recPhylo xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\" xmlns=\"http://www.recg.org\" >" << endl;

    bool NoSpeciesTree = false;
    if(Stree == NULL)
        NoSpeciesTree=true;
    else if(Stree->getRootNode() == NULL)
        NoSpeciesTree=true;

    if(!NoSpeciesTree)
    {
        WritePhyloXMLSpeTree( OUT , Stree );
    }


    ReconciledTree * RT;
    for( unsigned i = 0 ; i < recTreeList.size() ; i++ )
    {
        RT = &(recTreeList[i]) ;
        WritePhyloXMLRecTree(OUT , RT , -1);
    }
        


    OUT << "</recPhylo>" << endl;

}

/*
Takes:
    - string filename : name of the file to open

*/
void WritePhyloXML(string fileName, vector <ReconciledTree> &recTreeList , TreeTemplate<Node> * Stree)
{
    ofstream ofs;
    ofs.open(fileName.c_str(),ofstream::out);

    WritePhyloXML(ofs, recTreeList , Stree);
}