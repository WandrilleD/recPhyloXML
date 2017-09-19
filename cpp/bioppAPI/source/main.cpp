/*
Created by: Wandrille Duchemin

Last modified the: 18-09-2017
by: Wandrille Duchemin

*/



#include "ReconciledTree.h"
#include "recPhyloXMLIO.h"

void simpleTest()
{
    ReconciledTree RT();
    return; 
}


bool readingTest(string filename , vector <ReconciledTree>  &recTreeList, TreeTemplate<Node> & Stree)
{

    ReadRecPhyLoXMLFile(  recTreeList,  filename , Stree, true, 2);

    return true;
}

/////////////////////////////////////////////////
// Main 
/////////////////////////////////////////////////

int main(int args, char ** argv)
{
    if(args == 1)
    {
        simpleTest();
        cout << "performed a simple test. OK"<< endl;
        exit(0);
    }
    else
    {
        cout << "You gave me args, I will try to read the first one as reconciled trees in the recPhyloXML format." <<endl;
        
        string fileName(argv[1]);

        cout << "attempting to read from " <<  fileName <<endl;

        vector <ReconciledTree>  recTreeList;
        TreeTemplate<Node> * Stree = new TreeTemplate<Node>();
        bool worked = readingTest(fileName , recTreeList, *Stree);


        if(worked)
        {
            cout << "reading done."<< endl;

            bool NoSpeciesTree = false;
            if(Stree == NULL)
                NoSpeciesTree=true;
            else if(Stree->getRootNode() == NULL)
                NoSpeciesTree=true;

            if(NoSpeciesTree)
                cout << "did not ";

            cout << "read a species tree." << endl;

            cout << "read " << recTreeList.size() << " ReconciledTree"<<endl;
        }

        cout << "Now I will try to write the trees I just read in the recPhyloXML format." <<endl;
        ostream* OUT = &cout;
        ofstream fout;

        if(args > 2) {
            cout << "writing into the file given as second argument : " << argv[2] <<endl;
            fout.open(argv[2]);
            OUT = &fout;
        }

        WritePhyloXML(*OUT, recTreeList , Stree);
    }
 

    return 0;
}