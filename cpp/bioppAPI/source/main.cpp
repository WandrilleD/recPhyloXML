/*
Created by: Wandrille Duchemin

Last modified the: 15-09-2017
by: Wandrille Duchemin

*/



#include "ReconciledTree.h"


void simpleTest()
{
    ReconciledTree RT();
    return; 
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
    
    cout << "You gave me args, but I don't know what to do with them. try without anything." <<endl;
    
    return 0;
}