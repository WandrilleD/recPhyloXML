#ifndef XMLUTILS_H_
#define XMLUTILS_H_

/*

This file contains XML util functions

Created the: 09-02-2015
by: Wandrille Duchemin

Last modified the: 18-09-2017
by: Wandrille Duchemin

*/
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

using namespace std;

bool goToNextTag(ifstream& fileIN , string tag);

string goToNextOf(ifstream& fileIN , vector <string> &tags);

string getLineBaliseName(string line);
string InterpretLineXML(string line, map <string, string> * properties, string * value);



#endif