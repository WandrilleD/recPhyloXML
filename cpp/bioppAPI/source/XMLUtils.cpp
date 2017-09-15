/*

This file contains XML util functions

Created the: 09-02-2015
by: Wandrille Duchemin

Last modified the: 09-02-2016
by: Wandrille Duchemin

*/

#include "XMLUtils.h"




/*
Takes:
 - line (string): a line out of a phyloxml file

Returns:
 (string): name of the balise in the name. "" if there is no balise 
*/
string getLineBaliseName(string line)
{
    size_t baliseBegin = line.find("<");

    bool Closing = false;

    string baliseName;

    if(baliseBegin == string::npos) // character not found -> ignore
        return "";

    baliseBegin++;
    
    if(line[baliseBegin] == '/')
    {
        Closing = true;
        baliseBegin++;
    }

    size_t EndBaliseName = string::npos;
    size_t spacePos = line.find(" ",baliseBegin);
    if(spacePos != string::npos)
        EndBaliseName = spacePos;
    size_t greaterPos = line.find(">",baliseBegin);
    if((greaterPos != string::npos) && (( greaterPos < EndBaliseName) || (EndBaliseName == string::npos)))
        EndBaliseName = greaterPos;
    if(EndBaliseName == string::npos)
        return "";
    else
        baliseName = line.substr(baliseBegin,EndBaliseName - baliseBegin);

    return baliseName;
}

/*
 - line (string ) : line in an XLM file
 - properties (map <string, string> * ) : pointer to a map where keys will be attribute name and associated values will be the value of these attribute
 - value (string * ) : eventual value between the tags

 Returns:
    (string):

*/
string InterpretLineXML(string line, map <string, string> * properties, string * value)
{
    // 1. find the name of the balise

    size_t baliseBegin = line.find("<");

    string baliseName;

    if(baliseBegin == string::npos) // no balise
        return "";

    baliseBegin++;

    size_t EndBaliseName = string::npos;
    size_t spacePos = line.find(" ",baliseBegin);
    if(spacePos != string::npos)
        EndBaliseName = spacePos;
    size_t greaterPos = line.find(">",baliseBegin);
    if((greaterPos != string::npos) && (( greaterPos < EndBaliseName) || (EndBaliseName == string::npos)))
        EndBaliseName = greaterPos;


    if(EndBaliseName == string::npos)
        return "";// no balise
    else
        baliseName = line.substr(baliseBegin,EndBaliseName - baliseBegin);  

    //balise name found
    //2. look for properties
    size_t currentIndex = EndBaliseName;

    currentIndex++;

    while(currentIndex < greaterPos)
    {
        //reading a property
        string PropertyName = "";
        string PropertyValue = "";

        while(line[currentIndex] != '=')
        {
            PropertyName.push_back(line[currentIndex]);
            currentIndex++;
        }
        currentIndex++;
        if(line[currentIndex] == '"') // if first character of the value is a ", then we go to the next one
        {
            currentIndex++;
            while(line[currentIndex] != '"') // reading the value
            {
                PropertyValue.push_back(line[currentIndex]);
                currentIndex++;
            }
            currentIndex++;     
        }
        else // go to the next > or  ' '
        {
            while((line[currentIndex] != '>') && (line[currentIndex] != ' ') ) // reading the value
            {
                PropertyValue.push_back(line[currentIndex]);
                currentIndex++;
            }
        }
        //PropertyValue set, we are either on a ' ' or a '>'
        currentIndex++; // going over that character for next round

        (*properties)[PropertyName] = PropertyValue;

        //cout << "added property " <<  PropertyName << " -> " << PropertyValue << endl;
    }
    //we are after the > character

    size_t lowerpos = line.find("<",currentIndex);

    if(( lowerpos == string::npos)||( lowerpos == currentIndex)) //no value 
        return baliseName;

    //extracting value
    *value = *value + line.substr(currentIndex, lowerpos - currentIndex);

//  cout << line.substr(currentIndex, lowerpos - currentIndex) << endl;
//  cout << *value << endl;

    return baliseName;
}