/*

This file contains XML util functions

Created the: 09-02-2015
by: Wandrille Duchemin

Last modified the: 18-09-2017
by: Wandrille Duchemin

*/

#include "XMLUtils.h"

/*
returns true if target is in ref
*/
bool isIn(string target, vector <string> &ref)
{
    for(unsigned i = 0 ; i < ref.size() ; i++)
    {
        if(target.compare(ref[i]) == 0)
            return true;
    }
    return false;
}


/*
Takes:
 - ifstream& fileIN : file the data comes from
 - vector<string> & tags : lost of tag to find the next of

Returns:
 (string) : tag that was found
            "" if none wasas found

*/
string goToNextOf(ifstream& fileIN , vector <string> &tags)
{
 if( fileIN.eof() )
     return "";

 string line;
 getline( fileIN, line );

 map <string, string> * properties = new map <string,string>;
 string * value = new string(""); 

 string Tname = InterpretLineXML(line,properties,value);

 while(!isIn(Tname, tags))
 {
    if(fileIN.eof())
         return "";

    getline( fileIN, line );
    Tname = InterpretLineXML(line,properties,value);
 }

 return Tname;
}



/*
Takes:
 - ifstream& fileIN : file the data comes from
 - string tag : tag to find the next of

Returns:
 (bool) : true if the next tag was found
*/
bool goToNextTag(ifstream& fileIN , string tag)
{
 if( fileIN.eof() )
     return false;

 string line;
 getline( fileIN, line );

 map <string, string> * properties = new map <string,string>;
 string * value = new string(""); 

 string Tname = InterpretLineXML(line,properties,value);


 while(Tname.compare(tag) != 0) 
 {
    if(fileIN.eof())
         return false;

    getline( fileIN, line );
    Tname = InterpretLineXML(line,properties,value);
 }

 return true;
}




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