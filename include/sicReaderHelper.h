
#include <iostream>
#include <string>

using namespace std;


#include "parsifal.h"


// This callback function is called when the reader encounters a start tag.
// It checks for the chromatogram tag of a particular identifier. If found, it 
// raises a flag which will be used for redirecting the contents of that tag
// in to a variable for storage.
int StartElement(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName, LPXMLVECTOR atts);

// This callback function is called when the reader encounters start tag.
// It checks for the EXTRACTED_ION_CHROMATOGRAMS tag and once found stores all
// the attributes.
int StartElementForHeader(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName, LPXMLVECTOR atts);

// Callback to track all the end tags.
int EndElement(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName);

// Callback to store all the contents of a particular tag. It starts storing
// once a particular chromatogram is found.
int Characters(void *UserData, const XMLCH *chars, int cbChars);

// Dummy function for error handling.
// Future enhancements: Make it handle XML validation.
void ErrorHandler(LPXMLPARSER parser);

// Used for internal bit streaming.
int cstream(BYTE *buf, int cBytes, int *cBytesActual, void *inputData);

// The function has to be executed before calling any accessor functions.
// This will process the file and populate header variables.
int processHeader( const char *cpFileName );

// Given an identifier it returns the chromatogram contents.
// Returns "" if not found.
string getChro( const char *cpFileName, int iChromaId );

// Utility (Accessor) functions.
string getMSFile();

int getFirstId();

int getLastId();

string getProgram();

string getVersion();

