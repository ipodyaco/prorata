
#include "sicReaderHelper.h"


// Variable used to trace the Chromatogram tag while reading
// the file serially.
int iEntered = 0;

// Variables to hold header info
string sMSfile;
int iFirstId;
int iLastId;
string sProgram;
string sVersion;

// Variable to hold Choramatogram identifier.
string sChroId;
int iChroId;

// Variable to hold the actuall contents of the Chromatogram
// tag.
string sChroData;


int cstream(BYTE *buf, int cBytes, int *cBytesActual, void *inputData)
{
	*cBytesActual = fread(buf, 1, cBytes, (FILE*)inputData);
	return (*cBytesActual < cBytes);
}

int StartElementForHeader(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName, LPXMLVECTOR atts)
{
	if ( strcmp( (char*) qName, "EXTRACTED_ION_CHROMATOGRAMS" ) == 0 )
	{
		if (atts->length) 
		{
			int i;
			LPXMLRUNTIMEATT att;

			for (i=0; i<atts->length; i++) {
				att = (LPXMLRUNTIMEATT) XMLVector_Get(atts, i);

				string sQName = (char* )att->qname;
				//cout << "qname = " << sQName << endl;

				if ( sQName ==  "MSfile" )
				{
					sMSfile = (char *)att->value;
					continue;
				}

				if ( sQName ==  "first_identifier" )
				{
					iFirstId = atoi( (char *)att->value );
					continue;
				}

				if ( sQName ==  "last_identifier" )
				{
					iLastId = atoi( (char *)att->value );
					continue;
				}

				if ( sQName ==  "program" )
				{
					sProgram = (char *)att->value;
					continue;
				}

				if ( sQName ==  "version" )
				{
					sVersion = (char *)att->value;
					continue;
				}
			}
		}
	}
	return XML_ABORT;
}

int StartElement(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName, LPXMLVECTOR atts)
{
  	int *depthPtr = (int *)UserData;
	
	if ( strcmp( (char*) qName, "CHROMATOGRAM" ) == 0 )
	{
		if (atts->length) 
		{
			LPXMLRUNTIMEATT att;

			att = (LPXMLRUNTIMEATT) XMLVector_Get(atts, 0);

			if( (strcmp( (char *)att->qname, "identifier" ) == 0 &&
					( atoi( (char*) att->value ) == iChroId ) ) )
			{
				iEntered = 1;
			}
		}
	}

	if ( iEntered )
	{
		int i;

		for (i = 0; i < *depthPtr; i++)
		{
			sChroData += "\t";
		}

		sChroData = sChroData + "<" + (char *) qName ;

		if (atts->length) 
		{
			LPXMLRUNTIMEATT att;
			
			for (i=0; i<atts->length; i++) 
			{
				att = (LPXMLRUNTIMEATT) XMLVector_Get(atts, i);
				sChroData = sChroData + " " + (char *) att->qname +
					"=\"" + (char *) att->value + "\"";
			}
		}
		sChroData = sChroData + ">\n";
		*depthPtr += 1;
	}

	return 0;

}

int Characters(void *UserData, const XMLCH *chars, int cbChars)
{
	int i = 0;
  	int *depthPtr = (int *)UserData;

	if ( iEntered )
	{

		for (int j = 0; j < (*depthPtr); j++)
		{
			sChroData += "\t";
		}

		while ( i < cbChars )
		{
			sChroData += chars[i];
			i++;
		}
		sChroData = sChroData + "\n";
	}
	return 0;	
}

int EndElement(void *UserData, const XMLCH *uri, const XMLCH *localName, const XMLCH *qName)
{	
	int *depthPtr = (int*)UserData;

	if ( iEntered )
	{
		*depthPtr -= 1;

		for (int j = 0; j < (*depthPtr); j++)
		{
			sChroData += "\t";
		}

		sChroData = sChroData + "</" + (char*) qName + ">";
		sChroData = sChroData + "\n";

		if ( strcmp( (char*) qName, "CHROMATOGRAM" ) == 0 )
		{
			iEntered = 0;
		}
	}

	return 0;
}

void ErrorHandler(LPXMLPARSER parser)
{
}

int processHeader( const char *cpFileName )
{

	LPXMLPARSER parser;		
	FILE *fXMLFile; 
	int iReturns;
	int depth = 0;

	sMSfile = "";
	iFirstId = 0;
	iLastId = 0;
	sProgram = "";
	sVersion = "";

	if ( ( fXMLFile = fopen( cpFileName, "r" ) ) == NULL )
	{
		fprintf( stderr, "Error: Cannot open %s\n", cpFileName );
		return 1;
	}

	if (!XMLParser_Create(&parser)) {
		fprintf( stderr, "Error: creating parser!\n" );
		return 1;
	}

	parser->errorHandler = ErrorHandler;
	parser->startElementHandler = StartElementForHeader;
	parser->UserData = &depth;

	iReturns =  XMLParser_Parse( parser, cstream, fXMLFile, 0 );

	if ( iReturns == XML_ABORT )
	{
		// Ignoring abort errors for now.
	}

	/*
	cout << "MSfile = " << sMSfile << endl;
	cout << "sFirstId = " << iFirstId << endl;
	cout << "sLastId = " << iLastId << endl;
	cout << "sProgram = " << sProgram << endl;
	cout << "sVersion = " << sVersion << endl;
	*/

	XMLParser_Free(parser);
	fclose( fXMLFile );

	return 0;

}

string getChro( const char *cpFileName, int iChromaId )
{

	LPXMLPARSER parser;		
	FILE *fXMLFile; 
	int iReturns;
	int depth = 0;

	if ( ( fXMLFile = fopen( cpFileName, "r" ) ) == NULL )
	{
		fprintf( stderr, "Error: Cannot open %s\n", cpFileName );
		return "";
	}

	if ( iChromaId <= 0 )
	{
		fprintf( stderr, "Error: Chromatogram identifier %d "
				"is not valid\n", iChroId );
		return "";
	}

	sChroData = "";

	/*
	if ( !cpChroId )
	{
		fprintf( stderr, "Error: Chromatogram identifier %s "
				"is not valid\n", cpChroId );
		return "";
	}
	*/

	iChroId = iChromaId;

	if (!XMLParser_Create(&parser)) {
		fprintf( stderr, "Error: creating parser!\n" );
		return "";
	}

	parser->errorHandler = ErrorHandler;
	parser->startElementHandler = StartElement;
	parser->charactersHandler = Characters;
	parser->endElementHandler = EndElement;	
	parser->UserData = &depth;


	iReturns =  XMLParser_Parse( parser, cstream, fXMLFile, 0 );

	if ( iReturns == XML_ABORT )
	{
		// Ignoring abort errors for now.
	}

	//cout << sChroData << endl;

	sChroId = "";

	XMLParser_Free(parser);
	fclose( fXMLFile );

	return sChroData;

}

// Utility (Accessor) functions.
string getMSFile()
{	
	return sMSfile;		
}

int getFirstId()
{	
	return iFirstId;	
}

int getLastId()
{
	return iLastId;		
}

string getProgram()
{	
	return sProgram;	
}

string getVersion()
{	
	return sVersion;	
}


#if 0
int main(int argc, char* argv[])
{	
	//LPXMLPARSER parser;		
	//int depth = 0;
	//FILE *fXMLFile; 
	char cpFileName[] = "sample.xml";
	//int iReturns;

	/*
	getHeaders( cpFileName );
	*/

	getChromatogram( cpFileName, 1 );
#if 0

	if (!XMLParser_Create(&parser)) {
		printf("Error creating parser!\n");
		return 1;
	}

	parser->errorHandler = ErrorHandler;
	//parser->startElementHandler = StartElement;
	parser->startElementHandler = StartElementForHeader;
	parser->charactersHandler = Characters;
	parser->endElementHandler = EndElement;	
	parser->UserData = &depth;


	if ( ( fXMLFile = fopen( cpFileName, "r" ) ) == NULL )
	{
		fprintf( stderr, "Error: Cannot open %s\n", cpFileName );
	}

	//if (!XMLParser_Parse(parser, cstream, stdin, 0))
	
	iReturns =  XMLParser_Parse( parser, cstream, fXMLFile, 0 );

	if ( iReturns == XML_ABORT )
	{
		// Ignoring abort errors for now.
	}
	/*
	if (!XMLParser_Parse(parser, cstream, fXMLFile, 0))
		printf("Error: %s\nLine: %d Col: %d\n", 
			parser->ErrorString, parser->ErrorLine, parser->ErrorColumn);
	*/
			
	XMLParser_Free(parser);

	//printf( "Choma Count = %d\n", iChromaCount );
	
#endif
	return 0;
}
#endif

