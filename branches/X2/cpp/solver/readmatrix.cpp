

#include "readmatrix.h"

template bool ReadMatrix::GetNextMatrix (vector< vector< int > > &  vvMatrix);
template bool ReadMatrix::GetNextMatrix (vector< vector< double > > &  vvMatrix);

ReadMatrix::ReadMatrix()
{

}


ReadMatrix::~ReadMatrix()
{

}

void ReadMatrix::InitializeParser ( const string & sMatrixFileName )
{
    MatrixFile_stream.open(sMatrixFileName.c_str());
    if (MatrixFile_stream.is_open() == false)
    {
        cerr<<"File "<<sMatrixFileName<<" doesn't exist!"<<endl;
        exit(1);
    }
    sLastLine = "";
}

template <typename  ValueType>
bool ReadMatrix::GetNextMatrix (vector< vector< ValueType > > &  vvMatrix)
{
    bool bSuccessGetMatrix = false;
    bool bValueArea = false; // after reading the first * line, it is true
    int iColumnNumber = 0, i;
    vector < ValueType > vZeroVector;
    vector < ValueType > vRowValues;
    string sCurrentLine;
    if (sLastLine != "")
    {
        if (sLastLine.at(0) == '@')
            iColumnNumber++;
    }
    while (!MatrixFile_stream.eof())
    {
        sCurrentLine.clear();
        getline(MatrixFile_stream, sCurrentLine);
        sLastLine = sCurrentLine;
//        cout << sCurrentLine <<"!"<<endl;
        if (sCurrentLine == "")
            continue;
        if (sCurrentLine.at(0) == '#')
            continue;
        if (sCurrentLine.at(0) == '@')
        {
            if (bValueArea)
                break;
            else
                iColumnNumber++;
        }
        if (sCurrentLine.at(0) == '+')
            if (bValueArea)
                break;
        if (sCurrentLine.at(0) == '*')
        {
            bValueArea = true;
            if (vZeroVector.size() == 0)
                for (i=0; i<iColumnNumber; i++)
                    vZeroVector.push_back(0);
            vRowValues = ReadRowValues( vZeroVector, sCurrentLine);
            vvMatrix.push_back(vRowValues);
        }
    }
    if (vvMatrix.size() > 0)
        bSuccessGetMatrix = true;
    return bSuccessGetMatrix;
}

template <class  ValueType>
vector<ValueType> ReadMatrix::ReadRowValues( vector<ValueType> vZeroVector, const string &  sCurrentLine)
{
    vector<ValueType> vRowValues;
    size_t pos, pos_evl, pos_valueend;
    string sLine, sFieldLine, sId, sValue;
    int i, iColumnId;
    ValueType realValue;

    vRowValues = vZeroVector;
    sLine = sCurrentLine;
    pos = sLine.find("#");
    if (pos != string::npos)
        sLine = sLine.substr(0, pos);
    
    TokenVector words(sLine, " \t\n\r");

    for (i=2; i< (int) words.size(); i++)
    {
        sFieldLine = words[i];
        pos_evl = sFieldLine.find("=");
        pos_valueend = sFieldLine.find_last_not_of("(");
        if (pos_evl == string::npos)
        {
            cerr<<"can't find = in "<<sFieldLine<<endl;
            exit(1);
        }
        sId = sFieldLine.substr(0, pos_evl);
        sValue = sFieldLine.substr(pos_evl+1,  pos_valueend - pos_evl);
        TransferStrToValue(sId, iColumnId);
        TransferStrToValue(sValue, realValue);
        vRowValues[iColumnId] = realValue;
    }

    return vRowValues;
}

template <class  NumberType>
void ReadMatrix::TransferStrToValue(const string sVal, NumberType & realValue)
{
    istringstream input_stream;
    input_stream.clear();
    input_stream.str(sVal);
    input_stream >> realValue;
    input_stream.clear();
}

void ReadMatrix::EndParser()
{
    MatrixFile_stream.clear();
    MatrixFile_stream.close();
}
