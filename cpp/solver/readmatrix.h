

#ifndef READMATRIX_H_
#define READMATRIX_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>

#include "tokenvector.h"


using namespace std;


class ReadMatrix
{

public:
    ReadMatrix();
    ~ReadMatrix();
    void InitializeParser( const string & sMatrixFileName );

    template <class  ValueType>
    bool GetNextMatrix(vector< vector< ValueType  > > & vvMatrix);
    void EndParser();

private:
    fstream MatrixFile_stream;
    string sLastLine; //the content of the lastLine

    template <class  ValueType>
    vector<ValueType> ReadRowValues( vector< ValueType  > vZeroVector, const string &  sCurrentLine);

    template <class  NumberType>
    void TransferStrToValue(const string sVal, NumberType & realValue);

};

#endif /* READMATRIX_H_  */
