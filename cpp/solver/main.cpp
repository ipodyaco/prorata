#include <iostream>
#include <string>
#include <cstdlib>

#include "readmatrix.h"

using namespace std;

template <class ValueType>
void PrintMatrix(const vector< vector< ValueType  > >  & vvMatrix)
{
    size_t i, j;
    cout<<"Begin to print matrix"<<endl;
    for (i=0; i<vvMatrix.size(); i++)
    {
        for (j=0; j<vvMatrix[i].size(); j++)
            cout<<vvMatrix[i][j]<<",";
        cout<<endl;
    }
    cout<<"end of matrix"<<endl;
}

int main(int argc, char **argv) 
{
    ReadMatrix myReadMatrix;
    vector < vector <int> >    vvEMatrix;
    vector < vector <double> > vvSMatrix;
    myReadMatrix.InitializeParser("test/testM.txt");

    myReadMatrix.GetNextMatrix (vvEMatrix);
    PrintMatrix(vvEMatrix);
    myReadMatrix.GetNextMatrix(vvSMatrix);
    PrintMatrix(vvSMatrix);
    myReadMatrix.EndParser();
    return 0;
}
