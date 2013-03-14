#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>
 
#include "CoinTime.hpp"
#include "CoinError.hpp"
 
#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "MyTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"
//#define REDIRECT

#include <iostream>
#include <string>
#include <cstdlib>

#include "readmatrix.h"

using namespace std;
using namespace Ipopt;
using namespace Bonmin;

template <class ValueType>
void PrintMatrix(const std::vector< std::vector< ValueType  > >  & vvMatrix)
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

bool GetNextCase(ReadMatrix & myReadMatrix, std::vector< std::vector<int> > & vvEMatrix, std::vector< std::vector<double> > & vvSMatrix)
{
    bool bSuccess;
    
    // memory is not really released by calling clear()
    vvEMatrix.clear();
    vvSMatrix.clear();
    bSuccess = myReadMatrix.GetNextMatrix (vvEMatrix);
    if (bSuccess)
    {
        PrintMatrix(vvEMatrix);
        bSuccess = myReadMatrix.GetNextMatrix(vvSMatrix);
        if (bSuccess)
            PrintMatrix(vvSMatrix);
    }
    return bSuccess;
}

void SolveCurrentCase(int PathNum,  const std::vector< std::vector<int> > & vvEMatrix, const std::vector< std::vector<double> > & vvSMatrix)
{
    WindowsErrorPopupBlocker();
    SmartPtr<MyTMINLP> tminlp = new MyTMINLP(PathNum, vvEMatrix, vvSMatrix); 
#ifdef REDIRECT
    FILE * fp = fopen("log.out","w");
    CoinMessageHandler handler(fp);
    BonminSetup bonmin(&handler);
#else
    BonminSetup bonmin;
#endif
    bonmin.initializeOptionsAndJournalist();
    //Register an additional option
    bonmin.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                  "yes",
                                  "no", "No, we don't.",
                                  "yes", "Yes, we do.",
                                  "A longer comment can be put here"); 

   // Here we can change the default value of some Bonmin or Ipopt option
   bonmin.options()->SetNumericValue("bonmin.time_limit", 500); //changes bonmin's time limit
//   bonmin.options()->SetStringValue("mu_oracle","loqo");
 
   //Here we read several option files
   bonmin.readOptionsFile("test/Mybonmin.opt");
   bonmin.readOptionsFile();// This reads the default file "bonmin.opt"
 
   // Options can also be set by using a string with a format similar to the bonmin.opt file
   bonmin.readOptionsString("bonmin.algorithm B-BB\n");
 
   // Now we can obtain the value of the new option
   int printSolution;
   bonmin.options()->GetEnumValue("print_solution", printSolution,"");
   if(printSolution == 1)
   {
        tminlp->printSolutionAtEndOfAlgorithm();
   } 

    //Now initialize from tminlp
    bonmin.initialize(GetRawPtr(tminlp)); 

    //Set up done, now let's branch and bound
    double time1 = CoinCpuTime();
    try 
    {
        Bab bb;
        cout<<"begin solver"<<endl;
        bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc 
        cout<<"end solver"<<endl;
    }
    catch(TNLPSolver::UnsolvedError *E) 
    {
        cerr<<"wrong1"<<endl;
        //There has been a failure to solve a problem with Ipopt.
        std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
    }
    catch(OsiTMINLPInterface::SimpleError &E) 
    {
        cerr<<"wrong2"<<endl;
        std::cerr<<E.className()<<"::"<<E.methodName()
                    <<std::endl
                    <<E.message()<<std::endl;
    }
    catch(CoinError &E)
    {
        cerr<<"wrong3"<<endl;
        std::cerr<<E.className()<<"::"<<E.methodName()
                  <<std::endl
                  <<E.message()<<std::endl;
    } 

}

int main(int argc, char **argv) 
{
    ReadMatrix myReadMatrix;
    std::vector < std::vector <int> >    vvEMatrix;
    std::vector < std::vector <double> > vvSMatrix;
    myReadMatrix.InitializeParser("test/testM.txt");
    if (GetNextCase(myReadMatrix, vvEMatrix, vvSMatrix))
    {
        SolveCurrentCase(2, vvEMatrix, vvSMatrix);
    }
    else
        myReadMatrix.EndParser();
    return 0;
}
