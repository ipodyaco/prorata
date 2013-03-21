#include <mpi.h>
#include <iostream>
#include <string>
#include <cstdlib>


#include "directoryStructure.h"

#define WORKTAG    1
#define DIETAG     2


char getSeparator()
{
#if _WIN32
		return '\\' ;
#else
		return '/' ;
#endif
}


void searchFastaFiles (vector<string> & vsFastaFilenames, const string & sWorkingDirectory)
{
    int i, iFileNum;
    DirectoryStructure working_dir(sWorkingDirectory);
    working_dir.setPattern(".fasta");
    working_dir.getFiles(vsFastaFilenames);
    working_dir.setPattern(".Fasta");
    working_dir.getFiles(vsFastaFilenames);
    working_dir.setPattern(".FASTA");
    working_dir.getFiles(vsFastaFilenames);
    working_dir.setPattern(".fa");
    working_dir.getFiles(vsFastaFilenames);
    working_dir.setPattern(".faa");
    working_dir.getFiles(vsFastaFilenames);
    iFileNum = (int) vsFastaFilenames.size();
    if (iFileNum == 0)
    {
	cerr<<"no fasta file in the working directory"<<endl;
	exit(1);
    }
    for (i=0; i< iFileNum; i++)
	vsFastaFilenames.at(i) = sWorkingDirectory + getSeparator() + vsFastaFilenames.at(i);
}


/* 
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

void initializeArguments(int argc, char **argv, vector<string> & vsFastaFilenames, 
			 string & sOutputDirectory, string & sprofileDatabasePath,  
			 string & sCommandFirstPart)
{
    int i;

    string sWorkingDirectory, ssoftwarePath, sformatOption, sthreadNumber;
    // Grab command line arguments
    vector<string> vsArguments;
    
    sWorkingDirectory  = "";
    sOutputDirectory   = ""; 
    ssoftwarePath      = ""; 
    sformatOption      = ""; 
    sthreadNumber      = "";
    sprofileDatabasePath = "";    

    while(argc--) 
	vsArguments.push_back(*argv++);
    for(i = 1; i <= (int)vsArguments.size()-1; i++)
	if(vsArguments[i] == "-w") 
		sWorkingDirectory = vsArguments[++i]; 
	else if (vsArguments[i] == "-s")
		ssoftwarePath = vsArguments[++i];
	else if (vsArguments[i] == "-o")
		sOutputDirectory = vsArguments[++i];
	else if (vsArguments[i] == "-f")
		sformatOption = "--" + vsArguments[++i];
	else if (vsArguments[i] == "-c")
		sthreadNumber = "--cpu " + vsArguments[++i];
	else if (vsArguments[i] == "-p")
		sprofileDatabasePath = vsArguments[++i];
	else if ((vsArguments[i] == "-h") || (vsArguments[i] == "--help"))
	{
		cout << "Usage: -w WorkingDirectory, -s software path, -o: output file directory," << endl;
		cout << "-f format option, e.g. notextw, -c thread number,"<<endl;
		cout << "-p profile database path." << endl;
		exit(0);
	}else 
	{
		cerr << "Unknown option " << vsArguments[i] << endl << endl; 
		exit(1);
	}
    if ((sWorkingDirectory == "") || (sOutputDirectory == "") || (ssoftwarePath == "") || (sprofileDatabasePath == ""))
    {
	cerr<<"miss necessary parameter(s)"<<endl;
	exit(1);
    }
	
    searchFastaFiles(vsFastaFilenames, sWorkingDirectory);
    sCommandFirstPart = ssoftwarePath + " "+ sthreadNumber + " " + sformatOption + " ";

}



void MasterProcess(const vector<string> & vsFastaFilenames)
{
    size_t i, workloadSize, iBounderOfProcess;
    int  currentWorkId; //unit id of vWorkLoad
    int  iNumberOfProcessors, iNumberOfSlaves;
    int  result;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &iNumberOfProcessors);     /* get number of processes */
    workloadSize = vsFastaFilenames.size();

    iNumberOfSlaves = iNumberOfProcessors -1;
    iBounderOfProcess = ((workloadSize <= (size_t)iNumberOfSlaves) ? workloadSize : (size_t)iNumberOfSlaves) ;
    for (i=1; i<=iBounderOfProcess; i++)
    {
        currentWorkId = i-1;
        MPI_Send(&currentWorkId, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
    }
    if ((int)workloadSize >  iNumberOfSlaves)
    {
        currentWorkId = iNumberOfSlaves;
        while (currentWorkId < (int) workloadSize)
        {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Send(&currentWorkId, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            currentWorkId ++;
        }
    }
    /* Tell all the slaves to exit by sending an empty message with the DIETAG. */
    for (i=1; i<= (size_t)iNumberOfSlaves; i++)
        MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
//  cout<<"Master process is done."<<endl;
}




void SlaveProcess(const vector<string> & vsFastaFilenames, string & sOutputDirectory, 
		  string & sprofileDatabasePath, string & sCommandFirstPart)
{
    MPI_Status status;
    int currentWorkId, myid;
    string currentWorkFastaFile, sFileName, sFileNameRoot, sCommand;
    size_t pos;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    while(true)
    {
        MPI_Recv(&currentWorkId, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            break;
        currentWorkFastaFile = vsFastaFilenames.at(currentWorkId);
	pos = currentWorkFastaFile.rfind(getSeparator());
	if (pos == string::npos)
	    sFileName = currentWorkFastaFile;
	else
	    sFileName = currentWorkFastaFile.substr(pos+1);
	pos = sFileName.rfind(".");
	if (pos == string::npos)
	    sFileNameRoot = sFileName;
	else
            sFileNameRoot = sFileName.substr(0, pos);
	sCommand = sCommandFirstPart + " -o "+sOutputDirectory+ getSeparator() +sFileNameRoot
			+".hmmsearch.output.txt " + sprofileDatabasePath + " " + currentWorkFastaFile;
	system(sCommand.c_str());
        //cout<<"slave id:"<<myid<<" command: "<<sCommand<<endl;
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    //cout<<"Slave process "<<myid<<" is done."<<endl;
}




int main(int argc, char **argv) 
{
    vector<string> vsFastaFilenames;
    int myid; 
    string sOutputDirectory, sprofileDatabasePath, sCommandFirstPart;
    MPI_Init(&argc,&argv);              /* starts MPI */
    initializeArguments(argc, argv, vsFastaFilenames, sOutputDirectory, sprofileDatabasePath, sCommandFirstPart);

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
    if (myid == 0 )
        MasterProcess(vsFastaFilenames);
    else
        SlaveProcess(vsFastaFilenames, sOutputDirectory, sprofileDatabasePath, sCommandFirstPart);
    MPI_Finalize();          /* let MPI finish up ... */
    //std::cout << "Hello, world!" << std::endl;
    return 0;
}
