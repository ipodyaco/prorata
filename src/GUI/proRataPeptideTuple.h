

#ifndef PRORATAPEPTIDETUPLE_H
#define PRORATAPEPTIDETUPLE_H

class ProRataPeptideTuple
{
	public:

		ProRataPeptideTuple(){}
		~ProRataPeptideTuple(){}


		QString getChroFilename(){ return  qsChroFilename; }
		QString getChroIdentifier() { return qsChroIdentifier; }
		QString getChroValidity(){ return  qsChroValidity; }
		QString getChroSequence(){ return  qsChroSequence; }
		double getChroLog2Ratio(){ return  dChroLog2Ratio; }
		double getEigenvalueRatio(){ return  dEigenvalueRatio; }
		int getNTerminalDistance(){ return  iNTerminalDistance; }
		int getCTermicalDistance(){ return  iCTerminalDistance; }
		int getChargeState(){ return  iChargeState; }

		void setChroFilename( QString qsTemp ) { qsChroFilename = qsTemp; }
		void setChroIdentifier( QString qsTemp ) { qsChroIdentifier = qsTemp; }
		void setChroValidity(  QString qsTemp ) { qsChroValidity = qsTemp; }
		void setChroSequence( QString qsTemp ) { qsChroSequence = qsTemp; }
		void setChroLog2Ratio( double dTemp ) { dChroLog2Ratio = dTemp; }
		void setEigenvalueRatio( double dTemp ) { dEigenvalueRatio = dTemp; }
		void setNTerminalDistance( int iTemp ) { iNTerminalDistance = iTemp; }
		void setCTerminalDistance( int iTemp ) { iCTerminalDistance = iTemp; }
		void setChargeState( int iTemp ) { iChargeState = iTemp; }


	private:

		QString qsChroFilename;
		QString qsChroIdentifier;
		QString qsChroValidity;
		QString qsChroSequence;
		double dChroLog2Ratio;
		double dEigenvalueRatio;
		int iNTerminalDistance;
		int iCTerminalDistance;
		int iChargeState;

};


#endif  //PRORATAPEPTIDETUPLE_H
