
#ifndef PRORATAPROTEINTUPLE_H
#define PRORATAPROTEINTUPLE_H

class ProRataProteinTuple
{
	public:

		ProRataProteinTuple(){}
		~ProRataProteinTuple(){}

		QString getLocus() { return qsLocus; }
		QString getDescription() { return qsDescription; }
		double getProteinLog2Ratio() { return dProteinLog2Ratio; }
		double getCIWidth() { return dCIWidth; }

		void setLocus( const QString &qsTemp ) { qsLocus = qsTemp; }
		void setDescription( const QString &qsTemp ) { qsDescription = qsTemp; }
		void setProteinLog2Ratio( double dTemp) { dProteinLog2Ratio = dTemp; }
		void setCIWidth( double dTemp ) { dCIWidth = dTemp; }

	private:

		QString qsLocus;
		QString qsDescription;
		double dProteinLog2Ratio;
		double dCIWidth;

};


#endif  //PRORATAPROTEINTUPLE_H
