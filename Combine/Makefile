CC=	g++
CFLAGS=	 -O2 
Prorata: DriverCombine.cpp projectInfo.cpp proteinCombined.h proteinDirectComparison.h proteinIndirectComparison.h proteinReplicate.h projectInfo.h proteinCombined.cpp proteinDirectComparison.cpp proteinIndirectComparison.cpp proteinReplicate.cpp
	$(CC) $(CFLAGS) -c projectInfo.cpp
	$(CC) $(CFLAGS) -c proteinReplicate.cpp
	$(CC) $(CFLAGS) -c DriverCombine.cpp
	$(CC) $(CFLAGS) -c proteinCombined.cpp
	$(CC) $(CFLAGS) -c proteinDirectComparison.cpp
	$(CC) $(CFLAGS) -c proteinIndirectComparison.cpp
	$(CC) $(CFLAGS) -o Combine DriverCombine.o projectInfo.o proteinCombined.o proteinDirectComparison.o proteinIndirectComparison.o proteinReplicate.o pca.cpp peptideLabelFree.cpp lubksb.o proteomeInfo.o ludcmp.o pcaUtility.o peptideRatio.cpp peakPicker.cpp savgol.o proteinInfo.o smoother.cpp nrutil.o peptideInfo.cpp proteinRatio.o chromatogram.cpp dtaSelectReader.cpp idData.cpp isotopologue.cpp msData.cpp mzReader.cpp scanindex.cpp sicInfo.cpp tokenvector.cpp directoryStructure.cpp proRataConfig.cpp tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp
clean:
	rm Combine
