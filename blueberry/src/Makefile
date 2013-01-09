CC=	g++
CFLAGS=	 -O2
cc=	gcc
Prorata: DriverProteomeInfo.cpp pca.cpp peptideLabelFree.cpp proteinRatio.h lubksb.c pca.h peptideLabelFree.h proteomeInfo.cpp ludcmp.c pcaUtility.c peptideRatio.cpp proteomeInfo.h peakPicker.cpp peptideRatio.h savgol.c nr.h peakPicker.h proteinInfo.cpp smoother.cpp nrutil.c peptideInfo.cpp proteinInfo.h smoother.h nrutil.h peptideInfo.h proteinRatio.cpp chromatogram.cpp dtaSelectReader.cpp idData.cpp isotopologue.cpp msData.cpp mzReader.cpp scanindex.cpp sicInfo.cpp tokenvector.cpp chromatogram.h dtaSelectReader.h idData.h isotopologue.h msData.h mzReader.h scanindex.h sicInfo.h tokenvector.h directoryStructure.cpp directoryStructure.h proRataConfig.cpp proRataConfig.h tinystr.cpp tinystr.h tinyxml.cpp tinyxmlerror.cpp tinyxml.h tinyxmlparser.cpp
	$(cc) -c pcaUtility.c
	$(cc) -c savgol.c
	$(cc) -c nrutil.c
	$(cc) -c ludcmp.c
	$(cc) -c lubksb.c
	$(CC) $(CFLAGS) -c proteomeInfo.cpp
	$(CC) $(CFLAGS) -c proteinInfo.cpp
	$(CC) $(CFLAGS) -c proteinRatio.cpp
	$(CC) $(CFLAGS) -c proRataConfig.cpp
	$(CC) $(CFLAGS) -c tinyxml.cpp
	$(CC) $(CFLAGS) -o Prorata DriverProteomeInfo.cpp pca.cpp peptideLabelFree.cpp lubksb.o proteomeInfo.o ludcmp.o pcaUtility.o peptideRatio.cpp peakPicker.cpp savgol.o proteinInfo.o smoother.cpp nrutil.o peptideInfo.cpp proteinRatio.o chromatogram.cpp dtaSelectReader.cpp idData.cpp isotopologue.cpp msData.cpp mzReader.cpp scanindex.cpp sicInfo.cpp tokenvector.cpp directoryStructure.cpp proRataConfig.cpp tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

clean:
	rm Prorata *.o

