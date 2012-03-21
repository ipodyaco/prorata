
CC=gcc
CPP=g++
CFLAGS=-I. -I../../include -g3 -ggdb
CPPFLAGSPARSIFALINC=-I../Libparsifal/include/libparsifal
CPPFLAGS=$(CFLAGS)
LIBNAME = libProRataPRatio.a

TARGET=smoother.o savgol.o pcaUtility.o pca.o sicReaderHelper.o proRataSICReader.o proteinInfo.o peakPicker.o proteomeInfo.o peptideRatio.o proteinRatio.o peptideInfo.o nrutil.o ludcmp.o lubksb.o

lib: $(TARGET) 
	

exec: $(TARGET)
	$(CPP) -o pratio $(TARGET) ../Tinyxml/*.o ../Common/*.o -o pratio ../Libparsifal/temp/lib/libparsifal.a

lib: $(TARGET) 
	ar rc $(LIBNAME) $(TARGET)
	ranlib $(LIBNAME)

DriverProteomeInfo.o:DriverProteomeInfo.cpp $(TARGET)
	$(CPP) $(CPPFLAGS) -c $<

smoother.o:smoother.cpp
	$(CPP) $(CPPFLAGS) -c $<

pca.o:pca.cpp
	$(CPP) $(CPPFLAGS) -c $<

sicReaderHelper.o:sicReaderHelper.cpp
	$(CPP) $(CPPFLAGS) $(CPPFLAGSPARSIFALINC) -c $<

proRataSICReader.o:proRataSICReader.cpp
	$(CPP) $(CPPFLAGS) $(CPPFLAGSPARSIFALINC) -c $<

#mainProRataSICReader.o:mainProRataSICReader.cpp
#	$(CPP) $(CPPFLAGS) -c $<

proteinInfo.o:proteinInfo.cpp
	$(CPP) $(CPPFLAGS) -c $<

peakPicker.o:peakPicker.cpp
	$(CPP) $(CPPFLAGS) -c $<

proteomeInfo.o:proteomeInfo.cpp
	$(CPP) $(CPPFLAGS) -c $<

peptideRatio.o:peptideRatio.cpp
	$(CPP) $(CPPFLAGS) -c $<

proteinRatio.o:proteinRatio.cpp
	$(CPP) $(CPPFLAGS) -c $<

peptideInfo.o:peptideInfo.cpp               
	$(CPP) $(CPPFLAGS) -c $<

savgol.o:savgol.c
	$(CC) $(CFLAGS) -c $<

pcaUtility.o:pcaUtility.c
	$(CC) $(CFLAGS) -c $<

nrutil.o:nrutil.c
	$(CC) $(CFLAGS) -c $<

ludcmp.o:ludcmp.c
	$(CC) $(CFLAGS) -c $<

lubksb.o:lubksb.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(TARGET)
