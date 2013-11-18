project(autostrainer)

set(PATH_TO_CORE "/home/yf1/software/SeqAn/seqan-library-1.4.1")

set(CMAKE_CXX_FLAGS " -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -I${PATH_TO_CORE}/include ")

add_executable(autostrainer ${flags}  main.cpp)

#install(TARGETS autostrainer RUNTIME DESTINATION bin)
