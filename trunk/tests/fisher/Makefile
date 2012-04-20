VERSION:= 0.0.0

SRCDIR:= ./

EXE_A:= fisher

SRC_A:= main.cpp fisher2.cpp gw_maths.cpp

OBJS_A:= $(SRC_A:.cpp=.o)

CC:= gcc		# The C compiler.
CFLAGS:= -O3 -Wall	
					# C compilation options which relate to
					# optimization or debugging (usually
    			# just -g or -O).  Usually this wouldn't
    			# include -I options to specify the
    			# include directories, because then you
    			# couldn't override it on the command line
    			# easily as in the above example.
#CXX:= icpc     	
CXX:= g++     	
					# The C++ compiler.  
					# (Sometimes "CPP" instead of CXX.)
#CXXFLAGS:= -O3 -Wall -wd981	-wd1572 -wd383 -wd1418 
CXXFLAGS:= -O3 -Wall 
					# C++ compilation options related to 
    			# optimization or debugging (-O or -g).



INCLUDEDIR:= -I/usr/local/include 
LIBDIR:= -L/usr/local/lib 
				# Contains libraries we need to link in.
LIBS:= -lgsl -lgslcblas -lm

help : 
	@echo " "
	@echo "$(EXE_A) version $(VERSION) source code"
	@echo " "
	@echo "Type ...          To ..."
	@echo "make all          Compile the program"
	@echo "make clean        Delete temporary files"
	@echo " "

all: $(EXE_A)  
$(EXE_A): $(OBJS_A)
	$(CXX) $(INCLUDEDIR) $(OBJS_A) -o $(EXE_A) $(LIBDIR) $(LIBS)

.cpp.o: $*.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDEDIR) -o $@ -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)


clean:
	rm -rf $(SRCDIR)*.o $(SRCDIR)*~
