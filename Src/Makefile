# Make file

SHELL = /bin/sh

# System-specific settings
ROOT =	lbnl_processor_latest
EXE =	$(ROOT)_exec
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)
INSTALLDIR = $(HOME)



CC =		g++
CCFLAGS =	-O -I./ -I$(HOME)/include/ -I/usr/local/include
DEPFLAGS =	-M
LINK =		g++
LINKFLAGS =	-O
USRLIB =	-L$(HOME)/lib/ -L/usr/local/lib/ -lgslcblas -lgsl
ARCHIVE =	ar
ARFLAGS =	-rc
SYSLIB =
SIZE =		size


#.DEFAULT:
#	@if [ ! -d Obj_$(EXE) ]; then mkdir Obj_$(EXE); fi
#	@cp -p $(SRC) $(INC) Obj_$(EXE)
#	@cp Makefile Obj_$(EXE)/Makefile
#	@cd Obj_$(EXE); \
#	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
#	##@if [ -d Obj_$(EXE) ]; then cd Obj_$(EXE); rm -f $(SRC) $(INC) Makefile*; fi


# Link Target

exe:	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE).out
	make clean_interms
	$(SIZE) $(EXE).out

# Library Target

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) lib$(EXE).a $(OBJ)
	make clean_interms
	$(SIZE) lib$(EXE).a

all:	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE).out
	$(ARCHIVE) $(ARFLAGS) lib$(EXE).a $(OBJ)
	$(SIZE) $(EXE).out
	$(SIZE) lib$(EXE).a
#make clean_interms
	make install


#Install
install:	$(OBJ)
	@if [ -d $(INSTALLDIR)/bin ]; then cp $(EXE).out $(INSTALLDIR)/bin;\
	else mkdir -p $(INSTALLDIR)/bin;cp $(EXE).out $(INSTALLDIR)/bin; fi
	@if [ -d $(INSTALLDIR)/lib ]; then cp lib$(EXE).a $(INSTALLDIR)/lib;\
	else mkdir -p $(INSTALLDIR)/lib;cp lib$(EXE).a $(INSTALLDIR)/lib; fi
	@if [ -d $(INSTALLDIR)/include/$(EXE) ]; then cp -R *.h $(INSTALLDIR)/include/$(EXE);\
	else mkdir -p $(INSTALLDIR)/include/$(EXE);cp *.h $(INSTALLDIR)/include/$(EXE);\
	fi


static: $(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -static $(EXE).out


clean_interms:
	rm -rf *.o
	rm -rf *.d
clean_all:
	rm -rf *.o
	rm -rf *.d
	rm -rf *.out
	rm -rf lib*.a

# Make tar archive
tar:
	@cd ../; tar cvzf $(ROOT)_Src.tar.gz \
	  Src/Make* Src/*.cpp Src/*.h Src/resources/

	@echo "Created $(ROOT)_src.tar.gz"

# Compilation rules


%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
