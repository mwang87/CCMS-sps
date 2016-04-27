################################################################################
# Definitions
################################################################################
BASE_DIR := ../../..
#-------------------------------------------------------------------------------
include $(BASE_DIR)/Makefiles/Makefile.defs

################################################################################
# Conditions
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.cond

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files - SPS_SERVER
################################################################################
SPS_SERVER_SRC = \
	ReportServerInterface.cpp \
	spsReportServer.cpp

SPS_SERVER_HDR = \
	ReportServerInterface.h \
	spsReportServer.h

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
  $(EXEC_DIR)/CommandLineParser.cpp \
  $(EXEC_DIR)/ParameterList.cpp  \
  $(COMMON_FILES)/abruijn.cpp \
  $(COMMON_FILES)/ms1.cpp \
	$(COMMON_FILES)/SpectrumPair.cpp \
	$(COMMON_FILES)/SpectrumPairSet.cpp

CPPLIB_HDR = \
  $(EXEC_DIR)/CommandLineParser.h \
  $(EXEC_DIR)/ParameterList.h  \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h

################################################################################
# Compile rules - specific
################################################################################
SPS_SERVER_OBJ=$(SPS_SERVER_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

spsReportServer: $(SPS_SERVER_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC)  $(LIBSPSREPORTS_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(SPS_SERVER_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC)  $(LIBSPSREPORTS_SRC) -o $@ $(LIBS)

################################################################################
# Compile rules - target
################################################################################
EXECUTABLES = \
	spsReportServer

all: $(EXECUTABLES)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################

