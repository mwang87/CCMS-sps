################################################################################
# Definitions
################################################################################
BASE_DIR := ../..
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
# Libraries
################################################################################
LIBRARIES = -lpthread -lspecnets -lcore -lSpsReports -lxerces-c -lhpdf -lpng -lz -lm

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	main_comparepsms \
	main_convertpsm \
	main_createdecoydb \
	main_findtaglen \
	main_havetag \
	util_parsimony \
	main_scorecontig \
	main_specdump
	#main_havetag

.PHONY: subdirs $(SUBDIRS) clean all

ifeq ($(type),debug)
all: $(EXECUTABLES)
else
all: $(EXECUTABLES)
	$(STRIP) $(EXECUTABLES)
endif

################################################################################
# Compile rules: Cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES)

################################################################################
# Compile rules: Executables
################################################################################
# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c-3.0 $(LIBS)

# Windows
else
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c-3-0 $(LIBS)
endif

main_comparepsms: main_comparepsms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_comparepsms.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_comparepsms $(LIBS) $(PWIZ_LIBS)

main_convertpsm: main_convertpsm.cpp  $(LIBRARIES_SRC)
	$(CXX) main_convertpsm.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_convertpsm $(LIBS) $(PWIZ_LIBS)

main_createdecoydb: main_createdecoydb.cpp  $(LIBRARIES_SRC)
	$(CXX) main_createdecoydb.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_createdecoydb $(LIBS) $(PWIZ_LIBS)

main_findtaglen: main_findtaglen.cpp  $(LIBRARIES_SRC)
	$(CXX) main_findtaglen.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_findtaglen $(LIBS) $(PWIZ_LIBS)

main_havetag: main_havetag.cpp  $(LIBRARIES_SRC)
	$(CXX) main_havetag.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_havetag $(LIBS) $(PWIZ_LIBS)

main_scorecontig: main_scorecontig.cpp  $(LIBRARIES_SRC)
	$(CXX) main_scorecontig.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_scorecontig $(LIBS) $(PWIZ_LIBS)

main_specdump: main_specdump.cpp  $(LIBRARIES_SRC)
	$(CXX) main_specdump.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_specdump $(LIBS) $(PWIZ_LIBS)

util_parsimony: util_parsimony.cpp  $(LIBRARIES_SRC)
	$(CXX) util_parsimony.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o util_parsimony $(LIBS) $(PWIZ_LIBS)

################################################################################
# Compile rules: Dependencies
################################################################################
# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
