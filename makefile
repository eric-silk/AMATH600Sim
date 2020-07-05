# -*- mode: makefile -*-

# This sample (GNU) Makefile can be used to compile PETSc applications with a single
# source file and can be easily modified to compile multi-file applications.
# It relies on pkg_config tool, and PETSC_DIR and PETSC_ARCH variables.
# Copy this file to your source directory as "Makefile" and modify as needed.
#
# For example - a single source file can be compiled with:
#
#  $ cd src/snes/tutorials/
#  $ make -f $PETSC_DIR/share/petsc/Makefile.user ex17
#
# The following variable must either be a path to petsc.pc or just "petsc" if petsc.pc
# has been installed to a system location or can be found in PKG_CONFIG_PATH.
PETSC_ARCH := arch-linux-c-debug
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm
LDLIBS += -lstdc++

RM=rm -f

# Many suffixes are covered by implicit rules, but you may need to write custom rules
# such as these if you use suffixes that do not have implicit rules.
# https://www.gnu.org/software/make/manual/html_node/Catalogue-of-Rules.html#Catalogue-of-Rules
% : %.F90
	$(LINK.F) -o $@ $^ $(LDLIBS)
%.o: %.F90
	$(COMPILE.F) $(OUTPUT_OPTION) $<
% : %.cxx
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cxx
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
% : %.cpp
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cpp
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
% : %.c
	$(LINK.c) -o $@ $^ $(LDLIBS)
%.o: %.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

SRC := my_power.c
SRC += my_pffunctions.c
SRC += my_PFReadData.c

OBJ := $(subst .c,.o,$(SRC))

all : my_power

print :
	@echo $(COMPILE.c)

my_power : $(OBJ)
	  $(LINK.c) -o $@ $^ $(LDLIBS)

.PHONY: clean
clean :
	$(RM) $(OBJ) my_power

print_vars:
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo FC=$(FC)
	@echo CFLAGS=$(CFLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo FFLAGS=$(FFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)
	@echo LINK.F=$(LINK.F)
	@echo 
	@echo LINK.cc=$(LINK.cc)

