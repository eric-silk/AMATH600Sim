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

print:
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

# For a multi-file case, suppose you have the source files a.F90, b.c, and c.cxx
# (with a program statement appearing in a.F90 or main() appearing in the C or
# C++ source).  This can be built by uncommenting the following two lines.
#
all : power_example

power_example : power_example.cpp
	$(CXX) $(LINK.F) -o $@ $^ $(LDLIBS)

