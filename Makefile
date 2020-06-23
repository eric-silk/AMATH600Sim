# Referencing:
# https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile
# https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
# Need to decide how "clever" I want to be...

CXX = c++
CXXFLAGS = -I.
RM=rm -f
CPPFLAGS=-g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs)

SRC = hellomake.cpp
SRC += hellofunc.cpp

OBJ=$(subst .cpp,.o,$(SRC))

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

all: hellomake

hellomake: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJ)

distclean: clean
	$(RM) *~ .depend

include .depend
