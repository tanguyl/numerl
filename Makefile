
#######################################################################
# Template of for compiling erlang files                              #
#######################################################################
# code to compile
SOURCE = numerl.erl mat.erl numerl_test.erl
NIF_SOURCE = numerl.c
#TESTS = numerl_tests.erl

#Where include files are stored ".hrl"
EFLAGS = -I../include -I/usr/include/gsl

#Compiles the code into a ebin dir. relative to the source dir. 
EBIN = ebin
SRC = src
TARGETS = $(SOURCE:%.erl=$(EBIN)/%.beam)  $(NIF_SOURCE:%.c=$(EBIN)/%_nif.so) 
#$(TESTS:%.erl=$(EBIN)/%.beam)


$(EBIN)/%.beam: $(SRC)/%.erl
	$(ERLHOME)/bin/erlc  -W -b beam -o $(EBIN) $(EFLAGS) $(WAIT) $<

$(EBIN)/%_nif.so: $(SRC)/%.c
	gcc -fpic -shared $(EFLAGS) -o $@ $< -lgsl -lgslcblas -lm 

 
all: $(TARGETS)

test: all
	cd ebin;\
	erl -s numerl_test test -eval 'init:stop()';\
	cd ..;\
   
clean:
	\rm -f $(TARGETS)


realclean: clean
	\rm -f \.* *~ *\% #*  *.beam *.so

.PHONY: test