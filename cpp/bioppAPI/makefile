CXX=g++
ODIR=obj
BDIR=bin
SDIR=source
IDIR=include

BPP_INCLUDE=$(HOME)/local/bpp/dev/include
BPP_LIB=$(HOME)/local/bpp/dev/lib

BOOST_INCLUDE=/usr/include
BOOST_LIB=/usr/lib



CPPFLAGS= -O3 -I$(IDIR) -I$(BPP_INCLUDE) -I$(BOOST_INCLUDE) -std=c++11
LDFLAGS=-L$(BPP_LIB) -L$(BOOST_LIB)
LIBS=-lm -lbpp-core -lbpp-seq -lbpp-phyl


RECPHYLOXMLBPPAPI_FILES = main.cpp recPhyloXMLIO.cpp ReconciledTree.cpp XMLUtils.cpp
RECPHYLOXMLBPPAPI_OBJS = $(subst .cpp,.o,$(RECPHYLOXMLBPPAPI_FILES))
RECPHYLOXMLBPPAPI_OBJS := $(addprefix $(ODIR)/, $(RECPHYLOXMLBPPAPI_OBJS))
RECPHYLOXMLBPPAPI_SRCS = $(addprefix $(SDIR)/, $(RECPHYLOXMLBPPAPI_FILES))

_EXES= testRecPhyloXMLBppAPI
EXES=$(patsubst %,$(BDIR)/%,$(_EXES))


all: $(EXES)

debug: CPPFLAGS = -O3 -I$(IDIR) -I$(BPP_INCLUDE) -I$(BOOST_INCLUDE) -g -std=c++11
debug: all


$(BDIR)/testRecPhyloXMLBppAPI: $(RECPHYLOXMLBPPAPI_OBJS)
	mkdir -p bin 
	$(CXX) -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LIBS)


#general rules

$(ODIR)/%.o: $(SDIR)/%.cpp
	mkdir -p $(ODIR)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -c -MD -o $@ $<
	@cp $(ODIR)/$*.d $(ODIR)/$*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $(ODIR)/$*.d >> $(ODIR)/$*.P; \
	rm -f $(ODIR)/$*.d

$(BDIR)/%: $(ODIR)/%.o
	mkdir -p bin 
	$(CXX) -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LIBS)



all: $(EXES)


clean:
		rm -f $(ODIR)/*.o *~ $(SDIR)/*~ core  
		rm -rf $(ODIR) 
