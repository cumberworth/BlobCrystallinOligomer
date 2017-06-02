# Autodependency recipe from
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# But doesn't seem to be working...

OPTLEVEL = -g
BUILDDIR = build
PREFIX = ../../
TARGET = blobCrystallinOligomer
TARGETDIR = bin
SRCDIR = src
INCLUDEDIR = include
EXTERNALHEADERS = -Iinclude -I/usr/include/eigen3
EXTERNALLIBS = -lboost_program_options -lboost_mpi -lboost_serialization

vpath %.h $(INCLUDEDIR)/
vpath %.cpp $(SRCDIR)/

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(subst .cpp,.o,$(SOURCES))
OBJECTS := $(subst $(SRCDIR),$(BUILDDIR),$(OBJECTS))

CPP = mpicxx
CPPFLAGS = -I$(INCLUDEDIR) $(EXTERNALHEADERS) $(OPTLEVEL)
LDFLAGS = $(EXTERNALLIBS) $(OPTLEVEL)

# For compiling on Dexter (using local Boost installation (I think that's why this is needed))
#CPPFLAGS = -I/home/amc226/include $(EXTERNALHEADERS) $(OPTLEVEL)
#LDFLAGS = -L/home/amc226/lib $(EXTERNALLIBS) $(OPTLEVEL)

DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.cpp = $(CPP) $(DEPFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: $(TARGETDIR)/$(TARGET)

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	$(CPP) $(LDFLAGS) -o $@ $^ 

$(BUILDDIR)/%.o: %.cpp
$(BUILDDIR)/%.o: %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(sources))))

.PHONY: clean install
clean:
	rm $(BUILDDIR)/*.o

install:
	cp $(TARGET) $(PREFIX)/bin/$(TARGET)
