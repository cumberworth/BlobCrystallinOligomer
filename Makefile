# Autodependency recipe from
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# But doesn't seem to be working...

OPTLEVEL = -O3
BUILDDIR = build
PREFIX = ../..
TARGET = blobCrystallinOligomer
TESTTARGET = blobCrystallinOligomer_test
TARGETDIR = bin
SRCDIR = src
TESTDIR = test
INCLUDEDIR = include
EXTERNALHEADERS =
EXTERNALLIBS = -lboost_program_options
# If compiling with custom Boost installation,
# add -I/home/amc226/include to $(EXTERNALHEADERS)
# add -L/home/amc226/lib to $(EXTERNALLIBS)

vpath %.h $(INCLUDEDIR)/
vpath %.cpp $(SRCDIR)/ $(TESTDIR)/

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(subst .cpp,.o,$(SOURCES))
OBJECTS := $(subst $(SRCDIR),$(BUILDDIR),$(OBJECTS))

TESTSOURCES := $(wildcard $(TESTDIR)/*.cpp)
TESTOBJECTS := $(subst .cpp,.o,$(TESTSOURCES))
TESTOBJECTS := $(subst $(TESTDIR)/,$(BUILDDIR)/,$(TESTOBJECTS))
TESTOBJECTS += $(filter-out $(BUILDDIR)/main.o, $(OBJECTS))

CPP = g++
CPPFLAGS = -I$(INCLUDEDIR) $(EXTERNALHEADERS) $(OPTLEVEL) --std=c++14 -fPIC
LDFLAGS = $(EXTERNALLIBS) $(OPTLEVEL)

DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.cpp = $(CPP) $(DEPFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: $(TARGETDIR)/$(TARGET)

test: $(TARGETDIR)/$(TESTTARGET)
	$(TARGETDIR)/$(TESTTARGET)

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	$(CPP) $(LDFLAGS) -o $@ $^

$(TARGETDIR)/$(TESTTARGET): $(TESTOBJECTS)
	$(CPP) $(LDFLAGS) -o $@ $^

$(BUILDDIR)/%.o: %.cpp
$(BUILDDIR)/%.o: %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(sources))))

.PHONY: clean install archive
clean:
	rm $(BUILDDIR)/*.o
	rm .d/*.d

install:
	cp $(TARGETDIR)/$(TARGET) $(PREFIX)/bin/$(TARGET)

archive:
	ar -r lib/$(TARGET).a $(BUILDDIR)/*.o
