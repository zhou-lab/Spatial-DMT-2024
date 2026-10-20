CC ?= cc
CFLAGS = -W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result -O3
CLIB = -lz -lm
CF_OPTIMIZE = 1

OS := $(shell uname)
ifeq ($(OS),  Darwin)
	CFLAGS += -Wno-unused-function
endif

INCLUDE = include

########################
### different modes ####
########################

PROG = dupsifter

.PHONY : all profile debug build

all: build

build: exportcf $(PROG)

profile: CF_OPTIMIZE := 0
profile: CFLAGS += -pg
profile: CFLAGS := $(filter-out -O3,$(CFLAGS))
profile: build

debug: CF_OPTIMIZE := 0
debug: CFLAGS += -g
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: build

exportcf:
	$(eval export CF_OPTIMIZE)

#####################
##### libraries #####
#####################

PKG_CONFIG ?= pkg-config
HTSLIB_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
HTSLIB_LIBS ?= $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)
HTSLIB_INCLUDEDIR ?= $(shell $(PKG_CONFIG) --variable=includedir htslib 2>/dev/null)
HTSLIB_LIBDIR ?= $(shell $(PKG_CONFIG) --variable=libdir htslib 2>/dev/null)
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS = -lhts
endif
ifneq ($(strip $(HTSLIB_LIBDIR)),)
HTSLIB_LIBS += -Wl,-rpath,$(HTSLIB_LIBDIR)
endif
ifneq ($(strip $(HTSLIB_INCLUDEDIR)),)
HTSLIB_CFLAGS += -I$(HTSLIB_INCLUDEDIR)/htslib
endif

# Main program
dupsifter: dupsifter.o Makefile
	$(CC) $(CFLAGS) dupsifter.o -o $@.tmp $(HTSLIB_LIBS) $(CLIB)
	mv $@.tmp $@

dupsifter.o: dupsifter.c khashl.h refcache.h util.h version.h Makefile
	$(CC) -c $(CFLAGS) $(HTSLIB_CFLAGS) dupsifter.c -o $@.tmp
	mv $@.tmp $@

# Clean
.PHONY: clean
clean:
	rm -f dupsifter dupsifter.tmp dupsifter.o dupsifter.o.tmp

purge: clean

## clean to make a release zip
.PHONY: release
release:
	rm -rf release-source.zip dupsifter-release
	git clone . dupsifter-release
	make -C dupsifter-release cleanse
	zip -r release-source.zip dupsifter-release
	rm -rf dupsifter-release

# removes git history, for release internal use
cleanse : purge
	rm -f **/*.o .gitignore
	rm -rf .git
