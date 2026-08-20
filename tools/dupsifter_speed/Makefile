CC = gcc
CFLAGS = -W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result -O3
CLIB = -lz -lm -llzma -lbz2 -lcurl
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

HTSLIB_ROOT ?= /mnt/isilon/zhoulab/labsoftware/samtools/samtools-1.17/htslib-1.17
LHTSLIB_INCLUDE = $(HTSLIB_ROOT)/htslib
LHTSLIB = $(HTSLIB_ROOT)/libhts.a

# Main program
LIBS = $(LHTSLIB)
dupsifter: $(LIBS) dupsifter.o
	$(CC) $(CFLAGS) dupsifter.o -o $@.tmp -I$(LHTSLIB_INCLUDE) $(LIBS) $(CLIB)
	mv $@.tmp $@

dupsifter.o: dupsifter.c khashl.h refcache.h util.h version.h
	$(CC) -c $(CFLAGS) dupsifter.c -o $@.tmp -I$(LHTSLIB_INCLUDE)
	mv $@.tmp $@

# Clean
.PHONY: clean
clean:
	rm -f dupsifter dupsifter.tmp dupsifter.o dupsifter.o.tmp

purge: clean
	make -C $(LHTSLIB_DIR) clean

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
