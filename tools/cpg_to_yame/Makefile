CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Werror -std=gnu11
PKG_CONFIG ?= pkg-config
HTSLIB_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
HTSLIB_LIBS ?= $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)
HTSLIB_LIBDIR ?= $(shell $(PKG_CONFIG) --variable=libdir htslib 2>/dev/null)
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS = -lhts
endif
ifneq ($(strip $(HTSLIB_LIBDIR)),)
HTSLIB_LIBS += -Wl,-rpath,$(HTSLIB_LIBDIR)
endif
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS) -lz -lm -lpthread

.PHONY: all clean test

all: cpg_to_yame

test: cpg_to_yame
	./test.sh

cpg_to_yame: cpg_to_yame.o Makefile
	$(CC) $(CFLAGS) cpg_to_yame.o -o $@.tmp $(LDLIBS)
	mv $@.tmp $@

cpg_to_yame.o: cpg_to_yame.c Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -f cpg_to_yame cpg_to_yame.tmp cpg_to_yame.o cpg_to_yame.o.tmp
