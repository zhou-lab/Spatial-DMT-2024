CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -std=gnu11
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
LDLIBS += $(HTSLIB_LIBS) -lz -lm

.PHONY: all clean

all: split_bam_by_cb_hts

split_bam_by_cb_hts: split_bam_by_cb_hts.o Makefile
	$(CC) $(CFLAGS) split_bam_by_cb_hts.o -o $@.tmp $(LDLIBS)
	mv $@.tmp $@

split_bam_by_cb_hts.o: split_bam_by_cb_hts.c Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -f split_bam_by_cb_hts split_bam_by_cb_hts.tmp \
		split_bam_by_cb_hts.o split_bam_by_cb_hts.o.tmp
