CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -std=gnu11
HTSLIB_ROOT ?= /mnt/isilon/zhoulab/labsoftware/samtools/samtools-1.17/htslib-1.17
CPPFLAGS += -I$(HTSLIB_ROOT)
LDLIBS += $(HTSLIB_ROOT)/libhts.a -lz -lm -llzma -lbz2 -lcurl

.PHONY: all clean

all: split_bam_by_cb_hts

split_bam_by_cb_hts: split_bam_by_cb_hts.o
	$(CC) $(CFLAGS) $^ -o $@.tmp $(LDLIBS)
	mv $@.tmp $@

split_bam_by_cb_hts.o: split_bam_by_cb_hts.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -f split_bam_by_cb_hts split_bam_by_cb_hts.tmp \
		split_bam_by_cb_hts.o split_bam_by_cb_hts.o.tmp
