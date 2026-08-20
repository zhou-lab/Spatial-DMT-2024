CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Werror -std=gnu11
HTSLIB_ROOT ?= /mnt/isilon/zhoulab/labsoftware/samtools/samtools-1.17/htslib-1.17
CPPFLAGS += -I$(HTSLIB_ROOT)
LDLIBS += $(HTSLIB_ROOT)/libhts.a -lz -lm -llzma -lbz2 -lcurl -lpthread

.PHONY: all clean test

all: cpg_to_yame

test: cpg_to_yame
	./test.sh

cpg_to_yame: cpg_to_yame.o
	$(CC) $(CFLAGS) $^ -o $@.tmp $(LDLIBS)
	mv $@.tmp $@

cpg_to_yame.o: cpg_to_yame.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -f cpg_to_yame cpg_to_yame.tmp cpg_to_yame.o cpg_to_yame.o.tmp
