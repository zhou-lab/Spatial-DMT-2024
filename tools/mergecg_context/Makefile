CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Werror -std=gnu11
LDLIBS += -lm

.PHONY: all clean

all: mergecg_context

mergecg_context: mergecg_context.o
	$(CC) $(CFLAGS) $^ -o $@.tmp $(LDLIBS)
	mv $@.tmp $@

mergecg_context.o: mergecg_context.c
	$(CC) $(CFLAGS) -c $< -o $@.tmp
	mv $@.tmp $@

clean:
	rm -f mergecg_context mergecg_context.tmp \
		mergecg_context.o mergecg_context.o.tmp
