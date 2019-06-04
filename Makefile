CC        = gcc
EXEC      = tapyr
CFLAGS    = -Wall -Wextra -Wunused -D_FILE_OFFSET_BITS=64
CGDB      = -g -ggdb -dH -DGDB
CDEBUG    = -DDEBUG
COPTIMIZE = -Wuninitialized -O9 -fomit-frame-pointer
CLIBS     = -lm

CSRCS     = $(wildcard *.c)
CHDRS     = $(wildcard *.h)
TXTS      = $(wildcard *.txt README* LICENSE*)
SCRIPTS   = $(wildcard Makefile* *.sh)

NAME	:= "TAPyR"
VERSION	:= $(shell sed -n 's/.*VERSION \"\(.*\)\".*/\1/p' version.h)
CPUARCH	:= $(shell uname -m)

ifeq ($(MAKECMDGOALS),gdb)
	CFLAGS += $(CGDB)
	EXEC   := $(addsuffix -gdb, $(EXEC))
else ifeq ($(MAKECMDGOALS),debug)
	CFLAGS += $(CGDB) $(CDEBUG)
	EXEC   := $(addsuffix -debug, $(EXEC))
else
	CFLAGS += $(COPTIMIZE)
endif

.PHONY: all clean pack

all: clean bin

gdb: bin

debug: bin

bin:
	@echo :: Compiling \"$(NAME) v$(VERSION)\" \($(CPUARCH)\) ...
	$(CC) $(CFLAGS) $(CSRCS) -o $(EXEC) $(CLIBS)
	@echo :: Done

clean:
	@echo :: Cleaning up ...
	@rm -f $(EXEC) $(EXEC).tar.gz

pack:
	@echo :: Packing files ...
	tar -cvzhf $(EXEC).tar.gz $(CSRCS) $(CHDRS) $(TXTS) $(SCRIPTS)
	@echo :: Done
