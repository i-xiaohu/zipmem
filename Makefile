CC=			gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2
AR=			ar
DFLAGS=		-DHAVE_PTHREAD
PROG=		build/bwamem
INCLUDES=	bwalib bwamem cstl FM_index
LIBS=		-lm -lz -lpthread
ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

all:$(PROG)

bwalib_src = $(wildcard bwalib/*.c)
bwalib_obj = $(patsubst %.c, %.o, $(bwalib_src))
bwamem_src = $(wildcard bwamem/*.c)
bwamem_obj = $(patsubst %.c, %.o, $(bwamem_src))
cstl_src = $(wildcard cstl/*.c)
cstl_obj = $(patsubst %.c, %.o, $(cstl_src))
FM_index_src = $(wildcard FM_index/*.c)
FM_index_obj = $(patsubst %.c, %.o, $(FM_index_src))

$(PROG):$(bwalib_obj) $(bwamem_obj) $(cstl_obj) $(FM_index_obj)
	$(CC) $(CFLAGS) $(DFLAGS) -o $@ $^ $(LIBS)

$.o:%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(PROG) $(bwalib_obj) $(bwamem_obj) $(cstl_obj) $(FM_index_obj)