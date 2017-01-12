SOURCES=driver.c linterp_many.c readtable.c linterp_some.c \
	findtemp.c linterp_for_temp.c testing.c nuc_eos.c
INCLUDES=nuc_eos.h
HDF5DIR=/opt/local
HDF5INCS=-I$(HDF5DIR)/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5
OBJECTS=$(SOURCES:.c=.o )

CC=gcc
CFLAGS=-g -O0 -std=gnu99
EXTRALIBS=-lm

driver: $(OBJECTS) $(INCLUDES)
	$(CC) $(CFLAGS) -o driver $(OBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) $(HDF5INCS) -c $< -o $@

clean:
	rm -f *.o
