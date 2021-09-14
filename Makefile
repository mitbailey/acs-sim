CC=gcc
EDCFLAGS= -std=gnu11 -O2 -Wall -I./
EDLDFLAGS= -lpthread -lm

COBJS=bessel.o \
	datavis.o

all: $(COBJS)
	$(CC) -o acs-datagen.out $(COBJS) $(EDLDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(EDCFLAGS)

.PHONY: clean

clean:
	rm -vf *.o
	rm -vf *.out