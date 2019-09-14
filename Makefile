#!/bin/make

CC = icc
OPTS = -O2 -mkl=sequential  -D__SERIAL__ -D__FLOAT64__ -D__MKL__ -D__CHECKPOINT__ 

BINDIR=bin
SRCDIR=src

BINSRC = $(BINDIR)/weak_allee/ABCASMC_stochAllee.c $(BINDIR)/weak_allee/ABCASMC_detAllee.c $(BINDIR)/weak_allee/ABCAPCSMC_Allee.c $(BINDIR)/weak_allee/ABCAMMSMC_Allee.c $(BINDIR)/scratch_assay/ABCASMC_detFKPP.c $(BINDIR)/scratch_assay/ABCASMC_stochFKPP.c $(BINDIR)/scratch_assay/ABCAPCSMC_FKPP.c  $(BINDIR)/scratch_assay/ABCAMMSMC_FKPP.c

SRC = $(SRCDIR)/mcl.c $(SRCDIR)/SSAL.c $(SRCDIR)/abc/dabcrs.c $(SRCDIR)/abc/dabcasmc.c $(SRCDIR)/abc/dabcapcsmc.c $(SRCDIR)/abc/dabcammsmc.c $(SRCDIR)/sim/dalrws.c $(SRCDIR)/sim/drkf45s.c $(SRCDIR)/sim/dbtcsfpecs.c $(SRCDIR)/util/choldc.c $(SRCDIR)/util/cholfs.c $(SRCDIR)/util/mvprod.c $(SRCDIR)/util/copyDataset.c $(SRCDIR)/util/dmcint.c $(SRCDIR)/util/dmcintcv.c $(SRCDIR)/util/dmcintd.c $(SRCDIR)/util/dmcintv.c $(SRCDIR)/util/suarngs.c $(SRCDIR)/util/durngus.c $(SRCDIR)/util/durngns.c $(SRCDIR)/util/durngmvns.c $(SRCDIR)/util/durngpmfs.c

OBJS = $(SRC:.c=.o)
BINOBJS=$(BINSRC:.c=.o)
BIN = $(BINSRC:.c=)
INCDIR = ./include
INC = -I $(INCDIR) 
LIBS = -limf 

.SUFFIXES: .c .o

all: $(BIN) 
	@echo Binaries $(BIN) created!

.c.o: 
	$(CC) $(OPTS) -c $< -o $@ $(INC) 

$(BIN): $(BINOBJS) $(OBJS) 
	$(CC) $(OPTS) -o $@  $(INC) $@.o $(OBJS) $(LIBS)
clean:
	set nonomatch; rm -f $(BIN) $(BINOBJS) $(OBJS) 
