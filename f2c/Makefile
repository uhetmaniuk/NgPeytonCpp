CC = clang -g
AR = ar rv
RANLIB = ranlib
LD = clang

LIBS = $(METIS) $(LAPACK) $(BLAS)

OBJS = \
assmb.o \
btree2.o \
betree.o \
bfinit.o \
blkfct.o \
blkfc2.o \
blkslv.o \
chordr.o \
chlsup.o \
dscal.o \
epost2.o \
etordr.o \
etpost.o \
etree.o \
fcnthn.o \
fsup1.o \
fsup2.o \
fnsplt.o \
fntsiz.o \
genmmd.o \
getrhs.o \
gtimer.o \
igathr.o \
inpnv.o \
invinv.o \
ldindx.o \
mmdint.o \
mmdelm.o \
mmdupd.o \
mmdnum.o \
mmpy.o \
mmpy1.o \
mmpy2.o \
mmpy4.o \
mmpy8.o \
mmpyi.o \
ordmmd.o \
ordnat.o \
pchol.o \
sfinit.o \
smxpy1.o \
smxpy2.o \
smxpy4.o \
smxpy8.o \
symfct.o \
symfc2.o \
flo2ho.o \
ilo2ho.o \
C2Finterface.o

all: lap2d5pt.x

.c.o:
	$(CC) $(CFLAGS) -std=c99 -pedantic -c $<

lap2d5pt.x: $(OBJS) lap2d5pt.o
	$(LD) $(LDFLAGS) $(OBJS) lap2d5pt.o $(LIBS) -o $@

clean:
	rm -f *.o lap2d5pt.x
