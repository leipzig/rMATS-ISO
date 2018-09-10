CC      =	gcc
CFLAGS  =	-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
LDFLAGS =
DFLAGS  =	-g -Wall
HTSLIB_DIR = ./htslib-1.4.1
HTSLIB  =   $(HTSLIB_DIR)/libhts.a
LIB     =	$(HTSLIB) -lm -lz -lpthread
COMP_LIB=	-lz
INCLUDE = -I $(HTSLIB_DIR)


BIN_DIR =	.
SRC_DIR =   ./src

HTS_ALL =   hts_all
SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/IsoModule
GDB_DEBUG   =   $(BIN_DIR)/gdb_iso_module
NOR_DEBUG   =   $(BIN_DIR)/debug_iso_module
DMARCRO 	=	-D __DEBUG__

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

#all:       $(SOURCE) $(BIN) 
all:       $(HTS_ALL) $(BIN) 
gdb_iso_module: $(SOURCE) $(GDB_DEBUG) 
debug_iso_module: $(SOURCE) $(NOR_DEBUG)


$(HTS_ALL):
	cd $(HTSLIB_DIR); make;
	
# /usr/lib/../lib64/crt1.o: In function `_start':
# (.text+0x20): undefined reference to `main'
# collect2: error: ld returned 1 exit status
$(BIN): $(OBJS)
		$(CC) -shared $(LDFLAGS) $(OBJS) -o $@ $(LIB)

$(GDB_DEBUG):
		$(CC) $(DFLAGS) $(INCLUDE) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(NOR_DEBUG):
		$(CC) $(CFLAGS) $(INCLUDE) $(SOURCE) $(DMARCRO) -o $@ $(LIB)


clean:
		rm -f $(SRC_DIR)/*.o $(BIN) $(GDB_DEBUG) $(NOR_DEBUG)

clean_debug:
		rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(NOR_DEBUG) 
