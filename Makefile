SRC = sim.cpp
TARGET = $(SRC:.cpp=)

OS_NAME:=$(shell uname -s | tr A-Z a-z)
ifeq ($(OS_NAME),darwin)
STDINCDIR := -I/opt/local/include
STDLIBDIR := -L/opt/local/lib
else
STDINCDIR := 
STDLIBDIR := 
endif

LD = $(shell root-config --ld)
CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g -Wall -O2

$(TARGET): $(SRC)
	$(LD) $(CPPFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)