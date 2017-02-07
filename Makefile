appname := fastNGSadmixV4

CC := gcc
CXX := g++

FLAGS := -O3 -lz -lpthread

srcfiles := $(shell find . -iname "*V2.c")
objects  := $(patsubst %.c, %.o, $(srcfiles))
hfiles  := $(patsubst %.c, %.h, $(srcfiles))

all: $(appname)

$(appname): fastNGSadmixV4.cpp
	$(CXX) fastNGSadmixV4.cpp $(srcfiles) $(FLAGS) -o $(appname)

clean:
	rm  -f $(objects) $(appname)

