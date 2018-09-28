appname := fastNGSadmix

CC := gcc
CXX := g++

FLAGS := -O3 -lz

srcfiles := $(shell find -maxdepth 1 -iname "*V3.cpp")
objects  := $(patsubst %.c, %.o, $(srcfiles))
hfiles  := $(patsubst %.c, %.h, $(srcfiles))

all: $(appname)

$(appname): $(appname).cpp
	$(CXX) $(appname).cpp $(srcfiles) $(FLAGS) -o $(appname)

clean:
	rm  -f $(objects) $(appname)

