appname := fastNGSadmix

CC := gcc
CXX := g++

FLAGS := -O3 -lz -lpthread

srcfiles := $(shell find -maxdepth 1 -iname "*V2.c")
objects  := $(patsubst %.c, %.o, $(srcfiles))
hfiles  := $(patsubst %.c, %.h, $(srcfiles))

all: $(appname)

$(appname): $(appname).cpp
	$(CXX) $(appname).cpp $(srcfiles) $(FLAGS) -o $(appname)

clean:
	rm  -f $(objects) $(appname)

