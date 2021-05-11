UNAME := $(shell uname)

SRCS = Span.cpp Data.cpp Segment.cpp Fraction.cpp
OBJS = $(SRCS:.cpp=.o)


CC = g++
LDFLAGS= #-fopenmp
CFLAGS= -O3 -std=c++11 #-fopenmp

ifeq ($(UNAME), Linux)
LDFLAGS += -fopenmp
CFLAGS += -fopenmp
endif

all : span.run 

span.run: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)


%.o : %.cpp
	 $(CC) $(CFLAGS) -c $*.cpp

clean:
	-rm *.o span.run

depend:
	makedepend -Y -- $(CFLAGS) -- $(SRCS)
# DO NOT DELETE

Segment.o: Segment.h Data.h defs.h Fraction.h
Data.o: Data.h defs.h Fraction.h
Span.o: Data.h Segment.h Wall_Time.h Fraction.h
Fraction.o: Fraction.h defs.h
