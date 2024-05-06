CC = g++
CFLAGS = -std=c++17 -Wall
LIBS = -lm

SRCDIR = lib
BINDIR = bin

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o, $(SRCS))
MAIN = application

.PHONY: all clean

all: $(MAIN)

$(MAIN): $(OBJS) main.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

main.o: main.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(BINDIR)/*.o $(MAIN) main.o
