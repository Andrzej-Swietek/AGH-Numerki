CC = gcc
CFLAGS = -Wall -O3
LIBS = -lgsl -lgslcblas -lm

SRCS = main.c
OBJS = $(SRCS:.c=.o)
EXEC = output.o

# Python script for visualization
VIS_SCRIPT = visualize.py

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

visualize:
	python3 $(VIS_SCRIPT)


