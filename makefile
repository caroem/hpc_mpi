# the compiler: gcc for C program, define as g++ for C++
CC = mpicc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -std=c99 -lc -D _BSD_SOURCE 


# the build target executable:
TARGET = gameoflife

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

run: all
	mpirun -n 4 ./$(TARGET)

clean:
	$(RM) $(TARGET)
