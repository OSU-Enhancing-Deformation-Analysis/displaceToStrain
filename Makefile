# Makefile
#
# This is a simple makefile for building the strain_calc program.
#

# The compiler to use
CC = clang++

# The flags to pass to the compiler
CFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -g

# The directory containing the source code
BUILDDIR = build
SOURCES = strain_calc.cpp
OBJECTS = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(SOURCES))

# The target executable
TARGET = strain_calc

# The rule to build the target
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

# The rule to build each object file
$(BUILDDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean

# The rule to clean the build directory
clean:
	rm -rf $(BUILDDIR)