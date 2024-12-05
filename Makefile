# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -g -I${mkEigenInc}

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build

# Source and object files
SRCFILES = $(wildcard $(SRCDIR)/*.cpp)
OBJFILES = $(SRCFILES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Output executable
TARGET = fem1d.exe

# Default target
all: $(TARGET)

# Clean build (rebuild all objects files)
rebuild: clean all

# Link the object files to create the executable
$(TARGET): $(OBJFILES)
	$(CXX) $(OBJFILES) -o $@
# ./$(TARGET)

# Compile source files into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Clean up the build directory
clean:
	rm -f $(BUILDDIR)/*.o $(TARGET)

.PHONY: all clean
