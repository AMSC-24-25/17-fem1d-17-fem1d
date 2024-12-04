# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -g

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build

# Source and object files
SRCFILES = $(wildcard $(SRCDIR)/*.cpp)
OBJFILES = $(SRCFILES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Output executable
TARGET = my_project.exe

# Default target
all: clean $(TARGET)

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
