TARGET = bin/rt
OBJDIR = bin/objs
SRCDIR = src

CFLAG = -O2 -static -fopenmp 
OBJECTS = main.o

all: $(TARGET)

$(TARGET): $(OBJDIR)/$(OBJECTS)
	g++ -o $@ $(OBJDIR)/$(OBJECTS) 

.PHONY: clean

clean:
	@rm -rf $(OBJDIR)/* $(TARGET)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	g++ $(CFLAG) -c $< -o $@
