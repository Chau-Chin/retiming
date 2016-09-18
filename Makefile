TARGET = bin/rt
OBJDIR = bin/objs
SRCDIR = src

CFLAG = -O2 -static -fopenmp
OBJECTS = $(OBJDIR)/retimeDB.o $(OBJDIR)/main.o 

all: $(TARGET)

$(TARGET): $(OBJECTS)
	g++ -o $@ $(OBJECTS) 

.PHONY: clean

clean:
	@rm -rf $(OBJDIR)/*.o $(TARGET)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	g++ $(CFLAG) -c $< -o $@
