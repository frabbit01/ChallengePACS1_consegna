CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations -I../../pacs-examples/Examples/include

EXEC=challenge_jason
OBJECT=challenge_jason.o
SRC=challenge_jason.cpp

LDFLAGS ?= -L../../pacs-examples/Examples/lib
LIBS  ?= -lmuparser

#$(EXEC): $(OBJECT)
#	g++ $(OBJECT) -o $(EXEC) 
$(EXEC): %: %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

$(OBJECT): $(SRC)
	g++ $(SRC) -o $(OBJECT) -c $(CPPFLAGS)


clear:
	rm *.o $(EXEC)

clean:
	$(RM) *.o

distclean: clean
	$(RM) *~