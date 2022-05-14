override CFLAGS +=-Wall -Wextra -Wfatal-errors -std=c++14 -MD -MP -O3 -fconcepts 
CPP = g++
LIBS = -static -lpthread -lstdc++ -lboost_program_options -lboost_system -lboost_filesystem -lfmt
SRC = $(wildcard *.cpp)
HEADERS = $(wildcard *.h *.inl)
OBJS = $(SRC:.cpp=.o)
TARGET = commdis

all: $(TARGET)

release: CFLAGS += -DNDEBUG
release: $(TARGET)

pre.h.gch: pre.h
	$(CPP) -O3 -std=c++14 pre.h 
	
$(TARGET): $(OBJS)
	$(CPP) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBS)

-include $(SRC:.cpp=.d)

%.o: %.cpp pre.h.gch
	$(CPP) $(CFLAGS) -include pre.h -c $< -o $@

clean:
	rm -f *.o *.d $(TARGET) pre.h.gch
