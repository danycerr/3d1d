SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=branch_gen

CXXFLAGS=-std=c++14 
CXXFLAGS+= -lm 
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
