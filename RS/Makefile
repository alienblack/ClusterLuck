CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -pedantic -MMD -MP
HIGHS_ROOT ?= /Users/alien/Desktop/optil/heuristic/highs
HIGHS_INC ?= $(HIGHS_ROOT)/highs
HIGHS_BUILD_INC ?= $(HIGHS_ROOT)/build-linux-static
HIGHS_LIB ?= $(HIGHS_ROOT)/build-linux-static/lib/libhighs.a
LDLIBS ?= -lz -lpthread -ldl

SRC := $(shell find src -name '*.cpp' | sort)
OBJ := $(SRC:src/%.cpp=build/%.o)
DEP := $(OBJ:.o=.d)

.PHONY: all clean

all: build/rs_debug

build/rs_debug: $(OBJ)
	$(CXX) $(CXXFLAGS) $^ $(HIGHS_LIB) $(LDLIBS) -o $@

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Isrc -I$(HIGHS_INC) -I$(HIGHS_BUILD_INC) -c $< -o $@

build/rs_submit: CXXFLAGS += -static -s
build/rs_submit: clean $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(HIGHS_LIB) $(LDLIBS) -o $@

clean:
	rm -rf build

-include $(DEP)
