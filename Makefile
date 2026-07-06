CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -MMD -MP -static -s
BUILD := build
HIGHS_ROOT := $(CURDIR)/third_party/highs
HIGHS_BUILD := $(HIGHS_ROOT)/build-linux-static
HIGHS_LIB := $(HIGHS_BUILD)/lib/libhighs.a
CADICAL_ROOT := $(CURDIR)/third_party/cadical
CADICAL_LIB := $(CADICAL_ROOT)/build/libcadical.a

AGGRESSIVE_FLAGS := \
	-DSTRICT_EXACT_DEFAULT=1 \
	-DAUTO_TARGET_LOOP_DEFAULT=1 \
	-DALLOW_STABLE_EXCHANGE_DEFAULT=1 \
	-DALLOW_EMPIRICAL_EXACT_DEFAULT=1 \
	-DDEFAULT_TOTAL_SECONDS=1740.0

.PHONY: all clean

all: $(BUILD)/manytree_solver $(BUILD)/twotree_solver

$(BUILD):
	mkdir -p $@

$(CADICAL_LIB):
	cd $(CADICAL_ROOT) && ./configure
	$(MAKE) -C $(CADICAL_ROOT) -j1

$(BUILD)/manytree_solver: src_manytree/manytree_pool_packer.cpp $(CADICAL_LIB) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(AGGRESSIVE_FLAGS) -DUSE_CADICAL -I$(CADICAL_ROOT)/src $< $(CADICAL_LIB) -o $@

$(HIGHS_LIB):
	cmake -S $(HIGHS_ROOT) -B $(HIGHS_BUILD) -DBUILD_SHARED_LIBS=OFF -DFAST_BUILD=ON -DCMAKE_BUILD_TYPE=Release
	cmake --build $(HIGHS_BUILD) --target libhighs --parallel 1

$(BUILD)/twotree_solver: $(HIGHS_LIB) | $(BUILD)
	$(MAKE) -C RS build/rs_submit HIGHS_ROOT=$(HIGHS_ROOT) CXXFLAGS="$(CXXFLAGS)"
	cp RS/build/rs_submit $@

clean:
	rm -rf $(BUILD)
	$(MAKE) -C RS clean
