CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -MMD -MP -static -s
BUILD := build
EMBED := $(BUILD)/embed
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
	-DDEFAULT_TOTAL_SECONDS=1740.0 \
	-DDEFAULT_EXCHANGE_SECONDS=300.0 \
	-DDEFAULT_EXCHANGE_MAX_CANDIDATES=250000 \
	-DDEFAULT_EXCHANGE_MAX_REGIONS=100000 \
	-DDEFAULT_EXCHANGE_NODE_LIMIT=3000000LL \
	-DDEFAULT_HIGH_DIVERGENCE_PACK_SECONDS=180.0 \
	-DDEFAULT_EXACT_EXCHANGE_SECONDS=120.0 \
	-DDEFAULT_SMALL_SAVINGS_EXACT_SECONDS=300.0 \
	-DDEFAULT_SMALL_SAVINGS_COMPONENT_LIMIT=400000

.PHONY: all clean

all: solver_submit

$(BUILD):
	mkdir -p $@

$(EMBED): | $(BUILD)
	mkdir -p $@

$(CADICAL_LIB):
	cd $(CADICAL_ROOT) && ./configure
	$(MAKE) -C $(CADICAL_ROOT) -j1

$(BUILD)/manytree_solver: src_manytree/manytree_pool_packer.cpp $(CADICAL_LIB) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(AGGRESSIVE_FLAGS) -DUSE_CADICAL -I$(CADICAL_ROOT)/src $< $(CADICAL_LIB) -o $@

$(HIGHS_LIB):
	cmake -S $(HIGHS_ROOT) -B $(HIGHS_BUILD) -DBUILD_SHARED_LIBS=OFF -DFAST_BUILD=ON -DCMAKE_BUILD_TYPE=Release
	cmake --build $(HIGHS_BUILD) --target highs --parallel 1

$(BUILD)/twotree_solver: $(HIGHS_LIB) | $(BUILD)
	$(MAKE) -C RS h56_main_submit HIGHS_ROOT=$(HIGHS_ROOT) CXXFLAGS="$(CXXFLAGS)"
	cp RS/main_submit $@

$(EMBED)/manytree_solver.o: $(BUILD)/manytree_solver | $(EMBED)
	cp $< $(EMBED)/manytree_solver
	cd $(EMBED) && ld -r -b binary -o manytree_solver.o manytree_solver

$(EMBED)/twotree_solver.o: $(BUILD)/twotree_solver | $(EMBED)
	cp $< $(EMBED)/twotree_solver
	cd $(EMBED) && ld -r -b binary -o twotree_solver.o twotree_solver

solver_submit: combined_wrapper.cpp $(EMBED)/manytree_solver.o $(EMBED)/twotree_solver.o
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -rf $(BUILD)
	$(MAKE) -C RS clean
