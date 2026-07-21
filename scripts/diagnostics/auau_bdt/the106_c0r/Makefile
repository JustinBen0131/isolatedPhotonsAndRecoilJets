CXX ?= c++
ROOT_CONFIG ?= root-config
THE106_FRAMEWORK_PREFIX ?=
BUILD_DIR ?= build

CXXFLAGS ?= -O2 -g
COMMON_FLAGS := -std=c++11 -Wall -Wextra -Werror -fPIC
PHOOL_CPPFLAGS := -I$(THE106_FRAMEWORK_PREFIX)/include
PHOOL_LDFLAGS := -L$(THE106_FRAMEWORK_PREFIX)/lib -lphool
ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs)

.PHONY: all check clean require-prefix

all: require-prefix $(BUILD_DIR)/libTHE106C0R.so $(BUILD_DIR)/the106_c0r_replay $(BUILD_DIR)/the106_c0r_score_replay $(BUILD_DIR)/the106_c0r_freeze_inventory

require-prefix:
	@test -n "$(THE106_FRAMEWORK_PREFIX)" || \
	  (echo "THE106_FRAMEWORK_PREFIX must name the isolated C0-H2/C0-RH framework prefix" >&2; exit 2)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/libTHE106C0R.so: the106_c0r_diagnostic_writer.cc the106_c0r_diagnostic_writer.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(CXXFLAGS) $(PHOOL_CPPFLAGS) -shared \
	  the106_c0r_diagnostic_writer.cc $(PHOOL_LDFLAGS) -o $@

$(BUILD_DIR)/the106_c0r_replay: the106_c0r_diagnostic_writer.cc the106_c0r_replay.cc the106_c0r_diagnostic_writer.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(CXXFLAGS) $(PHOOL_CPPFLAGS) $(ROOT_CFLAGS) \
	  the106_c0r_diagnostic_writer.cc the106_c0r_replay.cc \
	  $(PHOOL_LDFLAGS) $(ROOT_LIBS) -o $@

$(BUILD_DIR)/the106_c0r_score_replay: the106_c0r_diagnostic_writer.cc the106_c0r_score_replay.cc the106_c0r_diagnostic_writer.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(CXXFLAGS) $(PHOOL_CPPFLAGS) $(ROOT_CFLAGS) \
	  the106_c0r_diagnostic_writer.cc the106_c0r_score_replay.cc \
	  $(PHOOL_LDFLAGS) $(ROOT_LIBS) -lTMVAUtils -lTMVA -o $@

$(BUILD_DIR)/the106_c0r_freeze_inventory: the106_c0r_diagnostic_writer.cc the106_c0r_freeze_inventory.cc the106_c0r_diagnostic_writer.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(CXXFLAGS) $(PHOOL_CPPFLAGS) $(ROOT_CFLAGS) \
	  the106_c0r_diagnostic_writer.cc the106_c0r_freeze_inventory.cc \
	  $(PHOOL_LDFLAGS) $(ROOT_LIBS) -o $@

$(BUILD_DIR)/test_the106_c0r_writer: the106_c0r_diagnostic_writer.cc test_the106_c0r_writer.cc the106_c0r_diagnostic_writer.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(CXXFLAGS) $(PHOOL_CPPFLAGS) \
	  the106_c0r_diagnostic_writer.cc test_the106_c0r_writer.cc \
	  $(PHOOL_LDFLAGS) -o $@

check: require-prefix $(BUILD_DIR)/test_the106_c0r_writer $(BUILD_DIR)/the106_c0r_replay $(BUILD_DIR)/the106_c0r_score_replay
	LD_LIBRARY_PATH="$(THE106_FRAMEWORK_PREFIX)/lib:$${LD_LIBRARY_PATH}" \
	  $(BUILD_DIR)/test_the106_c0r_writer
	LD_LIBRARY_PATH="$(THE106_FRAMEWORK_PREFIX)/lib:$${LD_LIBRARY_PATH}" \
	  python3 test_the106_c0r_replay.py \
	    $(BUILD_DIR)/test_the106_c0r_writer $(BUILD_DIR)/the106_c0r_replay
	LD_LIBRARY_PATH="$(THE106_FRAMEWORK_PREFIX)/lib:$${LD_LIBRARY_PATH}" \
	  python3 test_the106_c0r_score_replay.py \
	    $(BUILD_DIR)/test_the106_c0r_writer $(BUILD_DIR)/the106_c0r_score_replay \
	    "$$($(ROOT_CONFIG) --version)"
	LD_LIBRARY_PATH="$(THE106_FRAMEWORK_PREFIX)/lib:$${LD_LIBRARY_PATH}" \
	  python3 test_the106_c1_identity_order.py \
	    $(BUILD_DIR)/test_the106_c0r_writer $(BUILD_DIR)/the106_c0r_replay

clean:
	rm -rf $(BUILD_DIR)
