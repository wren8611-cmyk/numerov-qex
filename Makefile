# ============================================================
# User configuration
# ============================================================
include config.mk

# ============================================================
# Python configuration
# ============================================================
PYTHON = python3
PYTHON_INCLUDE := $(shell $(PYTHON) -c \
    "import sysconfig; print(sysconfig.get_paths()['include'])")
PYBIND_DIR := $(shell $(PYTHON) -c \
    "import pybind11; print(pybind11.get_include())")
EXT_SUFFIX := $(shell $(PYTHON) -c \
    "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

# ============================================================
# Compiler configuration
# ============================================================

CXXFLAGS = -O3 -std=c++17 -fPIC

INCLUDES = \
    -I$(PYTHON_INCLUDE) \
    -I$(PYBIND_DIR) \
    -I$(EIGEN_DIR) \
    -I$(GSL_DIR)/include \
    -Iinclude

# ============================================================
# Linker configuration
# ============================================================

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
else
    LDFLAGS = -shared
endif

LIBS = \
    -L$(GSL_DIR)/lib \
    -lgsl \
    -lgslcblas

# ============================================================
# Build configuration
# ============================================================

TARGET = numerov_qex/radialsolver$(EXT_SUFFIX)

SOURCES = src/radialsolver.cxx
HEADERS = $(wildcard include/*.hxx)

# ============================================================
# Targets
# ============================================================

.PHONY: all build clean rebuild

all: build

build: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) \
	    $(SOURCES) \
	    $(LDFLAGS) $(LIBS) \
	    -o $@

clean:
	rm -f numerov_qex/radialsolver*.so

rebuild:
	$(MAKE) clean
	$(MAKE) build