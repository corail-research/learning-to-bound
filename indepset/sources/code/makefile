# --- SYSTEM ---

CCC = g++

# --- FLAGS ---

CCOPT = -m64 -fPIC -fexceptions -DIL_STD

# --- OPTIMIZATION FLAGS ---

DEBUG_OPT = -DNDEBUG -O3

L_DIR=../../../learning-DD/models/misp-random/code
GNN_HOME=../../../learning-DD/graphnn
include $(GNN_HOME)/make_common
include_dirs = $(CUDA_HOME)/include $(MKL_ROOT)/include $(GNN_HOME)/include $(L_DIR)/include/dd $(L_DIR)/include/learning ./src

lib_dir = $(GNN_HOME)/build_cpuonly/lib
lib_dir2 = $(L_DIR)/build/lib
lib_dir3 = $(MKL_ROOT)/lib/intel64_lin
gnn_lib = $(lib_dir)/libgnn.a

CFLAGS = $(addprefix -I,$(include_dirs)) -ansi -g -std=c++11 -pedantic -Wall $(DEBUG_OPT) -c

CXXFLAGS += $(addprefix -I,$(include_dirs)) -g -Wno-unused-local-typedef -fPIC

LDFLAGS = -lm -lmkl_rt -ltbb 
#LDFLAGS = -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ltbb 

# ---- COMPILE  ----
SRC_DI    := src
SRC_DIR   := src $(L_DIR)/src/dd $(L_DIR)/src/learning
OBJ_DIR   := obj

SRC_DIRS  := $(shell find $(SRC_DIR) -type d)
OBJ_DIRS  := $(addprefix $(OBJ_DIR)/,$(SRC_DIRS))

SOURCES   := $(shell find $(SRC_DI) -name '*.cpp')
OBJ_FILES := $(addprefix $(OBJ_DIR)/, $(SOURCES:.cpp=.o))
OBJ_FILES2 := $(shell find $(lib_dir2) -name '*.o')

cpp_files = $(shell $(FIND) $(L_DIR)/src/learning -name "*.cpp" -printf "%P\n")
cpp_files += $(shell $(FIND) $(L_DIR)/src/dd -name "*.cpp" -printf "%P\n")

cxx_obj_files = $(subst .cpp,.o,$(cpp_files))

objs = $(addprefix $(L_DIR)/build/lib/,$(cxx_obj_files))
DEPS = $(objs:.o=.d)

vpath %.cpp $(SRC_DIRS)

# ---- TARGETS ----

EXECUTABLE=bb
LSO=$(L_DIR)/build/dll
LL=$(LSO)/learning_lib.so

all: $(LL) $(EXECUTABLE)

$(LL) : $(L_DIR)/src/learning_lib.cpp $(gnn_lib) $(objs)
	@echo LLSO
	$(dir_guard)
	$(CXX) -shared $(CXXFLAGS) -MMD -o $@ $(filter %.cpp %.o, $^) -L$(lib_dir) -lgnn $(LDFLAGS)

$(EXECUTABLE): makedir $(OBJ_FILES) $(gnn_lib) $(objs)
	$(CCC) $(OBJ_FILES) $(gnn_lib) $(LL) $(objs) $(LDFLAGS) $(CXXFLAGS) -o $@ -L$(lib_dir) -L$(LSO)
#$(L_DIR)/build/lib/
$(OBJ_DIR)/%.o: %.cpp 
	$(CCC) $(CFLAGS) -MMD $< -o $@ -L$(lib_dir)

makedir: $(OBJ_DIRS)

DEPS += $(target_dep)


$(L_DIR)/build/lib/%.o: $(L_DIR)/src/learning/%.cpp
	@echo LLO
	$(dir_guard)
	$(CXX) $(CXXFLAGS) -MMD -c -o $@ $(filter %.cpp, $^)

$(L_DIR)/build/lib/%.o: $(L_DIR)/src/dd/%.cpp
	$(dir_guard)
	$(CXX) $(CXXFLAGS) -MMD -c -o $@ $(filter %.cpp, $^)

$(OBJ_DIRS):
	@mkdir -p $@


-include $(DEPS)


clean:
	@rm -rf obj 
	@rm -rf build
	@rm -rf $(EXECUTABLE)
