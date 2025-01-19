CFLAGS = -O3 -std=c99 -Wall -Wextra -Wpedantic
LDFLAGS = -shared
LIBS = -lm
SRCS = tools.c integrator_whfast.c integrator_ias15.c integrator_rk_embedded.c integrator_simple.c acceleration_barnes_hut.c acceleration_fast_multipole.c acceleration.c common.c
OBJS = $(SRCS:.c=.o)

CUDA_FLAGS = -Xcompiler
CUDA_LIBS = -lcudart -lstdc++
CUDA_SRCS = acceleration_cuda.cu
CUDA_OBJS = $(CUDA_SRCS:.cu=.o)

ifeq ($(OS),Windows_NT)
    TARGET = c_lib.dll
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
        TARGET = c_lib.so
        CFLAGS += -fPIC
        CUDA_FLAGS += -fPIC
    else ifeq ($(UNAME_S), Darwin)
        TARGET = c_lib.dylib
        CFLAGS += -fPIC
        CUDA_FLAGS += -fPIC
    endif
endif

ifeq ($(USE_CUDA), 1)
    CFLAGS += -DUSE_CUDA
    OBJS += $(CUDA_OBJS)
    LIBS += $(CUDA_LIBS)
endif

.PHONY: all clean

all:
	@echo "Compiling C library..."
	@$(MAKE) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	nvcc $(CUDA_FLAGS) -c $< -o $@

clean:
	@echo "Removing everything except source files..."
    ifeq ($(OS),Windows_NT)
		del -f $(OBJS) $(CUDA_OBJS) $(TARGET)
    else
		rm -f $(OBJS) $(CUDA_OBJS) $(TARGET)
    endif