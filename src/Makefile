CFLAGS = -O3 -std=c99 -Wall -Wextra -Wpedantic
LDFLAGS = -shared
LIBS = -lm
SRCS = acceleration.c acceleration_barnes_hut.c error.c gravity_sim.c integrator_simple.c integrator_rk_embedded.c integrator_ias15.c integrator_whfast.c math_functions.c storing.c utils.c
OBJS = $(SRCS:.c=.o)

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

.PHONY: all clean

all:
	@echo "Compiling C library..."
	@$(MAKE) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo "Removing everything except source files..."
    ifeq ($(OS),Windows_NT)
		del -f $(OBJS) $(TARGET)
    else
		rm -f $(OBJS) $(TARGET)
    endif