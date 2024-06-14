CFLAGS = -shared -O3 -std=gnu99 -Wall -Wextra -Wpedantic

ifeq ($(OS),Windows_NT)
    c_lib.dll: c_lib.c
		@echo "Compiling C library..."
		$(CC) $(CFLAGS) c_lib.c -o c_lib.dll -lm
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CFLAGS += -fPIC -Wmisleading-indentation
        c_lib.so: c_lib.c
			@echo "Compiling C library..."
			$(CC) $(CFLAGS) c_lib.c -o c_lib.so -lm
    else ifeq ($(UNAME_S),Darwin)
        CFLAGS += -fPIC -Wmisleading-indentation
        c_lib.dylib: c_lib.c
			@echo "Compiling C library..."
			$(CC) $(CFLAGS) c_lib.c -o c_lib.dylib -lm
    endif
endif