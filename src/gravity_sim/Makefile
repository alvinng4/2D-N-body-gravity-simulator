CFLAGS = -shared -fPIC -O3 -std=gnu99 -Wall -Wextra -Wpedantic -Wmisleading-indentation

ifeq ($(OS),Windows_NT)
    c_lib.dll: c_lib.c
		@echo "Compiling C library..."
		$(CC) $(CFLAGS) c_lib.c -o c_lib.dll
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        c_lib.so: c_lib.c
			@echo "Compiling C library..."
			$(CC) $(CFLAGS) c_lib.c -o c_lib.so
    else ifeq ($(UNAME_S),Darwin)
        c_lib.dylib: c_lib.c
			@echo "Compiling C library..."
			$(CC) $(CFLAGS) c_lib.c -o c_lib.dylib
    endif
endif