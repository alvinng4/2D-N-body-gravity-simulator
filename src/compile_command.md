## Compilation
To compile the C libraries, you can run
```
make [CC="your C compiler"]
```
where `CC="your C compiler"` is optional.
After compilation, move the files into the `gravity_plot` or `gravity_sim` folder.

## Specifications on the given libraries files
In case you are interested, below are the commands I used to compile the libraries.
gnu99 instead of c99 is picked because testing on mac shows that gnu99 is faster.

On mac (gcc-14 (Homebrew GCC 14.1.0) 14.1.0):
```
gcc-14 -shared -fPIC -O3 -std=gnu99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.dylib
```

On windows (gcc.exe (Rev5, Built by MSYS2 project) 5.3.0):
```
gcc -shared -O3 -std=gnu99 -Wall -Wextra -Wpedantic c_lib.c -o c_lib.dll
```

On linux (gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0):
```
gcc -shared -fPIC -O3 -std=gnu99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.so
```