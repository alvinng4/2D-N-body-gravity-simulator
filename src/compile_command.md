## Compilation
Below is the command I used to compile the source files.
gnu99 instead of c99 is picked because testing on mac shows that gnu99 is faster.

On mac (gcc-13 (Homebrew GCC 13.2.0) 13.2.0):
```
gcc-13 -shared -fPIC -O3 -std=gnu99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.dylib
```

On windows (gcc.exe (Rev5, Built by MSYS2 project) 5.3.0):
```
gcc -shared -O3 -std=gnu99 -Wall -Wextra -Wpedantic c_lib.c -o c_lib.dll
```

On linux (gcc (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0):
```
gcc -shared -fPIC -O3 -std=gnu99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.so
```

After compilation, put the files into the `gravity_plot` or `gravity_sim` folder.
