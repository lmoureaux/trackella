You need CMake version 3.0 or higher (it's `cmake3` on `lxplus`) and ROOT 6.

Setup build directory:

```
mkdir trackella/build
```

Compile:

```
cd trackella/build
cmake ..
make -j$(nproc)
```

Run:

```
./find_doublets
./print_event_stats
```

If you don't run on the Parallella, you'll have to modify the input file in the
code.
