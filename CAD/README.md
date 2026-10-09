# NACA wing generator

A small C++17 program that will generate NACA airfoils, extrude them into a 3D wing,
and export the result for visualisation. It is the geometry ("out-of-CAD") building
block of the AER8875 aero-structural project.

**Current state:** project skeleton only. It builds, runs, and has one test.

## Project layout

```
.
├── CMakeLists.txt          # build description (read it, it's commented)
├── include/nacawing/       # public headers (.hpp): what the code offers
├── src/                    # implementation (.cpp) + main.cpp (the program)
└── tests/                  # small test programs run by CTest
```

Headers *declare* functions, `.cpp` files *define* them. Everything except `main.cpp`
goes into the `nacawing_core` library, so the program, the tests and future tools can
all reuse the same geometry code.

## 1. Install the toolchain

You need a C++17 compiler, CMake 3.16 or newer, and a build tool.

**Windows, option A (easiest):** install Visual Studio Community with the
*Desktop development with C++* workload. It includes the MSVC compiler, CMake and Ninja.
Open this folder with *File > Open > Folder*, then pick `naca_wing.exe` as the startup item.

**Windows, option B (VS Code):**

```powershell
winget install Kitware.CMake
winget install Microsoft.VisualStudio.2022.BuildTools --override "--add Microsoft.VisualStudio.Workload.VCTools --includeRecommended --passive"
```

Then work from *Developer PowerShell for VS*, and in VS Code install the
*C/C++* and *CMake Tools* extensions.

**macOS:**

```bash
xcode-select --install
brew install cmake ninja
```

**Linux (Debian/Ubuntu):**

```bash
sudo apt install build-essential cmake ninja-build gdb
```

(Fedora: `sudo dnf install gcc-c++ cmake ninja-build gdb`.)

**Check:** `cmake --version` and `c++ --version` (on Windows: `cl`).

## 2. Build and run

```bash
cmake -S . -B build        # configure: generate build files into build/
cmake --build build        # compile
./build/naca_wing          # run (Windows: .\build\Debug\naca_wing.exe or .\build\naca_wing.exe)
ctest --test-dir build     # run tests (add -C Debug with Visual Studio)
```

Expected output:

```
NACA wing generator v0.1.0
Toolchain OK: C++17 build is running.
```

For a debug build (single-config generators like Ninja or Makefiles):
`cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug`.

The `build/` folder is disposable: delete it any time and re-run the commands above.

## 3. Roadmap week 6

1. ~~Toolchain and project skeleton~~
2. NACA 4-digit airfoil equations (thickness distribution + camber line), 2D points
3. A basic primitive (a square) to test the geometry and export pipeline
4. Extrude the profile along the span into a 3D wing surface
5. Export (Tecplot `.dat` or `.vtk`) and visualise in ParaView or Tecplot
