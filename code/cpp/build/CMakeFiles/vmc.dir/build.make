# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/erikasan/Documents/comphys-2-project1/code/cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/erikasan/Documents/comphys-2-project1/code/cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/vmc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vmc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vmc.dir/flags.make

CMakeFiles/vmc.dir/main/main.o: CMakeFiles/vmc.dir/flags.make
CMakeFiles/vmc.dir/main/main.o: ../main/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/erikasan/Documents/comphys-2-project1/code/cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vmc.dir/main/main.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vmc.dir/main/main.o -c /home/erikasan/Documents/comphys-2-project1/code/cpp/main/main.cpp

CMakeFiles/vmc.dir/main/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vmc.dir/main/main.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erikasan/Documents/comphys-2-project1/code/cpp/main/main.cpp > CMakeFiles/vmc.dir/main/main.i

CMakeFiles/vmc.dir/main/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vmc.dir/main/main.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erikasan/Documents/comphys-2-project1/code/cpp/main/main.cpp -o CMakeFiles/vmc.dir/main/main.s

CMakeFiles/vmc.dir/main/main.o.requires:

.PHONY : CMakeFiles/vmc.dir/main/main.o.requires

CMakeFiles/vmc.dir/main/main.o.provides: CMakeFiles/vmc.dir/main/main.o.requires
	$(MAKE) -f CMakeFiles/vmc.dir/build.make CMakeFiles/vmc.dir/main/main.o.provides.build
.PHONY : CMakeFiles/vmc.dir/main/main.o.provides

CMakeFiles/vmc.dir/main/main.o.provides.build: CMakeFiles/vmc.dir/main/main.o


# Object files for target vmc
vmc_OBJECTS = \
"CMakeFiles/vmc.dir/main/main.o"

# External object files for target vmc
vmc_EXTERNAL_OBJECTS =

vmc: CMakeFiles/vmc.dir/main/main.o
vmc: CMakeFiles/vmc.dir/build.make
vmc: libmylib.a
vmc: CMakeFiles/vmc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikasan/Documents/comphys-2-project1/code/cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vmc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vmc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vmc.dir/build: vmc

.PHONY : CMakeFiles/vmc.dir/build

CMakeFiles/vmc.dir/requires: CMakeFiles/vmc.dir/main/main.o.requires

.PHONY : CMakeFiles/vmc.dir/requires

CMakeFiles/vmc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vmc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vmc.dir/clean

CMakeFiles/vmc.dir/depend:
	cd /home/erikasan/Documents/comphys-2-project1/code/cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikasan/Documents/comphys-2-project1/code/cpp /home/erikasan/Documents/comphys-2-project1/code/cpp /home/erikasan/Documents/comphys-2-project1/code/cpp/build /home/erikasan/Documents/comphys-2-project1/code/cpp/build /home/erikasan/Documents/comphys-2-project1/code/cpp/build/CMakeFiles/vmc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vmc.dir/depend

