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
CMAKE_SOURCE_DIR = /home/erikasan/Documents/comphys-2-project-2/code/cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/erikasan/Documents/comphys-2-project-2/code/cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/vmc_parallel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vmc_parallel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vmc_parallel.dir/flags.make

# Object files for target vmc_parallel
vmc_parallel_OBJECTS =

# External object files for target vmc_parallel
vmc_parallel_EXTERNAL_OBJECTS =

vmc_parallel: CMakeFiles/vmc_parallel.dir/build.make
vmc_parallel: libmylib.a
vmc_parallel: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
vmc_parallel: /usr/lib/x86_64-linux-gnu/libpthread.so
vmc_parallel: CMakeFiles/vmc_parallel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikasan/Documents/comphys-2-project-2/code/cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX executable vmc_parallel"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vmc_parallel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vmc_parallel.dir/build: vmc_parallel

.PHONY : CMakeFiles/vmc_parallel.dir/build

CMakeFiles/vmc_parallel.dir/requires:

.PHONY : CMakeFiles/vmc_parallel.dir/requires

CMakeFiles/vmc_parallel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vmc_parallel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vmc_parallel.dir/clean

CMakeFiles/vmc_parallel.dir/depend:
	cd /home/erikasan/Documents/comphys-2-project-2/code/cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikasan/Documents/comphys-2-project-2/code/cpp /home/erikasan/Documents/comphys-2-project-2/code/cpp /home/erikasan/Documents/comphys-2-project-2/code/cpp/build /home/erikasan/Documents/comphys-2-project-2/code/cpp/build /home/erikasan/Documents/comphys-2-project-2/code/cpp/build/CMakeFiles/vmc_parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vmc_parallel.dir/depend
