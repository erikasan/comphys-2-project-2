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
include CMakeFiles/gradientdescent.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/gradientdescent.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gradientdescent.dir/flags.make

# Object files for target gradientdescent
gradientdescent_OBJECTS =

# External object files for target gradientdescent
gradientdescent_EXTERNAL_OBJECTS =

gradientdescent: CMakeFiles/gradientdescent.dir/build.make
gradientdescent: libmylib.a
gradientdescent: CMakeFiles/gradientdescent.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikasan/Documents/comphys-2-project-2/code/cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX executable gradientdescent"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gradientdescent.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gradientdescent.dir/build: gradientdescent

.PHONY : CMakeFiles/gradientdescent.dir/build

CMakeFiles/gradientdescent.dir/requires:

.PHONY : CMakeFiles/gradientdescent.dir/requires

CMakeFiles/gradientdescent.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gradientdescent.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gradientdescent.dir/clean

CMakeFiles/gradientdescent.dir/depend:
	cd /home/erikasan/Documents/comphys-2-project-2/code/cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikasan/Documents/comphys-2-project-2/code/cpp /home/erikasan/Documents/comphys-2-project-2/code/cpp /home/erikasan/Documents/comphys-2-project-2/code/cpp/build /home/erikasan/Documents/comphys-2-project-2/code/cpp/build /home/erikasan/Documents/comphys-2-project-2/code/cpp/build/CMakeFiles/gradientdescent.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gradientdescent.dir/depend
