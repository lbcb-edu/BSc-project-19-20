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
CMAKE_SOURCE_DIR = /home/ivona/FER/Projekt/pink

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivona/FER/Projekt/pink

# Include any dependencies generated for this target.
include CMakeFiles/pink_mapper.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pink_mapper.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pink_mapper.dir/flags.make

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o: CMakeFiles/pink_mapper.dir/flags.make
CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o: src/pink_mapper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivona/FER/Projekt/pink/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o -c /home/ivona/FER/Projekt/pink/src/pink_mapper.cpp

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivona/FER/Projekt/pink/src/pink_mapper.cpp > CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.i

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivona/FER/Projekt/pink/src/pink_mapper.cpp -o CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.s

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.requires:

.PHONY : CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.requires

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.provides: CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.requires
	$(MAKE) -f CMakeFiles/pink_mapper.dir/build.make CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.provides.build
.PHONY : CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.provides

CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.provides.build: CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o


# Object files for target pink_mapper
pink_mapper_OBJECTS = \
"CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o"

# External object files for target pink_mapper
pink_mapper_EXTERNAL_OBJECTS =

pink_mapper: CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o
pink_mapper: CMakeFiles/pink_mapper.dir/build.make
pink_mapper: libpink_alignment.a
pink_mapper: vendor/bioparser/lib/libz.a
pink_mapper: CMakeFiles/pink_mapper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivona/FER/Projekt/pink/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable pink_mapper"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pink_mapper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pink_mapper.dir/build: pink_mapper

.PHONY : CMakeFiles/pink_mapper.dir/build

CMakeFiles/pink_mapper.dir/requires: CMakeFiles/pink_mapper.dir/src/pink_mapper.cpp.o.requires

.PHONY : CMakeFiles/pink_mapper.dir/requires

CMakeFiles/pink_mapper.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pink_mapper.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pink_mapper.dir/clean

CMakeFiles/pink_mapper.dir/depend:
	cd /home/ivona/FER/Projekt/pink && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivona/FER/Projekt/pink /home/ivona/FER/Projekt/pink /home/ivona/FER/Projekt/pink /home/ivona/FER/Projekt/pink /home/ivona/FER/Projekt/pink/CMakeFiles/pink_mapper.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pink_mapper.dir/depend

