# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build

# Include any dependencies generated for this target.
include CMakeFiles/assignment.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/assignment.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignment.dir/flags.make

CMakeFiles/assignment.dir/spring.cpp.o: CMakeFiles/assignment.dir/flags.make
CMakeFiles/assignment.dir/spring.cpp.o: ../spring.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/assignment.dir/spring.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/assignment.dir/spring.cpp.o -c /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/spring.cpp

CMakeFiles/assignment.dir/spring.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment.dir/spring.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/spring.cpp > CMakeFiles/assignment.dir/spring.cpp.i

CMakeFiles/assignment.dir/spring.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment.dir/spring.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/spring.cpp -o CMakeFiles/assignment.dir/spring.cpp.s

CMakeFiles/assignment.dir/spring.cpp.o.requires:
.PHONY : CMakeFiles/assignment.dir/spring.cpp.o.requires

CMakeFiles/assignment.dir/spring.cpp.o.provides: CMakeFiles/assignment.dir/spring.cpp.o.requires
	$(MAKE) -f CMakeFiles/assignment.dir/build.make CMakeFiles/assignment.dir/spring.cpp.o.provides.build
.PHONY : CMakeFiles/assignment.dir/spring.cpp.o.provides

CMakeFiles/assignment.dir/spring.cpp.o.provides.build: CMakeFiles/assignment.dir/spring.cpp.o

CMakeFiles/assignment.dir/assignment.cpp.o: CMakeFiles/assignment.dir/flags.make
CMakeFiles/assignment.dir/assignment.cpp.o: ../assignment.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/assignment.dir/assignment.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/assignment.dir/assignment.cpp.o -c /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/assignment.cpp

CMakeFiles/assignment.dir/assignment.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment.dir/assignment.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/assignment.cpp > CMakeFiles/assignment.dir/assignment.cpp.i

CMakeFiles/assignment.dir/assignment.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment.dir/assignment.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/assignment.cpp -o CMakeFiles/assignment.dir/assignment.cpp.s

CMakeFiles/assignment.dir/assignment.cpp.o.requires:
.PHONY : CMakeFiles/assignment.dir/assignment.cpp.o.requires

CMakeFiles/assignment.dir/assignment.cpp.o.provides: CMakeFiles/assignment.dir/assignment.cpp.o.requires
	$(MAKE) -f CMakeFiles/assignment.dir/build.make CMakeFiles/assignment.dir/assignment.cpp.o.provides.build
.PHONY : CMakeFiles/assignment.dir/assignment.cpp.o.provides

CMakeFiles/assignment.dir/assignment.cpp.o.provides.build: CMakeFiles/assignment.dir/assignment.cpp.o

# Object files for target assignment
assignment_OBJECTS = \
"CMakeFiles/assignment.dir/spring.cpp.o" \
"CMakeFiles/assignment.dir/assignment.cpp.o"

# External object files for target assignment
assignment_EXTERNAL_OBJECTS =

bin/assignment: CMakeFiles/assignment.dir/spring.cpp.o
bin/assignment: CMakeFiles/assignment.dir/assignment.cpp.o
bin/assignment: CMakeFiles/assignment.dir/build.make
bin/assignment: /usr/lib/x86_64-linux-gnu/libGLEW.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGL.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGLEW.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libjpeg.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGL.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libjpeg.so
bin/assignment: CMakeFiles/assignment.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/assignment"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignment.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignment.dir/build: bin/assignment
.PHONY : CMakeFiles/assignment.dir/build

CMakeFiles/assignment.dir/requires: CMakeFiles/assignment.dir/spring.cpp.o.requires
CMakeFiles/assignment.dir/requires: CMakeFiles/assignment.dir/assignment.cpp.o.requires
.PHONY : CMakeFiles/assignment.dir/requires

CMakeFiles/assignment.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignment.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignment.dir/clean

CMakeFiles/assignment.dir/depend:
	cd /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build /v/filer4b/v38q001/hari2018/CS354_GRAPHICS/FINAL_PROJECT/build/CMakeFiles/assignment.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assignment.dir/depend

