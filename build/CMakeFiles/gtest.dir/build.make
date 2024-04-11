# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dmz/Baidu/convex_generation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dmz/Baidu/convex_generation/build

# Utility rule file for gtest.

# Include any custom commands dependencies for this target.
include CMakeFiles/gtest.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/gtest.dir/progress.make

CMakeFiles/gtest: CMakeFiles/gtest-complete

CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-install
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-mkdir
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-download
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-update
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-patch
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-configure
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-build
CMakeFiles/gtest-complete: gtest/src/gtest-stamp/gtest-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'gtest'"
	/usr/bin/cmake -E make_directory /home/dmz/Baidu/convex_generation/build/CMakeFiles
	/usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/CMakeFiles/gtest-complete
	/usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-done

gtest/src/gtest-stamp/gtest-build: gtest/src/gtest-stamp/gtest-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'gtest'"
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && $(MAKE)
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && /usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-build

gtest/src/gtest-stamp/gtest-configure: gtest/tmp/gtest-cfgcmd.txt
gtest/src/gtest-stamp/gtest-configure: gtest/src/gtest-stamp/gtest-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'gtest'"
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && /usr/bin/cmake "-GUnix Makefiles" /home/dmz/Baidu/convex_generation/build/gtest/src/gtest
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && /usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-configure

gtest/src/gtest-stamp/gtest-download: gtest/src/gtest-stamp/download-gtest.cmake
gtest/src/gtest-stamp/gtest-download: gtest/src/gtest-stamp/gtest-urlinfo.txt
gtest/src/gtest-stamp/gtest-download: gtest/src/gtest-stamp/gtest-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (download, verify and extract) for 'gtest'"
	cd /home/dmz/Baidu/convex_generation/build/gtest/src && /usr/bin/cmake -P /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/download-gtest.cmake
	cd /home/dmz/Baidu/convex_generation/build/gtest/src && /usr/bin/cmake -P /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/verify-gtest.cmake
	cd /home/dmz/Baidu/convex_generation/build/gtest/src && /usr/bin/cmake -P /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/extract-gtest.cmake
	cd /home/dmz/Baidu/convex_generation/build/gtest/src && /usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-download

gtest/src/gtest-stamp/gtest-install: gtest/src/gtest-stamp/gtest-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "No install step for 'gtest'"
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && /usr/bin/cmake -E echo_append
	cd /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-build && /usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-install

gtest/src/gtest-stamp/gtest-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'gtest'"
	/usr/bin/cmake -Dcfgdir= -P /home/dmz/Baidu/convex_generation/build/gtest/tmp/gtest-mkdirs.cmake
	/usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-mkdir

gtest/src/gtest-stamp/gtest-patch: gtest/src/gtest-stamp/gtest-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'gtest'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-patch

gtest/src/gtest-stamp/gtest-update: gtest/src/gtest-stamp/gtest-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dmz/Baidu/convex_generation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No update step for 'gtest'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch /home/dmz/Baidu/convex_generation/build/gtest/src/gtest-stamp/gtest-update

gtest: CMakeFiles/gtest
gtest: CMakeFiles/gtest-complete
gtest: gtest/src/gtest-stamp/gtest-build
gtest: gtest/src/gtest-stamp/gtest-configure
gtest: gtest/src/gtest-stamp/gtest-download
gtest: gtest/src/gtest-stamp/gtest-install
gtest: gtest/src/gtest-stamp/gtest-mkdir
gtest: gtest/src/gtest-stamp/gtest-patch
gtest: gtest/src/gtest-stamp/gtest-update
gtest: CMakeFiles/gtest.dir/build.make
.PHONY : gtest

# Rule to build all files generated by this target.
CMakeFiles/gtest.dir/build: gtest
.PHONY : CMakeFiles/gtest.dir/build

CMakeFiles/gtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gtest.dir/clean

CMakeFiles/gtest.dir/depend:
	cd /home/dmz/Baidu/convex_generation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dmz/Baidu/convex_generation /home/dmz/Baidu/convex_generation /home/dmz/Baidu/convex_generation/build /home/dmz/Baidu/convex_generation/build /home/dmz/Baidu/convex_generation/build/CMakeFiles/gtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gtest.dir/depend

