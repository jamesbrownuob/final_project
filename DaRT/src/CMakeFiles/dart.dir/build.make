# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /software/spack/linux-rocky8-broadwell/gcc-12.3.0/cmake-3.27.9-s6cv/bin/cmake

# The command to remove a file.
RM = /software/spack/linux-rocky8-broadwell/gcc-12.3.0/cmake-3.27.9-s6cv/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /user/home/hm20670/DaRT_RBE/DaRT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /user/home/hm20670/DaRT_RBE/DaRT/src

# Include any dependencies generated for this target.
include CMakeFiles/dart.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/dart.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/dart.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dart.dir/flags.make

CMakeFiles/dart.dir/dart.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/dart.cc.o: /user/home/hm20670/DaRT_RBE/DaRT/dart.cc
CMakeFiles/dart.dir/dart.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dart.dir/dart.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/dart.cc.o -MF CMakeFiles/dart.dir/dart.cc.o.d -o CMakeFiles/dart.dir/dart.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/dart.cc

CMakeFiles/dart.dir/dart.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/dart.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/dart.cc > CMakeFiles/dart.dir/dart.cc.i

CMakeFiles/dart.dir/dart.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/dart.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/dart.cc -o CMakeFiles/dart.dir/dart.cc.s

CMakeFiles/dart.dir/ActionInitialization.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/ActionInitialization.cc.o: ActionInitialization.cc
CMakeFiles/dart.dir/ActionInitialization.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/dart.dir/ActionInitialization.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/ActionInitialization.cc.o -MF CMakeFiles/dart.dir/ActionInitialization.cc.o.d -o CMakeFiles/dart.dir/ActionInitialization.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/ActionInitialization.cc

CMakeFiles/dart.dir/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/ActionInitialization.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/ActionInitialization.cc > CMakeFiles/dart.dir/ActionInitialization.cc.i

CMakeFiles/dart.dir/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/ActionInitialization.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/ActionInitialization.cc -o CMakeFiles/dart.dir/ActionInitialization.cc.s

CMakeFiles/dart.dir/CommandLineParser.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/CommandLineParser.cc.o: CommandLineParser.cc
CMakeFiles/dart.dir/CommandLineParser.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/dart.dir/CommandLineParser.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/CommandLineParser.cc.o -MF CMakeFiles/dart.dir/CommandLineParser.cc.o.d -o CMakeFiles/dart.dir/CommandLineParser.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/CommandLineParser.cc

CMakeFiles/dart.dir/CommandLineParser.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/CommandLineParser.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/CommandLineParser.cc > CMakeFiles/dart.dir/CommandLineParser.cc.i

CMakeFiles/dart.dir/CommandLineParser.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/CommandLineParser.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/CommandLineParser.cc -o CMakeFiles/dart.dir/CommandLineParser.cc.s

CMakeFiles/dart.dir/DetectorConstruction.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/DetectorConstruction.cc.o: DetectorConstruction.cc
CMakeFiles/dart.dir/DetectorConstruction.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/dart.dir/DetectorConstruction.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/DetectorConstruction.cc.o -MF CMakeFiles/dart.dir/DetectorConstruction.cc.o.d -o CMakeFiles/dart.dir/DetectorConstruction.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorConstruction.cc

CMakeFiles/dart.dir/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/DetectorConstruction.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorConstruction.cc > CMakeFiles/dart.dir/DetectorConstruction.cc.i

CMakeFiles/dart.dir/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/DetectorConstruction.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorConstruction.cc -o CMakeFiles/dart.dir/DetectorConstruction.cc.s

CMakeFiles/dart.dir/DetectorMessenger.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/DetectorMessenger.cc.o: DetectorMessenger.cc
CMakeFiles/dart.dir/DetectorMessenger.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/dart.dir/DetectorMessenger.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/DetectorMessenger.cc.o -MF CMakeFiles/dart.dir/DetectorMessenger.cc.o.d -o CMakeFiles/dart.dir/DetectorMessenger.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorMessenger.cc

CMakeFiles/dart.dir/DetectorMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/DetectorMessenger.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorMessenger.cc > CMakeFiles/dart.dir/DetectorMessenger.cc.i

CMakeFiles/dart.dir/DetectorMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/DetectorMessenger.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/DetectorMessenger.cc -o CMakeFiles/dart.dir/DetectorMessenger.cc.s

CMakeFiles/dart.dir/EventAction.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/EventAction.cc.o: EventAction.cc
CMakeFiles/dart.dir/EventAction.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/dart.dir/EventAction.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/EventAction.cc.o -MF CMakeFiles/dart.dir/EventAction.cc.o.d -o CMakeFiles/dart.dir/EventAction.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/EventAction.cc

CMakeFiles/dart.dir/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/EventAction.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/EventAction.cc > CMakeFiles/dart.dir/EventAction.cc.i

CMakeFiles/dart.dir/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/EventAction.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/EventAction.cc -o CMakeFiles/dart.dir/EventAction.cc.s

CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o: G4RadioactiveDecay.cc
CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o -MF CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o.d -o CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecay.cc

CMakeFiles/dart.dir/G4RadioactiveDecay.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/G4RadioactiveDecay.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecay.cc > CMakeFiles/dart.dir/G4RadioactiveDecay.cc.i

CMakeFiles/dart.dir/G4RadioactiveDecay.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/G4RadioactiveDecay.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecay.cc -o CMakeFiles/dart.dir/G4RadioactiveDecay.cc.s

CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o: G4RadioactiveDecayPhysics.cc
CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o -MF CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o.d -o CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecayPhysics.cc

CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecayPhysics.cc > CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.i

CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/G4RadioactiveDecayPhysics.cc -o CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.s

CMakeFiles/dart.dir/PhysicsList.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/PhysicsList.cc.o: PhysicsList.cc
CMakeFiles/dart.dir/PhysicsList.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/dart.dir/PhysicsList.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/PhysicsList.cc.o -MF CMakeFiles/dart.dir/PhysicsList.cc.o.d -o CMakeFiles/dart.dir/PhysicsList.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/PhysicsList.cc

CMakeFiles/dart.dir/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/PhysicsList.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/PhysicsList.cc > CMakeFiles/dart.dir/PhysicsList.cc.i

CMakeFiles/dart.dir/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/PhysicsList.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/PhysicsList.cc -o CMakeFiles/dart.dir/PhysicsList.cc.s

CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o: PrimaryGeneratorAction.cc
CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o -MF CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o.d -o CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/PrimaryGeneratorAction.cc

CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/PrimaryGeneratorAction.cc > CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.i

CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/PrimaryGeneratorAction.cc -o CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.s

CMakeFiles/dart.dir/RunAction.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/RunAction.cc.o: RunAction.cc
CMakeFiles/dart.dir/RunAction.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/dart.dir/RunAction.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/RunAction.cc.o -MF CMakeFiles/dart.dir/RunAction.cc.o.d -o CMakeFiles/dart.dir/RunAction.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/RunAction.cc

CMakeFiles/dart.dir/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/RunAction.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/RunAction.cc > CMakeFiles/dart.dir/RunAction.cc.i

CMakeFiles/dart.dir/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/RunAction.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/RunAction.cc -o CMakeFiles/dart.dir/RunAction.cc.s

CMakeFiles/dart.dir/SteppingAction.cc.o: CMakeFiles/dart.dir/flags.make
CMakeFiles/dart.dir/SteppingAction.cc.o: SteppingAction.cc
CMakeFiles/dart.dir/SteppingAction.cc.o: CMakeFiles/dart.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/dart.dir/SteppingAction.cc.o"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dart.dir/SteppingAction.cc.o -MF CMakeFiles/dart.dir/SteppingAction.cc.o.d -o CMakeFiles/dart.dir/SteppingAction.cc.o -c /user/home/hm20670/DaRT_RBE/DaRT/src/SteppingAction.cc

CMakeFiles/dart.dir/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dart.dir/SteppingAction.cc.i"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /user/home/hm20670/DaRT_RBE/DaRT/src/SteppingAction.cc > CMakeFiles/dart.dir/SteppingAction.cc.i

CMakeFiles/dart.dir/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dart.dir/SteppingAction.cc.s"
	/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-sknc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /user/home/hm20670/DaRT_RBE/DaRT/src/SteppingAction.cc -o CMakeFiles/dart.dir/SteppingAction.cc.s

# Object files for target dart
dart_OBJECTS = \
"CMakeFiles/dart.dir/dart.cc.o" \
"CMakeFiles/dart.dir/ActionInitialization.cc.o" \
"CMakeFiles/dart.dir/CommandLineParser.cc.o" \
"CMakeFiles/dart.dir/DetectorConstruction.cc.o" \
"CMakeFiles/dart.dir/DetectorMessenger.cc.o" \
"CMakeFiles/dart.dir/EventAction.cc.o" \
"CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o" \
"CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o" \
"CMakeFiles/dart.dir/PhysicsList.cc.o" \
"CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/dart.dir/RunAction.cc.o" \
"CMakeFiles/dart.dir/SteppingAction.cc.o"

# External object files for target dart
dart_EXTERNAL_OBJECTS =

dart: CMakeFiles/dart.dir/dart.cc.o
dart: CMakeFiles/dart.dir/ActionInitialization.cc.o
dart: CMakeFiles/dart.dir/CommandLineParser.cc.o
dart: CMakeFiles/dart.dir/DetectorConstruction.cc.o
dart: CMakeFiles/dart.dir/DetectorMessenger.cc.o
dart: CMakeFiles/dart.dir/EventAction.cc.o
dart: CMakeFiles/dart.dir/G4RadioactiveDecay.cc.o
dart: CMakeFiles/dart.dir/G4RadioactiveDecayPhysics.cc.o
dart: CMakeFiles/dart.dir/PhysicsList.cc.o
dart: CMakeFiles/dart.dir/PrimaryGeneratorAction.cc.o
dart: CMakeFiles/dart.dir/RunAction.cc.o
dart: CMakeFiles/dart.dir/SteppingAction.cc.o
dart: CMakeFiles/dart.dir/build.make
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4Tree.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4FR.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4GMocren.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4visHepRep.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4RayTracer.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4VRML.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4ToolsSG.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4vis_management.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4modeling.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4interfaces.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4persistency.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4error_propagation.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4readout.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4physicslists.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4run.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4event.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4tracking.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4parmodels.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4processes.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4digits_hits.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4track.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4particles.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4geometry.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4materials.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4graphics_reps.so
dart: libgit_version.a
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/xerces-c-3.2.5-kv2e/lib/libxerces-c.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4analysis.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/zlib-ng-2.1.6-reu4/lib/libz.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/expat-2.6.2-aaht/lib/libexpat.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4intercoms.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4global.so
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/geant4-11.1.3-wdly/lib64/libG4ptl.so.2.3.3
dart: /software/spack/linux-rocky8-broadwell/gcc-12.3.0/clhep-2.4.7.1-dgls/lib/libCLHEP-2.4.7.1.so
dart: CMakeFiles/dart.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable dart"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dart.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dart.dir/build: dart
.PHONY : CMakeFiles/dart.dir/build

CMakeFiles/dart.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dart.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dart.dir/clean

CMakeFiles/dart.dir/depend:
	cd /user/home/hm20670/DaRT_RBE/DaRT/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /user/home/hm20670/DaRT_RBE/DaRT /user/home/hm20670/DaRT_RBE/DaRT /user/home/hm20670/DaRT_RBE/DaRT/src /user/home/hm20670/DaRT_RBE/DaRT/src /user/home/hm20670/DaRT_RBE/DaRT/src/CMakeFiles/dart.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/dart.dir/depend

