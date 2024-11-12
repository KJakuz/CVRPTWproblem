CXX = g++
CXXFLAGS = -std=c++20 -Iinclude

# Set the build mode (debug or release)
MODE ?= debug

# Check if MODE is set to either debug or release
check_mode:
ifeq ($(filter $(MODE),debug release),)
	$(error MODE must be set to either debug or release)
else
# Null operation prevents message "Nothing to be done" when MODE is set proprerly
	@: 
endif

# Adjust compiler flags based on the build mode
ifeq ($(MODE),release)
    CXXFLAGS += -O2  # Optimization for release mode
else
    CXXFLAGS += -g  # Debug symbols and warnings for debug mode
endif

# Directories and paths
BUILD_DIR = target
SRC_DIR = src
TARGET_DIR = $(BUILD_DIR)/$(MODE)
OBJ_DIR = $(TARGET_DIR)/object
TARGET = $(TARGET_DIR)/final

# List all source files in the src directory
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)

# Generate a list of object files in the obj directory
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Default target builds the executable
build: check_mode $(TARGET)

# Rule to compile object files from source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link object files into the final executable
$(TARGET): $(OBJ_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) $(LDFLAGS) -o $@


# Clean target removes all generated files
clean:
	rm -rf $(OBJ_DIR) $(BUILD_DIR)

# Compile and Run
run: $(TARGET)
	./$(TARGET) test.txt

# Remove generated files and rebuild the project
rebuild: clean build