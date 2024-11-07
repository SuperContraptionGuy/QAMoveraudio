BUILD_DIR := build
SRC_DIR := src

TARGETS := \
	$(BUILD_DIR)/qam \
	$(BUILD_DIR)/qamDecoder

CC := gcc

SRCS := \
	$(SRC_DIR)/avg.c \
	$(SRC_DIR)/avg_complex.c \
	$(SRC_DIR)/pid.c

OBJS := $(addprefix $(BUILD_DIR)/,$(notdir $(SRCS:.c=.o)))
DEPS := $(OBJS:.o:.d) $(TARGET_DEPS)

DEPFLAGS = -MT "$@" -MMD -MP -MF "$(BUILD_DIR)/$*.d"

CFLAGS := \
	-Wall \
	-Wextra \
	-Wstrict-prototypes \
	-flto \
	-g3 \
	-std=gnu18 \
	-I$(SRC_DIR) \
	# -O3

LDFLAGS := \
	-lm

.PHONY: all clean
all: $(TARGETS)

$(TARGETS): $(TARGET_OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo [ CC ] $@
	@$(CC) -x c $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/qam: $(BUILD_DIR)/qam.o $(OBJS)
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/qamDecoder: $(BUILD_DIR)/qamDecoder.o $(OBJS)
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)

$(DEPS): ;
include $(wildcard $(DEPS))
