BUILD_DIR := build

CC := gcc

SRCS := qam.c qamDecoder.c
OBJS := $(addprefix $(BUILD_DIR)/,$(SRCS:.c=.o))
DEPS := $(OBJS:.o:.d)

DEPFLAGS = -MT "$@" -MMD -MP -MF "$(BUILD_DIR)/$*.d"

CFLAGS := \
	#-O3 \
	-Wall \
	-Wextra \
	-Wstrict-prototypes \
	-flto \
	-g3 \
	-std=gnu18

LDFLAGS := \
	-lm

.PHONY: all qam qamDecoder clean
all: qam qamDecoder

qam: $(BUILD_DIR)/qam

qamDecoder: $(BUILD_DIR)/qamDecoder

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	@echo [ CC ] $@
	@$(CC) -x c $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/qam: $(BUILD_DIR)/qam.o
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/qamDecoder: $(BUILD_DIR)/qamDecoder.o
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)

$(DEPS): ;
include $(wildcard $(DEPS))
