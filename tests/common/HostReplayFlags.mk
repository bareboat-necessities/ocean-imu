# The fixed-draw replay protocol uses AVX2 Eigen packets on Linux x86-64.
# -march=native silently selects AVX-512 on some hosted runners, changing
# floating-point reductions and the long-running nonlinear replay. Keep the
# host regression arithmetic independent of the runner CPU. Firmware flags
# are unaffected. Override HOST_REPLAY_ARCH_FLAGS for a diagnostic build;
# non-default arithmetic is not the fixed-draw acceptance protocol.
HOST_REPLAY_TARGET := $(shell $(CC) -dumpmachine)
ifeq ($(UNAME_S),Linux)
ifneq ($(filter x86_64%,$(HOST_REPLAY_TARGET)),)
    HOST_REPLAY_ARCH_FLAGS ?= -march=x86-64-v3
else
    HOST_REPLAY_ARCH_FLAGS ?= -march=native
endif
else ifeq ($(UNAME_S),Darwin)
    HOST_REPLAY_ARCH_FLAGS ?= -march=native
else
    HOST_REPLAY_ARCH_FLAGS ?=
endif
CXXFLAGS = $(BASEFLAGS) $(HOST_REPLAY_ARCH_FLAGS)
