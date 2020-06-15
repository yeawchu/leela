#
# Author: Yeaw Chu Lee
# Description: A standard compliant C++ Makefile for single- and multi-target sources
#

#
# Parameter description
#
# MSG		Welcome message
#
# SRC_DIRS	Source directory(s)
# BLD_DIR	Build directory
# BIN_DIR	Binary directory
# DAT_DIR	Data and resource directory
#
# EXE           Single executable
# APP_SRCS	Full path to application source(s) where "int main()" is defined
# TST_SRCS	Full path to test source(s) where "int main()" is defined
# SRC_EXT	Source file extension, e.g. cpp, cc, c, cxx
# OBJ_EXT	Object file extension, e.g. o, obj, object
# DEP_EXT	Dependency file extension, e.g. d, dep, depend
# BIN_EXT	Binary executable extension, e.g. out, exe
# APP_SRC_EXT	Application source file extension, e.g. main.cpp
# TST_SRC_EXT	Test source file extension, e.g. test.cpp
#
# INC_DIRS	Include directory(s), to be prefix by -I flag
# LIB_DIRS	Library directory(s), to be prefix by -L flag
# LIBS		Library source(s), to be prefix by -l flag
#
# PPFLAGS	Pre-processor flags, e.g. -DXXX

#
# The following are the minimum pre-requisite parameters:
#
# 1. SRC_DIRS
# 2. APP_SRCS or APP_SRC_EXT or TST_SRC_EXT
#

SRC_DIRS := leela
SRC_EXT := cpp
DAT_DIR := leela/config
EXE := leela

#
#
#

#
# DO NOT CROSS THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
#

#
#
#

#
# Version
#

VER := 0.1.3

#
# Default values
#

MSG ?= C++ Makefile

SRC_DIRS ?=
BLD_DIR ?= build
BIN_DIR ?= bin
DAT_DIR ?=
EXE ?=

APP_SRCS ?=
TST_SRCS ?=
SRC_EXT ?= cpp
OBJ_EXT ?= o
DEP_EXT ?= d
BIN_EXT ?= out
APP_SRC_EXT ?= main.$(SRC_EXT)
TST_SRC_EXT ?= test.$(SRC_EXT)

INC_DIRS ?=
LIB_DIRS ?=
LIBS ?=
PPFLAGS ?=

#
# Print message
#

$(info Makefile version $(VER))
$(info Copyright (C) 2019 Yeaw Chu Lee)
$(info )
$(info $(MSG))
$(info )

#
# Pre-requisite parameter checks
#

ifeq ($(SRC_DIRS),)
ERR_MSG += SRC_DIRS
endif

ifeq ($(APP_SRCS),)
ifneq ($(ERR_MSG),)
ERR_MSG += and
endif
ERR_MSG += APP_SRCS
endif

ifeq ($(strip $(SRC_DIRS) $(APP_SRCS)),)
$(error $(ERR_MSG) pre-requisite parameter not defined)
endif

#
# Lists
#

COM_SRCS := $(foreach i,$(SRC_DIRS),$(shell if [ -d $(i) ]; then find $(i) -type f -name "*.$(SRC_EXT)" ! -name "*.$(APP_SRC_EXT)" ! -name "*.$(TST_SRC_EXT)" ! -name "$(notdir $(APP_SRCS))" ! -name "$(notdir $(TST_SRCS))"; fi))
APP_SRCS += $(foreach i,$(SRC_DIRS),$(shell if [ -d $(i) ]; then find $(i) -type f -name "*.$(APP_SRC_EXT)"; fi))
TST_SRCS += $(foreach i,$(SRC_DIRS),$(shell if [ -d $(i) ]; then find $(i) -type f -name "*.$(TST_SRC_EXT)"; fi))
DAT_SRCS := $(foreach i,$(DAT_DIR),$(shell if [ -d $(i) ]; then find $(i) -type f -name "*"; fi))
COM_OBJS := $(COM_SRCS:%.$(SRC_EXT)=$(BLD_DIR)/%.$(OBJ_EXT))
APP_OBJS := $(APP_SRCS:%.$(SRC_EXT)=$(BLD_DIR)/%.$(OBJ_EXT))
TST_OBJS := $(TST_SRCS:%.$(SRC_EXT)=$(BLD_DIR)/%.$(OBJ_EXT))
DEPS := $(strip $(COM_OBJS:%.$(OBJ_EXT)=%.$(DEP_EXT)) $(APP_OBJS:%.$(OBJ_EXT)=%.$(DEP_EXT)) $(TST_OBJS:%.$(OBJ_EXT)=%.$(DEP_EXT)))
APPS := $(APP_SRCS:%.$(SRC_EXT)=$(BIN_DIR)/%.$(BIN_EXT)) $(BIN_DIR)/$(EXE)
TSTS := $(TST_SRCS:%.$(SRC_EXT)=$(BIN_DIR)/%.$(BIN_EXT))
DATS := $(DAT_SRCS:$(DAT_DIR)/%=$(BIN_DIR)/%)
CLN_FILES := $(strip $(filter $(COM_OBJS) $(APP_OBJS) $(TST_OBJS) $(DEPS) $(APPS) $(TSTS) $(DATS),$(foreach i,$(strip $(BLD_DIR) $(BIN_DIR)),$(shell if [ -d $(i) ]; then find $(i) -type f -name "*"; fi))))
CLN_DIRS := $(shell for i in $(foreach a,$(strip $(BLD_DIR) $(BIN_DIR)),$(shell if [ -d $(a) ]; then find $(a) -type d -name "*"; fi)); do echo $$i; done | sort -r)

#
# Compilation flags and options
#

CXXFLAGS += -std=c++17 -march=native -O3 -pipe -flto
CXXFLAGS += -Wall -Wextra
CPPFLAGS += $(addprefix -I,$(INC_DIRS)) $(PPFLAGS)
OUTPUT_OPTIONS += -MMD -MP
LDFLAGS += $(addprefix -L,$(LIB_DIRS))
LDLIBS += $(addprefix -l,$(LIBS))

#
# Special built-in targets
#

.DELETE_ON_ERROR:
.SECONDEXPANSION:
.SECONDARY:
.PHONY: all app test data distclean clean help

#
# Rules
#

##Available rules:
##
##all:         Build default app, test and data targets
##app:         Build application targets
##test:        Build test targets
##data:        Copy data or resource files to binary directory
##clean:       Remove all build and binary files
##distclean:   Remove all build and binary files and directories
##help:        Show help options

all: app test data

app: $(APPS)

test: $(TSTS)

data: $(DATS)

$(BIN_DIR)/$(EXE): $(COM_OBJS) | $$(@D)/.
	$(info LINK     $^ -> $@)
	$(LINK.cpp) $(OUTPUT_OPTIONS) $^ $(LDLIBS) -o $@

$(BIN_DIR)/%.$(BIN_EXT): $(COM_OBJS) $(BLD_DIR)/%.$(OBJ_EXT) | $$(@D)/.
	$(info LINK     $^ -> $@)
	$(LINK.cpp) $(OUTPUT_OPTIONS) $^ $(LDLIBS) -o $@

$(BLD_DIR)/%.$(OBJ_EXT): %.$(SRC_EXT) | $$(@D)/.
	$(info COMPILE  $< -> $@)
	@ $(COMPILE.cpp) $(OUTPUT_OPTIONS) -o $@ $<

$(DATS): $$(patsubst $(BIN_DIR)/%,$(DAT_DIR)/%,$$@) | $$(@D)/.
	$(info CP       $< -> $@)
	@ cp $< $@

%/.:
	$(info MKDIR    $(@D))
	@ mkdir -p $(@D)

distclean: clean
ifneq ($(CLN_DIRS),)
	$(info RMDIR    $(CLN_DIRS))
	@ rmdir $(CLN_DIRS)
endif

clean:
ifneq ($(CLN_FILES),)
	$(info RM       $(CLN_FILES))
	@ $(RM) $(CLN_FILES)
endif

help: makefile
	@ sed -n 's/^##//p' $<

include $(wildcard $(DEPS))
