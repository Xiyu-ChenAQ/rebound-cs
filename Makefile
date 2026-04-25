# Cosmic Stars — unified Makefile
#
# Targets:
#   make              - build librebound + libcs
#   make librebound   - build REBOUND shared library only
#   make libcs        - build CS static library
#   make test         - build and run test suite
#   make demo         - build and run physics demos
#   make clean        - remove all build artifacts
#
# Requirements:
#   Linux/macOS : gcc + ar + make
#   Windows     : MSVC (run vcvars64.bat first) + make

include src/Makefile.defs

# --- directories ---
SRCDIR = src
CSDIR  = cs
BLDDIR = build

# --- source files ---
CS_SRC   = $(CSDIR)/cs_simulation.c $(CSDIR)/cs_gr.c \
           $(CSDIR)/cs_radiation.c $(CSDIR)/cs_harmonics.c \
           $(CSDIR)/cs_tides.c $(CSDIR)/cs_solarmass.c
SRC_ALL  = $(wildcard $(SRCDIR)/*.c)
TST_SRC  = $(BLDDIR)/test_cs.c
DEMO_SRC = $(BLDDIR)/demo_physics.c
DMAS_SRC = $(BLDDIR)/demo_mass.c

# --- object files (all go to build/) ---
CS_OBJ   = $(patsubst $(CSDIR)/%.c,$(BLDDIR)/%.$(OBJFILEEXT),$(CS_SRC))
SRC_OBJ  = $(patsubst $(SRCDIR)/%.c,$(BLDDIR)/%.$(OBJFILEEXT),$(SRC_ALL))
TST_OBJ  = $(patsubst $(BLDDIR)/%.c,$(BLDDIR)/%.$(OBJFILEEXT),$(TST_SRC))
DEMO_OBJ = $(patsubst $(BLDDIR)/%.c,$(BLDDIR)/%.$(OBJFILEEXT),$(DEMO_SRC))
DMAS_OBJ = $(patsubst $(BLDDIR)/%.c,$(BLDDIR)/%.$(OBJFILEEXT),$(DMAS_SRC))

# --- library / executable names ---
ifeq ($(OS), Windows_NT)
    LIBCS   = $(BLDDIR)/cs.lib
    TST_EXE = $(BLDDIR)/test_cs.exe
    DEMO_EXE= $(BLDDIR)/demo_physics.exe
    DMAS_EXE= $(BLDDIR)/demo_mass.exe
else
    LIBCS   = $(BLDDIR)/libcs.a
    TST_EXE = $(BLDDIR)/test_cs
    DEMO_EXE= $(BLDDIR)/demo_physics
    DMAS_EXE= $(BLDDIR)/demo_mass
endif

CINCL = -I. -I$(SRCDIR)

# ====================================================================
# default target
# ====================================================================
.PHONY: all librebound libcs test demo demo_mass clean

all: librebound libcs

# ====================================================================
# REBOUND shared library (delegates to src/)
# ====================================================================
librebound:
	$(MAKE) -C $(SRCDIR)
	@$(LINKORCOPYLIBREBOUNDMAIN)

# ====================================================================
# CS static library
# ====================================================================
libcs: $(LIBCS)

ifeq ($(OS), Windows_NT)
$(LIBCS): $(CS_OBJ)
	@echo "  LIB $@"
	lib /nologo /OUT:$@ $(CS_OBJ)
else
$(LIBCS): $(CS_OBJ)
	@echo "  AR  $@"
	ar rcs $@ $(CS_OBJ)
endif

# ====================================================================
# test suite
# ====================================================================
test: $(SRC_OBJ) $(LIBCS) $(TST_EXE)
	@echo "  RUN  $(notdir $(TST_EXE))"
	cd $(BLDDIR) && ./$(notdir $(TST_EXE))

$(TST_EXE): $(TST_OBJ) $(CS_OBJ) $(SRC_OBJ)
ifeq ($(OS), Windows_NT)
	@echo "  LINK $@"
	$(CC) /nologo $(SRC_OBJ) $(CS_OBJ) $(TST_OBJ) /Fe:$@
else
	@echo "  LINK $@"
	$(CC) $(OPT) $(SRC_OBJ) $(CS_OBJ) $(TST_OBJ) $(LIB) -o $@
endif

# ====================================================================
# physics demos
# ====================================================================
demo: $(SRC_OBJ) $(LIBCS) $(DEMO_EXE)
	@echo "  RUN  $(notdir $(DEMO_EXE))"
	cd $(BLDDIR) && ./$(notdir $(DEMO_EXE))

$(DEMO_EXE): $(DEMO_OBJ) $(CS_OBJ) $(SRC_OBJ)
ifeq ($(OS), Windows_NT)
	@echo "  LINK $@"
	$(CC) /nologo $(SRC_OBJ) $(CS_OBJ) $(DEMO_OBJ) /Fe:$@
else
	@echo "  LINK $@"
	$(CC) $(OPT) $(SRC_OBJ) $(CS_OBJ) $(DEMO_OBJ) $(LIB) -o $@
endif

demo_mass: $(SRC_OBJ) $(LIBCS) $(DMAS_EXE)
	@echo "  RUN  $(notdir $(DMAS_EXE))"
	cd $(BLDDIR) && ./$(notdir $(DMAS_EXE))

$(DMAS_EXE): $(DMAS_OBJ) $(CS_OBJ) $(SRC_OBJ)
ifeq ($(OS), Windows_NT)
	@echo "  LINK $@"
	$(CC) /nologo $(SRC_OBJ) $(CS_OBJ) $(DMAS_OBJ) /Fe:$@
else
	@echo "  LINK $@"
	$(CC) $(OPT) $(SRC_OBJ) $(CS_OBJ) $(DMAS_OBJ) $(LIB) -o $@
endif

# ====================================================================
# compile pattern rules (output → build/)
# ====================================================================
ifeq ($(OS), Windows_NT)
$(BLDDIR)/%.$(OBJFILEEXT): $(CSDIR)/%.c
	@if not exist "$(BLDDIR)" mkdir "$(BLDDIR)"
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) $< -Fo$@

$(BLDDIR)/%.$(OBJFILEEXT): $(SRCDIR)/%.c
	@if not exist "$(BLDDIR)" mkdir "$(BLDDIR)"
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) $< -Fo$@

$(BLDDIR)/%.$(OBJFILEEXT): $(BLDDIR)/%.c
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) $< -Fo$@
else
$(BLDDIR)/%.$(OBJFILEEXT): $(CSDIR)/%.c
	@mkdir -p $(BLDDIR)
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) -o $@ $<

$(BLDDIR)/%.$(OBJFILEEXT): $(SRCDIR)/%.c
	@mkdir -p $(BLDDIR)
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) -o $@ $<

$(BLDDIR)/%.$(OBJFILEEXT): $(BLDDIR)/%.c
	@echo "  CC  $<"
	$(CC) -c $(OPT) $(PREDEF) $(CINCL) -o $@ $<
endif

# ====================================================================
# clean
# ====================================================================
clean:
	@echo "Cleaning build artifacts..."
	-$(RM) $(BLDDIR)/*.$(OBJFILEEXT)
	-$(RM) $(LIBCS)
	-$(RM) $(TST_EXE) $(DEMO_EXE) $(DMAS_EXE)
	-$(RM) $(BLDDIR)/*.lib $(BLDDIR)/*.exp
	$(MAKE) -C $(SRCDIR) clean
