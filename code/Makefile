#@file    Makefile
#@brief   Makefile for compare cutting planes project
#@author  Christopher Hojny


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       ${SCIP_DEVELOP_PATH}

#----------------------------------------------------------------------------
# include default project Makefile from SCIP
#----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
# -include $(SCIPDIR)/make/local/make.moskito
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		=	1.0
CONTINUE	=	false
LOCK		= 	false
TEST    	=       mipdev2-solvable
SETTINGS    	=   	default
TIME		=	3600
NODE		=	2100000000
MEM 		=	20000
DISPFREQ	=	10000
PERMUTE         =       0
ONLYPRE		=       false
SETCUTOFF	=	false
SEEDS       =   0
SETUSESOL   =   false

# determine git hash
GITDEF 		:= 	-D'SYMGITHASH="$(shell git describe --always | sed 's/^.*-g//')"'

USRCFLAGS	+=      -Iinclude -Ilib $(GITDEF)
USRCXXFLAGS	+=      -Iinclude -Ilib $(GITDEF)

#-----------------------------------------------------------------------------
# Programs
#-----------------------------------------------------------------------------

SBCS		=	sbcs
SBCSOBJ		=	main.o permutation.o prop_symmetry.o readArguments.o cons_symretope.o cons_orbisack.o cons_orbitope.o cons_symresack.o
SBCSSRC		=	main.cpp permutation.c prop_symmetry.c readArguments.c cons_symretope.c cons_orbisack.c cons_orbitope.c cons_symresack.c
SBCSOBJFILES	=	$(addprefix $(OBJDIR)/,$(SBCSOBJ))
SBCSSRCFILES	=	$(addprefix $(SRCDIR)/,$(SBCSSRC))
SBCSDEP		=	$(SRCDIR)/depend.$(SBCS)
SBCSFILE	=	$(BINDIR)/$(SBCS).$(BASE).$(LPS)$(EXEEXTENSION)

# add GMP-C++ bindings
ifeq ($(GMP),true)
LDFLAGS   +=  -lgmpxx
endif


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(SBCSFILE) $(SBCSOBJFILES)
endif

sbcs:		$(SBCSFILE)

.PHONY: all
all:            $(SCIPDIR) $(SBCSFILE)


# ----------- tests --------------
.PHONY: test testcluster
test:
	cd check; \
	$(SHELL) ./check.sh $(TEST) $(SBCSFILE) $(SETTINGS) $(notdir $(SBCSFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(SETCUTOFF) $(ONLYPRE) $(SEEDS);

testcluster:
	cd check; \
	$(SHELL) ./check_cluster.sh $(TEST) $(SBCSFILE) $(SETTINGS) \
	$(notdir $(SBCSFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) \
	$(DISPFREQ) $(CONTINUE) $(VERSION) \
	$(LPS) $(QUEUE) $(QUEUETYPE) $(PPN) $(CLIENTTMPDIR) $(NOWAITCLUSTER) $(PERMUTE) $(SETCUTOFF) $(SEEDS) $(SETUSESOL);


# ----------- others -------------
.PHONY: doc
doc:
		cd doc; $(DOXY) scip_sbcs.dxy

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.cpp src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: clean
clean:
		-rm -rf $(OBJDIR)/*

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(SBCSSRCFILES) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(SBCSDEP)'

# include dependencies
-include	$(SBCSDEP)



# link targets
$(SBCSFILE): $(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(SBCSOBJFILES)
		@echo "-> linking $@"
ifdef LINKCCSCIPALL
		-$(CXX) $(SBCSOBJFILES) $(LINKCCSCIPALL) -o $@
else
		-$(CXX) $(SBCSOBJFILES) -L$(SCIPDIR)/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) -l$(NLPILIB) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) -o $@
endif

#--------------------------- lint -------------------------------------------

.PHONY: lint
lint:	$(SBCSSRCFILES)
	-rm -f lint.out
	$(SHELL) -ec 'for i in $^; \
		do \
		echo $$i; \
		$(LINT) -i lint co-gcc.lnt +os\(lint.out\) -u -zero \
		-Iinclude -Ilib $(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE -DSYMGITHASH="" $$i; \
		done'

# ----------------------------------------------------------------------------
$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(CC_c) $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@

# --- EOF ---------------------------------------------------------------------
