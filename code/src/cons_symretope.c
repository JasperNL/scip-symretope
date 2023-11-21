/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_symretope.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for symmetry handling constraints for symretopes
 * @author Christopher Hojny
 * @author Jasper van Doornmalen
 *
 * The type of constraints of this constraint handler is described in cons_symretope.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_symresack.h"
#include "scip/cons_setppc.h"
#include "cons_symretope.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include <ctype.h>
#include <string.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "symretope"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on symretopes"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

#define DEFAULT_FORCECONSCOPY     FALSE /**< whether symretope constraints should be forced to be copied to sub SCIPs */
#define DEFAULT_SYMRETOPEPEEK      TRUE /**< Execute peek-procedure for symretope propagation */
#define DEFAULT_SYMRETOPEMAXORDER 10000 /**< Maximal group order for symretope. */
#define DEFAULT_SYMRETOPEMAXORDERNVARS 5000000 /**< Maximal group order * number of vars in symretope support. */
#define DEFAULT_SEPAALLVIOLPERMS   TRUE /**< Whether all violating permutations should be separated, or only the first */
#define DEFAULT_PROBINGPEEK       FALSE /**< Whether peeking should be done during probing. */

/* event handler properties */
#define EVENTHDLR_SYMRETOPE_NAME    "symretope"
#define EVENTHDLR_SYMRETOPE_DESC    "mark symretope constraint for propagation"

/* Constants to store fixings */
#define FIXED0    1                     /* When a variable is fixed to 0. */
#define FIXED1    2                     /* When a variable is fixed to 1. */
#define UNFIXED   0                     /* When a variable is neither fixed to 0 or to 1. */
#define FIXEDMAX  4                     /* An upper bound on the bitwise or's of FIXED0, FIXED1 and UNFIXED. */

/* A macro for checking if a variable was fixed before a bound-change index */
#define ISFIXED(x, bdchgidx)   (SCIPvarGetUbAtIndex(x, bdchgidx, FALSE) - SCIPvarGetLbAtIndex(x, bdchgidx, FALSE) < 0.5)

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   maxnvars;           /**< maximal number of variables in a symresack constraint */
   SCIP_Bool             forceconscopy;      /**< whether symresack constraints should be forced to be copied to sub SCIPs */
   SCIP_Bool             symretopepeek;      /**< whether symretope constraints should test unfixed variables whether feasibility applies if they are fixed to 0 or 1 */
   int                   maxorder;           /**< maximal group order for symretope */
   int                   maxordernvars;      /**< maximal group order for symretope multiplied with group support size */
   SCIP_Bool             sepaallviolperms;   /**< Whether a separating inequality should be added only for one violated symresack (FALSE) or for all violating symresacks (TRUE) */
   SCIP_Bool             probingpeek;        /**< Whether peeking should be done during probing. */
   SCIP_EVENTHDLR*       eventhdlr;          /**< An event handler for deciding whether a constraint must be propagated. */
};

enum SCIP_SymretopeGraphNodeType
{
   SYMRETOPE_ROOT = 0,
   SYMRETOPE_COND = 1,
   SYMRETOPE_NECC = 2,
};
typedef enum SCIP_SymretopeGraphNodeType SCIP_SYMRETOPEGRAPHNODETYPE;

typedef struct SCIP_SymretopeGraphNode_ SCIP_SYMRETOPEGRAPHNODE;
struct SCIP_SymretopeGraphNode_
{
   SCIP_SYMRETOPEGRAPHNODE*     predecessor;
   SCIP_SYMRETOPEGRAPHNODE*     successor1;
   SCIP_SYMRETOPEGRAPHNODE*     successor2;
   SCIP_SYMRETOPEGRAPHNODETYPE  nodetype;
   int                          fixing;
};
typedef struct SCIP_SymretopeGraphNode_ SCIP_SymretopeGraphNode;

/** constraint data for symresack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   SCIP_PERMUTATION*     permutation;        /**< the constraint permutation */
   SCIP_Bool             ismodelcons;        /**< whether the symresack is a model constraint */
#ifdef SCIP_DEBUG
   int                   debugcnt;           /**< counter to store number of added cover inequalities */
#endif
   int                   nperms;             /**< The number of permutations of the group to consider */
   SCIP_Bool             execprop;           /**< Whether we should propagate for this constraint. */
   SCIP_Bool*            affectedentries;    /**< For each variable, whether it is interesting to check. */
   SCIP_EVENTDATA*       vareventdata;       /**< Variable data for each event. */
};

/** Eventhandler data */
struct SCIP_EventData
{
   int varid;                                /**< Variable ID of the variable for which this event is added. */
   SCIP_CONSDATA* consdata;                  /**< Pointer to the associated constraint data */
};

struct SCIP_SymretopeVirtualFixings
{
   int* entrystack;                          /**< A stack of entries that are not unfixed in the virtual fixings. */
   int* entrylookup;                         /**< A lookup structure stating which entry is fixed. */
   int nvirtualfixings;                      /**< The number of virtual fixings */
   #ifndef NDEBUG
   int nvars;                                /**< The number of variables, i.e. the array sizes. Debug only. */
   #endif
};

typedef struct SCIP_SymretopeVirtualFixings SCIP_SYMRETOPEVIRTUALFIXINGS;

struct SCIP_SymretopeGraph
{
   SCIP_SymretopeGraphNode*     permgraphroots; /**< Stores all the tree-roots */
   SCIP_SymretopeGraphNode**    permgraphleaves; /**< Stores pointers to all the leaves */
   SCIP_SymretopeGraphNode*     permgraphs; /**< Stores the internal nodes */
   int*                         permpows; /**< Stores the powers of the permutations */
   SCIP_Bool*                   permsinqueue; /**< Takes permutation index, yields if it's already queued. */
   int*                         permsqueue; /**< The first permsqueuesize elements are the indices of the permutations in the queue. */
   int                          permsqueuesize; /**< The size of the permutation queue */
   int*                         permindices; /**< For each permutation the index state */
   #ifndef NDEBUG
   int                          nvars;    /**< The number of variables to account for the internal nodes. Debug only. */
   int                          maxnperms; /**< The maximal value of nperms that is supported. Debug only. */
   #endif
};
typedef struct SCIP_SymretopeGraph SCIP_SYMRETOPEGRAPH;

struct SCIP_FixingQueue
{
   int* fixinginqueue; /* Takes values UNFIXED (0), FIXED0 (1), FIXED1 (2) or FIXED0|FIXED1 (3) per index */
   int* fixingqueue; /* The first fixingqueuesize elements are the fixings that have to be applied. */
   int* fixingpermpows; /* For conflict analysis: What is the permutation (given by the power) 'causing' the fixing? */
   int fixingqueuesize; /* The number of fixings in the queue. */
};
typedef struct SCIP_FixingQueue SCIP_FIXINGQUEUE;

/*
 * Local methods
 */

/*
 * For virtual fixings
 */
static
SCIP_RETCODE allocVirtualFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMRETOPEVIRTUALFIXINGS** vf,        /**< pointer to virtual fixings data */
   int                   nvars               /**< number of variables */
)
{
   assert( scip != NULL );
   assert( vf != NULL );
   assert( nvars >= 0 );

   SCIP_CALL( SCIPallocBuffer(scip, vf) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*vf)->entrystack), nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &((*vf)->entrylookup), nvars) );
   (*vf)->nvirtualfixings = 0;
   #ifndef NDEBUG
   (*vf)->nvars = nvars;
   #endif
   return SCIP_OKAY;
}

static
void setVirtualFixing(
   SCIP_SYMRETOPEVIRTUALFIXINGS*  vf,        /**< pointer to virtual fixings data */
   int                            entry,     /**< the entry of the fixing */
   int                            value      /**< the fixing value, which may not be UNFIXED. */
)
{
   int* curval;

   assert( vf != NULL );
   assert( vf->nvirtualfixings >= 0 );
   assert( entry >= 0 );
   assert( entry < vf->nvars );
   assert( value != UNFIXED );

   /* Lookup the entry */
   curval = &(vf->entrylookup[entry]);

   if ( *curval == UNFIXED )
      /* This is not a virtual fixing yet. Mark it as such. */
      vf->entrystack[vf->nvirtualfixings++] = entry;
   *curval |= value;
}

static
int getVirtualFixing(
   SCIP_SYMRETOPEVIRTUALFIXINGS*  vf,        /**< pointer to virtual fixings data */
   int                            entry      /**< the entry of the fixing */
)
{
   assert( vf != NULL );
   assert( entry >= 0 );
   return vf->entrylookup[entry];
}

static
void clearVirtualFixings(
   SCIP_SYMRETOPEVIRTUALFIXINGS* vf          /**< pointer to virtual fixings data */
)
{
   assert( vf != NULL );
   assert( vf->nvirtualfixings >= 0 );
   while (vf->nvirtualfixings > 0)
   {
      /* Get the entry of the last element in the stack, and set the fixing hereof to UNFIXED. */
      vf->entrylookup[vf->entrystack[--(vf->nvirtualfixings)]] = UNFIXED;
   }
}

static
void copyVirtualFixings(
   SCIP_SYMRETOPEVIRTUALFIXINGS* from,
   SCIP_SYMRETOPEVIRTUALFIXINGS* to
)
{
   int i;
   int entry;
   assert( from != NULL );
   assert( to != NULL );
   clearVirtualFixings(to);
   for (i = 0; i < from->nvirtualfixings; ++i)
   {
      entry = from->entrystack[i];
      setVirtualFixing(to, entry, from->entrylookup[entry]);
   }
}

static
void freeVirtualFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMRETOPEVIRTUALFIXINGS** vf         /**< pointer to virtual fixings data */
)
{
   assert( vf != NULL );
   assert( (*vf)->nvirtualfixings >= 0 );
   clearVirtualFixings(*vf);

   SCIPfreeCleanBufferArray(scip, &((*vf)->entrylookup));
   SCIPfreeBufferArray(scip, &((*vf)->entrystack));
   SCIPfreeBuffer(scip, vf);
}


/*
 * For implication graphs
 */
static
SCIP_RETCODE allocSymretopeGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMRETOPEGRAPH** symgraph,           /**< pointer to symretope graph data */
   int                   nvars,              /**< number of variables */
   int                   maxnperms           /**< maximal number of nperms the graph structure should support. */
)
{
   assert( scip != NULL );
   assert( symgraph != NULL );
   assert( nvars >= 0 );
   assert( maxnperms >= 0 );

   SCIP_CALL( SCIPallocBuffer(scip, symgraph) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permpows), maxnperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permgraphroots), maxnperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permgraphleaves), 2 * maxnperms) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &((*symgraph)->permgraphs), 2 * nvars * maxnperms ) );

   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permsqueue), maxnperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permsinqueue), maxnperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*symgraph)->permindices), maxnperms) );

   #ifndef NDEBUG
   (*symgraph)->nvars =  nvars;
   (*symgraph)->maxnperms = maxnperms;
   (*symgraph)->permsqueuesize = -1;
   #endif

   return SCIP_OKAY;
}

static
void freeSymretopeGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMRETOPEGRAPH** symgraph            /**< pointer to symretope graph data */
)
{
   #ifndef NDEBUG
   int i;
   SCIP_SymretopeGraphNode* node;
   #endif

   assert( scip != NULL );
   assert( symgraph != NULL );

   #ifndef NDEBUG
   /* Sanity check: Make sure this memory is completely nulled by now. */
   for (i = 0; i < 2 * (*symgraph)->nvars * (*symgraph)->maxnperms; ++i)
   {
      node = &((*symgraph)->permgraphs[i]);
      assert( node->fixing == 0 );
      assert( node->nodetype == 0 );
      assert( node->predecessor == NULL );
      assert( node->successor1 == NULL );
      assert( node->successor2 == NULL );
   }
   #endif

   SCIPfreeBufferArray(scip, &((*symgraph)->permindices));
   SCIPfreeBufferArray(scip, &((*symgraph)->permsinqueue));
   SCIPfreeBufferArray(scip, &((*symgraph)->permsqueue));
   SCIPfreeCleanBufferArray(scip, &((*symgraph)->permgraphs));
   SCIPfreeBufferArray(scip, &((*symgraph)->permgraphleaves));
   SCIPfreeBufferArray(scip, &((*symgraph)->permgraphroots));
   SCIPfreeBufferArray(scip, &((*symgraph)->permpows));
   SCIPfreeBuffer(scip, symgraph);
}


/*
 * For fixing queues
 */
static
SCIP_RETCODE allocFixingQueue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FIXINGQUEUE**    fq,                 /**< pointer to the fixing queue data */
   int                   nvars               /**< number of variables */
)
{
   SCIP_CALL( SCIPallocBuffer(scip, fq) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &((*fq)->fixinginqueue), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*fq)->fixingqueue), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &((*fq)->fixingpermpows), nvars) );
   (*fq)->fixingqueuesize = 0;

   return SCIP_OKAY;
}

static
void freeFixingQueue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FIXINGQUEUE**    fq                  /**< pointer to the fixing queue data */
)
{
   SCIPfreeBufferArray(scip, &((*fq)->fixingpermpows));
   SCIPfreeBufferArray(scip, &((*fq)->fixingqueue));
   SCIPfreeCleanBufferArray(scip, &((*fq)->fixinginqueue));
   SCIPfreeBuffer(scip, fq);
}

/** frees a symretope constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to symretope constraint data */
   SCIP_CONSHDLR*        conshdlr            /**< symretope constraint handler */
   )
{
   SCIP_PERMUTATION* permutation;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( *consdata != NULL );
   assert( conshdlr != NULL );

   nvars = (*consdata)->nvars;

   if ( nvars == 0 )
   {
      assert( (*consdata)->vars == NULL );
      assert( (*consdata)->permutation == NULL );

      SCIPfreeBlockMemory(scip, consdata);

      return SCIP_OKAY;
   }

   /* free event handler */
   if ( SCIPisTransformed(scip) )
   {
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[i],
            SCIP_EVENTTYPE_VARCHANGED,
            conshdlrdata->eventhdlr, &(*consdata)->vareventdata[i], -1) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->affectedentries), nvars );
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->vareventdata), nvars );
   }

   /* free permutation stuff */
   permutation = (*consdata)->permutation;

   assert( permutation != NULL );
   assert( permutation->perm != NULL );
   SCIPfreePermutationContents(scip, permutation, TRUE);
   SCIPfreeBlockMemory(scip, &permutation);

   for (i = 0; i < nvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nvars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/*
 * Event handler callback method
 */

/** exec the event handler for handling variable bound changes
 *
 * Propagation is complete, so one does not need to propagate for symretopes if the variable bounds of the affected
 * constraints remain unchanged.
 *
 */
static
SCIP_DECL_EVENTEXEC(eventExec)
{
   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_SYMRETOPE_NAME) == 0 );
   assert( event != NULL );

   /* If the variable is impactful for the constraint (i.e. in the last propagator run it affected the outcome),
    * Then mark the propagator to run again. */
   if ( !eventdata->consdata->execprop && eventdata->consdata->affectedentries[eventdata->varid] )
   {
      eventdata->consdata->execprop = TRUE;
   }
   return SCIP_OKAY;
}

/** creates symretope constraint data
 *
 *  If the input data contains non-binary variables or fixed
 *  points, we delete these variables in a preprocessing step.
 */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< symretope constraint handler */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       inputvars,          /**< input variables of the constraint handler */
   int                   inputnvars,         /**< input number of variables of the constraint handler*/
   int*                  inputperm,          /**< input permutation of the constraint handler */
   SCIP_Bool             ismodelcons         /**< whether the symretope is a model constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** vars;
   int* indexcorrection;
   int* perm;
   int naffectedvariables;
   int i;
   int j = 0;
   SCIP_PERMUTATION* permutation;

   assert( consdata != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

#ifdef SCIP_DEBUG
   (*consdata)->debugcnt = 0;
#endif

   (*consdata)->ismodelcons = ismodelcons;

   /* COMMENT: You need to catch the case inputnvars == 0, cf. merge request 2660 in SCIP */
   /* count the number of binary variables which are affected by the permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &indexcorrection, inputnvars) );
   indexcorrection[0] = -1;
   for (i = 0; i < inputnvars; ++i)
   {
      if ( inputperm[i] != i && SCIPvarIsBinary(inputvars[i]) )
      {
         if ( i == 0 )
            indexcorrection[i] = 0;
         else
            indexcorrection[i] = indexcorrection[i - 1] + 1;
      }
      else
      {
         if ( i > 0 )
            indexcorrection[i] = indexcorrection[i - 1];
      }
   }
   naffectedvariables = indexcorrection[inputnvars - 1] + 1;

   (*consdata)->nvars = naffectedvariables;

   /* Stop if we detect that the permutation fixes each binary point. */
   if ( naffectedvariables == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indexcorrection);

      (*consdata)->vars = NULL;
      (*consdata)->nperms = 0;
      (*consdata)->nperms = 0;
      (*consdata)->permutation = NULL;
      return SCIP_OKAY;
   }

   /* remove fixed points from permutation representation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, naffectedvariables) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, naffectedvariables) );
   for (i = 0; i < inputnvars; ++i)
   {
      if ( i == 0 )
      {
         if ( indexcorrection[i] > -1 )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
      else
      {
         if ( indexcorrection[i] > indexcorrection[i - 1] )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
   }
   SCIPfreeBufferArrayNull(scip, &indexcorrection);

   SCIP_CALL( SCIPallocBlockMemory(scip, &permutation) );
   (*consdata)->permutation = permutation;

   SCIP_CALL( SCIPgetPermutation(scip, perm, naffectedvariables, permutation) );

   SCIPdebugMessage("Permutation: nvars=%d; ncycles=%d; order=%lld; ismonotone=%d; isordered=%d\n",
      permutation->nvars, permutation->ncycles, permutation->order, permutation->ismonotone, permutation->isordered);

   /* Specify the number of non-idenity permutations from the group to consider.
    * Note: the group order can be exponentially large;
    * If the group order is too large, then the number of considered permutations is limited.
    * TODO: Add warning, as conscheck can say "yes this is lexmax" while it is not.
    */
   assert( conshdlr != NULL );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   (*consdata)->nperms = MIN(permutation->order - 1, INT_MAX);
   if ( (*consdata)->nperms > conshdlrdata->maxorder
      || (*consdata)->nperms * (*consdata)->nvars > conshdlrdata->maxordernvars )
      SCIPwarningMessage(scip, "Symretope constraint will not capture all symmetries.\n");
   if ( conshdlrdata->maxorder > 0 && (*consdata)->nperms > conshdlrdata->maxorder )
   {
      (*consdata)->nperms = conshdlrdata->maxorder;
      SCIPwarningMessage(scip, "=> The symmetry group order %lld is larger than maxorder: %d. "
         "Restricting to %d permutations.\n",
         permutation->order, conshdlrdata->maxorder, (*consdata)->nperms);
   }
   if ( conshdlrdata->maxordernvars > 0 && (*consdata)->nperms * (*consdata)->nvars > conshdlrdata->maxordernvars )
   {
      (*consdata)->nperms = conshdlrdata->maxordernvars / (*consdata)->nvars ;
      /* In the extreme case that we have so many variables in the cycle that the integer division yields 0,
       * we should run at least for a single permutation.
       */
      if ( (*consdata)->nperms <= 0 )
         (*consdata)->nperms = 1;
      SCIPwarningMessage(scip, "=> The symmetry group order * cardinality of support (%lld * %d) "
         "is larger than maxordernvars: %d. Restricting to %d permutations.\n",
         permutation->order, (*consdata)->nvars, conshdlrdata->maxordernvars, (*consdata)->nperms);
   }

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      /* Make sure that all variables cannot be multiaggregated (cannot be handled by cons_symretope, since one cannot
       * easily eliminate single variables from a symretope constraint.
       */
      for (i = 0; i < naffectedvariables; ++i)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, vars[i], &vars[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[i]) );
      }

      /* Add events */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->vareventdata), naffectedvariables) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->affectedentries), naffectedvariables) );
      for (i = 0; i < naffectedvariables; ++i)
      {
         SCIP_EVENTDATA* vareventdata;
         vareventdata = &(*consdata)->vareventdata[i];
         vareventdata->varid = i;
         vareventdata->consdata = *consdata;
         SCIP_CALL( SCIPcatchVarEvent(scip, vars[i],
            SCIP_EVENTTYPE_VARCHANGED,
            conshdlrdata->eventhdlr, vareventdata, NULL) );
      }

      /* Mark that we want to propagate. */
      (*consdata)->execprop = TRUE;
   }
   else
   {
      (*consdata)->vareventdata = NULL;
      (*consdata)->affectedentries = NULL;
      (*consdata)->execprop = FALSE;
   }

   for (i = 0; i < naffectedvariables; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
   }
   (*consdata)->vars = vars;

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the ordering inequality for the pair \f$(1, \gamma^{-1}(1))\f$, i.e.,
 *  the inequality \f$-x_{1} + x_{\gamma^{-1}(1)} \leq 0\f$.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_PERMUTATION* permutation;
   int nvars;
   int cycleid;
   int cyclen;
   int* cycle;
   int k;
#ifdef SCIP_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   /* avoid stupid problems */
   if ( nvars <= 1 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   vars = consdata->vars;

   assert( consdata->permutation != NULL );
   permutation = consdata->permutation;

   /* Get cycle of the first variable */
   cycleid = permutation->varcycle[0];
   cyclen = permutation->cyclelengths[cycleid];
   cycle = permutation->cycles[cycleid];

   /* For every entry k in the cycle that is not entry 0, add x[0] <= x[k]; */
   for (k = 0; k < cyclen; ++k)
   {
      SCIP_ROW* row;

      if ( cycle[k] == 0 )
         continue;

      /* add ordering inequality */
      #ifdef SCIP_DEBUG
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_init_%s_%d", SCIPconsGetName(cons), k);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
      #else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
      #endif

      SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], -1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars[cycle[k]], 1.0) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/* Helper function:
 * During conflict analysis, list the variables that cause that "boundtype" cannot be tightened for "infervar".
 */
static
SCIP_RETCODE resolveSymretopeConflictVariables(
   SCIP* scip,                               /**< SCIP pointer */
   SCIP_VAR* infervar,                       /**< The variable for which conflict analysis is executed */
   SCIP_BOUNDTYPE boundtype,                 /**< Whether the lower or upper bound of infervar is set */
   SCIP_VAR** vars,                          /**< The array of variables */
   int nvars,                                /**< The number of variables */
   SCIP_PERMUTATION* permutation,            /**< Permutation information */
   int permpow,                              /**< Power of the permutation for which the conflict should follow */
   SCIP_BDCHGIDX* bdchgidx                   /**< Bound index at which the fixing of infervar is made. */
)
{
   int i;
   int j;
   #ifndef NDEBUG
   int infervarid;
   #endif
   int* virtualfixings;

   /* Store all fixings part of our conflict. */
   SCIP_CALL( SCIPallocBufferArray(scip, &virtualfixings, nvars) );
   #ifndef NDEBUG
   infervarid = -1;
   #endif
   for (i = 0; i < nvars; ++i)
   {
      virtualfixings[i] = UNFIXED;
      if ( vars[i] == infervar )
      {
         #ifndef NDEBUG
         infervarid = i;
         #endif
         /* If variable i gets fixed to 0, then upper bound is set.
          * If variable i gets fixed to 1, then the lower bound is set.
          * This is the case, because one cannot set the lower bound in an effective manner if we want to fix to 0,
          * because this is the global lower bound, already.
          */
         if ( boundtype == SCIP_BOUNDTYPE_UPPER )
            virtualfixings[i] = FIXED1;
         /* Same: If variable i gets fixed to 1, then the lower bound is set.
          * Infeasibility is found if we assume that it is 0.
          */
         else
         {
            assert( boundtype == SCIP_BOUNDTYPE_LOWER );
            virtualfixings[i] = FIXED0;
         }
      }
   }
   assert( infervar == NULL || infervarid >= 0 );

   for (i=0; i < nvars; ++i)
   {
      j = permGet(permutation, i, -permpow);
      assert( j >= 0 );
      assert( j < nvars );

      /* Ignore fixed points of permutation. */
      if ( i == j )
         continue;

      /* If infeasibility is found. */
      if ( virtualfixings[i] == FIXED0 && virtualfixings[j] == FIXED1 )
         break;

      /* If var i is fixed to 0, then var j should be fixed to 0 by propagation. */
      if ( virtualfixings[i] == FIXED0 )
      {
         /* Is j fixed to 1? Then it's infeasible. */
         if ( virtualfixings[j] == FIXED1 )
            break;

         if ( SCIPvarGetLbAtIndex(vars[j], bdchgidx, FALSE) > 0.5 )
         {
            SCIP_CALL( SCIPaddConflictLb(scip, vars[j], bdchgidx) );
            break;
         }

         /* Var j is not fixed to 1, so it gets propagated to 0. */
         virtualfixings[j] = FIXED0;
         continue;
      }

      /* If var j is fixed to 1, then var j is fixed to 1 by propagation. */
      if ( virtualfixings[j] == FIXED1 )
      {
         /* Is i fixed to 0? Then it's infeasible. */
         if ( virtualfixings[i] == FIXED0 )
            break;

         if ( SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            break;
         }

         /* Var i is not fixed to 0, so it gets propagated to 1. */
         virtualfixings[i] = FIXED1;
         continue;
      }

      /* (1, 0) is not possible. */
      assert( !(virtualfixings[i] == FIXED1 && virtualfixings[j] == FIXED0) );

      /* Remaining cases are (1, _), (_, 0) or (_, _). In these cases, we are not aware of the bounds (yet). */

      /* If var i is fixed to 0. Then we have (0, 01_), must be (0, 0). */
      if ( SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) < 0.5 )
      {
         assert( SCIPvarGetLbAtIndex(vars[i], bdchgidx, 0) < 0.5 );
         SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
         virtualfixings[i] = FIXED0;

         /* If it turns out to be (0, 1), break. */
         if ( SCIPvarGetLbAtIndex(vars[j], bdchgidx, FALSE) > 0.5 )
         {
            /* COMMENT: You don't need to set virtualfixings, because you break afterwards anyways */
            virtualfixings[j] = FIXED1;
            if ( vars[j] != infervar )
            {
               SCIP_CALL( SCIPaddConflictLb(scip, vars[j], bdchgidx) );
            }
            break;
         }
         /* Otherwise, it's implied. */
         else
            virtualfixings[j] = FIXED0;
      }

      /* If var j is fixed to 1. Then we have (01_, 1), must be (1, 1). */
      if ( SCIPvarGetLbAtIndex(vars[j], bdchgidx, 0) > 0.5 )
      {
         assert( SCIPvarGetUbAtIndex(vars[j], bdchgidx, 0) > 0.5 );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[j], bdchgidx) );
         virtualfixings[j] = FIXED1;

         /* If it turns out to be (0, 1), break. */
         if ( SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) < 0.5 )
         {
            /* COMMENT: You don't need to set virtualfixings, because you break afterwards anyways */
            virtualfixings[i] = FIXED0;
            if ( vars[i] != infervar )
            {
               SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            }
            break;
         }
         /* Otherwise, it's implied. */
         else
            virtualfixings[i] = FIXED1;
      }

      /* Now it must be a constant row, and not free. */
      assert( virtualfixings[i] != UNFIXED );
      assert( virtualfixings[j] != UNFIXED );
      assert( !(virtualfixings[i] == FIXED1 && virtualfixings[j] == FIXED0) );
      assert( virtualfixings[i] == virtualfixings[j] );
   }
   /* We must have seen "infervarid", unless we have no infervar. */
   assert( infervar == NULL || i >= infervarid || i >= permGet(permutation, infervarid, permpow) );
   /* The loop must always break somewhere. */
   assert( i < nvars );

   SCIPfreeBufferArray(scip, &virtualfixings);
   return SCIP_OKAY;
}

/** Helper function: Get the variable fixing */
static
int getVarFixing(
   SCIP_VAR** vars,                          /**< The array of variables */
   int varid,                                /**< The variable ID for which the fixing is sought after */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< The virtual fixings structure, or NULL if fixings are applied globally. */
   SCIP_Bool useproblembounds,               /**< Whether local bounds of the problem instance should be used. */
   SCIP_Bool* checkedentries                 /**< For each variable entry, whether the entry has been looked up */
)
{
   SCIP_VAR* var;
   assert( vars != NULL );
   /* if useproblembounds is FALSE, then virtualfixings must be defined. */
   assert( useproblembounds ? TRUE : virtualfixings != NULL );

   /* Mark entry as checked */
   if ( checkedentries != NULL )
      checkedentries[varid] = TRUE;

   /* First, check if the entry is fixed virtually */
   if ( virtualfixings != NULL )
   {
      #ifndef NDEBUG
      /* For debugging purposes, we do need var here. */
      var = vars[varid];
      assert( var != NULL );
      #endif

      switch (getVirtualFixing(virtualfixings, varid))
      {
         case FIXED0:
            assert( !useproblembounds || SCIPvarGetLbLocal(var) < 0.5 ); /* A zero-fixing must be possible */
            return FIXED0;
         case FIXED1:
            assert( !useproblembounds || SCIPvarGetUbLocal(var) > 0.5 ); /* A one-fixing must be possible */
            return FIXED1;
         case UNFIXED:
            break;
         default:
            assert( FALSE );  /* It must be FIXED0, FIXED1 or UNFIXED. */
            break;
      }
   }

   /* The entry is not fixed virtually. Check the bounds of the problem (the local bounds). */
   if ( useproblembounds )
   {
      var = vars[varid];
      assert( var != NULL );
      if ( SCIPvarGetLbLocal(var) > 0.5 )
      {
         assert( SCIPvarGetUbLocal(var) > 0.5 );
         assert( virtualfixings == NULL || (getVirtualFixing(virtualfixings, varid) & FIXED0) == 0 );
         if ( virtualfixings != NULL )
            setVirtualFixing(virtualfixings, varid, FIXED1);
         return FIXED1;
      }
      else if ( SCIPvarGetUbLocal(var) < 0.5 )
      {
         assert( SCIPvarGetLbLocal(var) < 0.5 );
         assert( virtualfixings == NULL || (getVirtualFixing(virtualfixings, varid) & FIXED1) == 0 );
         if ( virtualfixings != NULL )
            setVirtualFixing(virtualfixings, varid, FIXED0);
         return FIXED0;
      }
   }

   /* The only case remains: The variable is not fixed. */
   return UNFIXED;
}


/** Helper function: Set the variable fixing.
 * If @p virtualfixings is not NULL, then the fixings are applied on the problem.
 * Otherwise, the fixings are applied on the virtual fixings defined by @p virtualfixings .
 */
static
SCIP_RETCODE setVarFixing(
   SCIP* scip,                            /**< SCIP pointer */
   SCIP_CONS* cons,                       /**< Symretope constraint pointer */
   SCIP_VAR** vars,                       /**< The variable array pointer */
   int varid,                             /**< The index of the variable in the variable array that is fixed */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< The virtual fixings structure, or NULL */
   int fixing,                            /**< The fixing, which is either FIXED0 or FIXED1 */
   SCIP_Bool* infeasible,                 /**< To store whether infeasibility is found */
   SCIP_Bool* tightened,                  /**< To store whether the problem is tightened */
   int inferinfo                          /**< Information for conflict resolution. */
)
{
   assert( vars != NULL );
   assert( fixing == FIXED0 || fixing == FIXED1 );
   assert( infeasible != NULL );
   assert( tightened != NULL );

   if ( virtualfixings == NULL )
   {
      /* Virtual fixings are not defined, so we apply the fixings. */
      if ( fixing == FIXED0 )
      {
         assert( vars[varid] != NULL );
         SCIP_CALL( SCIPinferVarUbCons(scip, vars[varid], 0.0, cons, inferinfo, FALSE, infeasible, tightened) ); /*lint !e713*/
      }
      else if ( fixing == FIXED1 )
      {
         assert( vars[varid] != NULL );
         SCIP_CALL( SCIPinferVarLbCons(scip, vars[varid], 1.0, cons, inferinfo, FALSE, infeasible, tightened) ); /*lint !e713*/
      }
   }
   else
   {
      *tightened = (getVirtualFixing(virtualfixings, varid) & fixing) == 0;
      setVirtualFixing(virtualfixings, varid, fixing);
      *infeasible = getVirtualFixing(virtualfixings, varid) == (FIXED0 | FIXED1);
   }
   return SCIP_OKAY;
}


#ifndef NDEBUG
/** Test function: If setVarFixing is called without virtual fixings, test for violation of claimed permutation. */
static
SCIP_RETCODE setVarFixingTest(
   SCIP* scip,                               /**< SCIP pointer */
   SCIP_VAR** vars,                          /**< Pointer to variables array */
   int nvars,                                /**< Number of variables in the array */
   SCIP_PERMUTATION* permutation,            /**< Pointer to the permutation information */
   int varid,                                /**< The variable index that is fixed. */
   int fixing,                               /**< The value at which the variable gets fixed, either FIXED0 or FIXED1 */
   int permpow                               /**< The permutation power for which the variable fixing is found. */
)
{
   int i;
   int j;
   int fixi;
   int fixj;
   int* virtualfixings;

   assert( vars != NULL );
   assert( fixing == FIXED0 || fixing == FIXED1 );

   /* If inferinfo (permpow) is set < 0, then resprop is not executed. */
   if ( permpow < 0 )
      return SCIP_OKAY;

   /* Initialize virtual fixings as "everything unfixed" */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &virtualfixings, nvars) );

   /* Assume that we fix varid to the inverse fixing. */
   if ( fixing == FIXED0 )
      virtualfixings[varid] = FIXED1;
   else /* i.e. fixing == FIXED1, see assert above: fixing is either FIXED0 or FIXED1. */
      virtualfixings[varid] = FIXED0;

   for (i = 0; i < nvars; ++i)
   {
      j = permGet(permutation, i, -permpow);
      assert( vars[i] != NULL );
      assert( vars[j] != NULL );

      /* Ignore fixed points */
      if ( i == j )
         continue;

      /* Get fixing of i */
      fixi = virtualfixings[i];
      if ( fixi == UNFIXED )
      {
         assert( SCIPisLE(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) );
         if ( SCIPvarGetUbLocal(vars[i]) < 0.5 )
            fixi = FIXED0;
         else if ( SCIPvarGetLbLocal(vars[i]) > 0.5 )
            fixi = FIXED1;
      }

      /* Get fixing of j */
      fixj = virtualfixings[j];
      if ( fixj == UNFIXED )
      {
         assert( SCIPisLE(scip, SCIPvarGetLbLocal(vars[j]), SCIPvarGetUbLocal(vars[j])) );
         if ( SCIPvarGetUbLocal(vars[j]) < 0.5 )
            fixj = FIXED0;
         else if ( SCIPvarGetLbLocal(vars[j]) > 0.5 )
            fixj = FIXED1;
      }

      /* If (0, 1) is found: We're happy, this is indeed infeasible. */
      if ( fixi == FIXED0 && fixj == FIXED1 )
         break;

      /* If it could be made (1, 0), something's wrong. */
      assert( !( fixi == FIXED1 && fixj == FIXED0 ) );
      assert( !( fixi == UNFIXED && fixj == FIXED0 ) );
      assert( !( fixi == FIXED1 && fixj == UNFIXED ) );
      assert( !( fixi == UNFIXED && fixj == UNFIXED ) );

      if ( fixi == FIXED0 && fixj == UNFIXED )
         virtualfixings[j] = FIXED0;
      if ( fixj == FIXED1 && fixi == UNFIXED )
         virtualfixings[i] = FIXED1;
   }

   /* The loop must be terminated by the break-statement. */
   assert( i < nvars );

   SCIPfreeBufferArray(scip, &virtualfixings);
   return SCIP_OKAY;
}
#endif


/** Helper function: Add a fixing to the fixing queue. */
static
SCIP_RETCODE enqueueFixing(
   SCIP* scip,                               /**< SCIP pointer */
   SCIP_VAR** vars,                          /**< Pointer to variable array */
   SCIP_CONS* cons,                          /**< Pointer to symretope constraint */
   int fixing,                               /**< The encoded fixing that gets applied: varid + nvars * (0 or 1) */
   SCIP_PERMUTATION* permutation,            /**< The permutation information */
   int permpow,                              /**< The permutation power for which this fixing is found. */
   int nvars,                                /**< The number of variables */
   SCIP_FIXINGQUEUE* fq,                     /**< Pointer to the fixing queue */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< Pointer to the virtual fixings structure, or NULL if fixings are applied globally */
   SCIP_Bool* infeasible                     /**< Pointer to store whether infeasibility is found. */
   )
{
   int b;
   int i;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( cons != NULL );
   assert( fixing >= 0 && fixing < 2 * nvars );
   assert( nvars > 0 );
   assert( fq != NULL );
   assert( fq->fixingqueue != NULL );
   assert( fq->fixinginqueue != NULL );
   assert( fq->fixingqueuesize >= 0 );
   assert( infeasible != NULL );

   i = fixing % nvars;
   b = fixing >= nvars ? FIXED1 : FIXED0;
   assert( i >= 0 && i < nvars );
   if ( (fq->fixinginqueue[i] & b) == 0 )
   {
      /* Fixing is not queued yet.
       * Check if converse fixing is contained in queue (->infeasible), and otherwise add fixing to queue.
       */
      if ( (fq->fixinginqueue[i] | b) == (FIXED0 | FIXED1) )
      {
         /* Variable at index i must be fixed to 0 and 1 at the same time! A contradiction. */
         *infeasible = TRUE;

         /* Apply conflict analysis. Can only execute this when not peeking. */
         if ( virtualfixings == NULL && SCIPisConflictAnalysisApplicable(scip) )
         {
            int otherpermpow;
            assert( fq->fixingpermpows != NULL );
            otherpermpow = fq->fixingpermpows[i];

            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            if ( b == FIXED0 )
            {
               /* Add why permid wants to fix vars[i] to 0 (tightening UB),
                * and why otherpermid wants to fix vars[i] to 1.
                */
               SCIP_CALL( resolveSymretopeConflictVariables(scip, vars[i], SCIP_BOUNDTYPE_UPPER, vars, nvars,
                  permutation, permpow, NULL) );
               SCIP_CALL( resolveSymretopeConflictVariables(scip, vars[i], SCIP_BOUNDTYPE_LOWER, vars, nvars,
                  permutation, otherpermpow, NULL) );
            }
            else /* I.e. b == FIXED1 */
            {
               assert( b == FIXED1 );
               /* Add why permid wants to fix vars[i] to 1 (tightening LB),
                * and why otherpermid wants to fix vars[i] to 0.
                */
               SCIP_CALL( resolveSymretopeConflictVariables(scip, vars[i], SCIP_BOUNDTYPE_LOWER, vars, nvars,
                  permutation, permpow, NULL) );
               SCIP_CALL( resolveSymretopeConflictVariables(scip, vars[i], SCIP_BOUNDTYPE_UPPER, vars, nvars,
                  permutation, otherpermpow, NULL) );
            }

            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }
         return SCIP_OKAY;
      }
      /* Mark that the fixing of variable "i" is due to the permutation with power permpow. */
      fq->fixingpermpows[i] = permpow;
      /* Add this fixing to the queue */
      assert( fq->fixinginqueue[i] == UNFIXED );
      fq->fixinginqueue[i] |= b;
      fq->fixingqueue[(fq->fixingqueuesize)++] = fixing;
   }
   return SCIP_OKAY;
}

/** Given a node ( @p root ) in an implication tree, this function removes the subtree rooted in this node,
 * including the node itself.
 * With removing we mean that the nodes of the implication trees are set to NULL, along with the pointers to them.
 */
static
SCIP_RETCODE removeSubtree(
   SCIP* scip,                               /**< SCIP instance */
   SCIP_SymretopeGraphNode* root,            /**< The implication tree node whose subtree needs to be removed. */
   SCIP_SymretopeGraphNode** permgraphleaves /**< The array (of size two) containing both possible leaves of the
                                              * implication tree. By removing the subtree we also remove pointers
                                              * of these leaves. */
)
{
   /* Starting at "root", the root of the subtree, we remove it and go to each successor.
    * Do this until we cannot go any further. This is the "leaf" node (which should be handled separately).
    * Note: this function must also unlink root from its parent.
    */
   SCIP_SymretopeGraphNode* next;

   /* First, disconnect the "root" (which is the root of the subtree to be removed) from the tree. */
   /* By definition, the actual tree root has no predecessor, so only do this for internal nodes. */
   if ( root->nodetype != SYMRETOPE_ROOT )
   {
      if ( root == root->predecessor->successor1 )
      {
         root->predecessor->successor1 = root->predecessor->successor2;
         root->predecessor->successor2 = NULL;
      }
      else
      {
         assert( root == root->predecessor->successor2 );
         root->predecessor->successor2 = NULL;
      }
   }

   next = root;

   while(next != NULL)
   {
      root = next;
      if ( root->successor1 != NULL )
      {
         /* Assumption: This is a path. (That's not the case: Recurse?) */
         if ( root->successor2 != NULL )
         {
            SCIP_CALL( removeSubtree(scip, root->successor2, permgraphleaves) );
         }
         next = root->successor1;
      }
      else
      {
         next = root->successor2;
      }

      /* Now remove all attributes of "root". */
      root->fixing = 0;
      root->nodetype = 0;
      root->predecessor = NULL;
      root->successor1 = NULL;
      root->successor2 = NULL;
   }

   /* Now root is a leaf of the tree. Check if this is a "loose end", and if so remove. */
   if ( root == permgraphleaves[0] )
      permgraphleaves[0] = NULL;
   if ( root == permgraphleaves[1] )
      permgraphleaves[1] = NULL;

   return SCIP_OKAY;
}


static
SCIP_RETCODE applyFixings(
   SCIP* scip,                               /**< SCIP instance */
   SCIP_CONS* cons,                          /**< Symretope constraint instance */
   SCIP_VAR** vars,                          /**< Array of variables */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< Virtual fixings structure, or NULL if the fixings are applied to the problem */
   SCIP_PERMUTATION* permutation,            /**< Pointer to permutation information */
   SCIP_SymretopeGraphNode* permgraphs,      /**< Pointer to the data structure storing all implication tree */
   SCIP_SymretopeGraphNode** permgraphleaves,/**< Pointer to array storing all (at most 2) leaves per implication tree */
   int* permpows,                            /**< Array that stores which power belongs to which permutation index */
   int nvars,                                /**< The number of variables in the permutation */
   int nperms,                               /**< The total number of permutations to explore (does not include the identity) */
   SCIP_FIXINGQUEUE* fq,                     /**< Pointer to the fixing queue */
   SCIP_Bool* permsinqueue,                  /**< Pointer to the array for looking up which perms are in the queue */
   int* permsqueue,                          /**< Pointer to the permutation queue, which is implemented as a stack */
   int* permsqueuesize,                      /**< Pointer to integer storing the current number of permutations in the queue */
   int* ngen,                                /**< Pointer to store the number of propagations */
   SCIP_Bool* infeasible,                    /**< Pointer to store if infeasibility is found */
   SCIP_Bool* tightened                      /**< Pointer to store if the problem is tightened by the applied fixings */
)
{
   int fixingvalue;
   int fixingvaluealt;
   int fixingvarid;
   int j;
   int k;
   SCIP_SymretopeGraphNode* succ;
   SCIP_SymretopeGraphNode* node;
   SCIP_SymretopeGraphNode* twin;
   SCIP_SymretopeGraphNode* twinsucc;
   SCIP_SymretopeGraphNode* pred;
   SCIP_SymretopeGraphNode* permgraph;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars != NULL );
   assert( permutation != NULL );
   assert( permgraphs != NULL );
   assert( permgraphleaves != NULL );
   assert( permpows != NULL );
   assert( nvars >= 0 );
   assert( nperms >= 0 );
   assert( fq != NULL );
   assert( permsinqueue != NULL );
   assert( permsqueue != NULL );
   assert( permsqueuesize != NULL );
   assert( ngen != NULL );
   assert( infeasible != NULL );
   assert( tightened != NULL );

   while (fq->fixingqueuesize > 0)
   {
      /* Get the fixing encoding, which is i + nvars * b, where b \in {0, 1} */
      fixingvalue = fq->fixingqueue[--(fq->fixingqueuesize)];
      assert( fixingvalue >= 0 && fixingvalue < 2 * nvars );
      /* Extract fixing variable index (i) and the associated fixing value (which is b in the description above). */
      fixingvarid = fixingvalue % nvars;
      fixingvalue = fixingvalue >= nvars ? FIXED1 : FIXED0;
      fq->fixinginqueue[fixingvarid] &= ~fixingvalue;

      /* Apply this fixing. */
      SCIP_CALL( setVarFixing(scip, cons, vars, fixingvarid, virtualfixings, fixingvalue, infeasible, tightened,
         fq->fixingpermpows == NULL ? -1 : fq->fixingpermpows[fixingvarid]) );
      #ifndef NDEBUG
      if ( virtualfixings == NULL )
      {
         SCIP_CALL( setVarFixingTest(scip, vars, nvars, permutation, fixingvarid, fixingvalue,
            fq->fixingpermpows == NULL ? -1 : fq->fixingpermpows[fixingvarid]) );
      }
      #endif
      if ( *tightened )
         ++(*ngen);
      if ( *infeasible )
         return SCIP_OKAY;

      /* Update data structures */
      for (k=0; k < nperms; ++k)
      {
         permgraph = &permgraphs[2 * nvars * k];
         for (j = 0; j < 2; ++j)
         {
            /* COMMENT: To make reading the code easier and also to be more robust towards future code changes:
             * What about introducing functions SCIPgetPermGraph(k) that returns &permgraphs[2*nvars*k]
             * and SCIPgetPermGraphNode(fixvar, j) that returns &permgraph[2*fixvar+j]?
             * Then, it is clear what is supposed to happen.
             */
            node = &permgraph[2 * fixingvarid + j];

            /* If node does not exist. */
            if ( node->predecessor == NULL )
               continue;

            assert( node->fixing % nvars == fixingvarid );
            fixingvaluealt = node->fixing >= nvars ? FIXED1 : FIXED0;

            if ( fixingvalue == fixingvaluealt )
            {
               /* The fixing of the node is the same as the fixing that needs to be applied. */

               /* If "node" has a sibling (so it's conditional),
                  * then the subtree induced by that sibling is an infeasible subtree.
                  * Namely, if this conditional fixing is (i, 0), then the sibling is conditional with (j, 1),
                  * and this has a necessary fixing node as child, with fixing (i, 1).
                  */
               if ( node->predecessor->successor1 != NULL && node->predecessor->successor2 != NULL )
               {
                  assert( node->predecessor->successor1 != node->predecessor->successor2 );
                  assert( node->nodetype == SYMRETOPE_COND );

                  /* Select sibling of "node" */
                  twin = node->predecessor->successor1;
                  if ( twin == node )
                     twin = node->predecessor->successor2;

                  /* Remove whole subtree of the sibling. */
                  SCIP_CALL( removeSubtree(scip, twin, &permgraphleaves[2 * k]) );
               }

               /* After the previous part, "node" has no sibling. */
               /* COMMENT: Why is succ called succ? I would have expected it to be called pred, because it's a predecessor of node. */
               succ = node->predecessor;
               assert( node == succ->successor1 || node == succ->successor2 );
               assert( NULL == succ->successor1 || NULL == succ->successor2 );
               succ->successor1 = node->successor1;
               succ->successor2 = node->successor2;
               if ( node->successor1 != NULL )
                  node->successor1->predecessor = succ;
               if ( node->successor2 != NULL )
                  node->successor2->predecessor = succ;

               /* And remove "node" */
               node->fixing = 0;
               node->nodetype = 0;
               node->predecessor = NULL;
               node->successor1 = NULL;
               node->successor2 = NULL;

               /* If "node" is a leaf, then now "succ" is a leaf. */
               if ( node == permgraphleaves[2 * k + 0] )
                  permgraphleaves[2 * k + 0] = succ;
               if ( node == permgraphleaves[2 * k + 1] )
                  permgraphleaves[2 * k + 1] = succ;

               /* If "node" is directly after the root node, and if the successors are neccesary fixing nodes
                * then we must apply these fixings.
                */
               if ( succ->nodetype == SYMRETOPE_ROOT )
               {
                  if ( succ->successor1 != NULL && succ->successor1->nodetype == SYMRETOPE_NECC )
                  {
                     /* A fixing must be applied! */
                     SCIP_CALL( enqueueFixing(scip, vars, cons, succ->successor1->fixing, permutation, permpows[k],
                        nvars, fq, virtualfixings, infeasible) );
                     if ( *infeasible )
                        return SCIP_OKAY;
                  }
                  if ( succ->successor2 != NULL && succ->successor1->nodetype == SYMRETOPE_NECC )
                  {
                     /* A fixing must be applied! */
                     SCIP_CALL( enqueueFixing(scip, vars, cons, succ->successor2->fixing, permutation, permpows[k],
                        nvars, fq, virtualfixings, infeasible) );
                     if ( *infeasible )
                        return SCIP_OKAY;
                  }
               }
            }
            else
            {
               /* The fixing of "node" is the converse! */
               if ( node->nodetype == SYMRETOPE_NECC )
               {
                  /* This is now an infeasible subtree! */

                  /* Go to the first non-necessary fixing ancestor. */
                  /* COMMENT: Again succ vs. pred */
                  succ = node->predecessor;
                  assert( (succ->successor1 == node) ^ (succ->successor2 == node) );
                  SCIP_CALL( removeSubtree(scip, node, &permgraphleaves[2 * k]) );
                  while (succ->nodetype == SYMRETOPE_NECC)
                     succ = succ->predecessor;

                  if ( succ->nodetype == SYMRETOPE_ROOT )
                  {
                     /* There is a path from the root to "node", only visiting necessary fixing nodes.
                      * This is infeasible.
                      */
                     *infeasible = TRUE;

                     /* Conflict analysis: Infeasibility follows because we have fixed i to fixing,
                      * but permutation k says that we have to fix i to fixingalt, the converse.
                      * Only do such conflict analysis when not peeking.
                      */
                     if ( virtualfixings == NULL && SCIPisConflictAnalysisApplicable(scip) )
                     {
                        SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                        SCIP_CALL( resolveSymretopeConflictVariables(scip, NULL, SCIP_BOUNDTYPE_LOWER,
                           vars, nvars, permutation, permpows[k], NULL) );
                        SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
                     }

                     return SCIP_OKAY;
                  }
                  else
                  {
                     assert( succ->nodetype == SYMRETOPE_COND );
                     /* On the path from root to "node" there is a conditional fixing node. */

                     /* If this has no siblings, then replace by necessary fixing node with converse fixing */
                     pred = succ->predecessor;
                     assert( succ == pred->successor1 || succ == pred->successor2 );

                     /* Get his potential sibling. */
                     twin = pred->successor1;
                     if ( twin == succ )
                        twin = pred->successor2;

                     if ( twin == NULL )
                     {
                        /* Succ has no sibling. Replace by necessary fixing node with converse fixing. */
                        succ->nodetype = SYMRETOPE_NECC;
                        if ( succ->fixing >= nvars )
                           succ->fixing = succ->fixing - nvars;
                        else
                           succ->fixing = succ->fixing + nvars;
                        assert(succ->fixing >= 0 && succ->fixing < 2 * nvars);

                        /* Remove nodes after succ. */
                        if ( succ->successor1 != NULL )
                        {
                           SCIP_CALL( removeSubtree(scip, succ->successor1, &permgraphleaves[2 * k]) );
                        }
                        if ( succ->successor2 != NULL)
                        {
                           SCIP_CALL( removeSubtree(scip, succ->successor2, &permgraphleaves[2 * k]) );
                        }

                        /* If succ is now child of the root, then the fixing of succ must be applied. */
                        if ( pred->nodetype == SYMRETOPE_ROOT )
                        {
                           SCIP_CALL( enqueueFixing(scip, vars, cons, succ->fixing, permutation, permpows[k], nvars,
                              fq, virtualfixings, infeasible) );
                        }
                     }
                     else
                     {
                        /* Succ has a sibling, and it is twin. */
                        assert( twin != succ );
                        /* Twin must have exactly one successor. */
                        assert( (twin->successor1 != NULL) ^ (twin->successor2 != NULL) );

                        /* Get the successor of the twin. */
                        twinsucc = twin->successor1;
                        if ( twinsucc == NULL )
                           twinsucc = twin->successor2;
                        assert( twinsucc != NULL );

                        /* Now, by the setup, the fixing of twinsucc must be the converse of the fixing of succ. */
                        assert( abs(twinsucc->fixing - succ->fixing) == nvars );

                        /* Remove succ, and move twinsucc one place towards the root (so before twin). */
                        /* 1. Remove subtree rooted in succ */
                        SCIP_CALL( removeSubtree(scip, succ, &permgraphleaves[2 * k]) );

                        /* 2. Place twinsucc one place forward.
                         * The path becomes pred -> twinsucc -> twin -> successors of twinsucc.
                         */
                        assert( (twinsucc->successor1 == NULL) || (twinsucc->successor2 == NULL) );
                        assert( (twin->successor1 == NULL) ^ (twin->successor2 == NULL) );
                        assert( (twin->successor1 == twinsucc) ^ (twin->successor2 == twinsucc) );
                        assert( (twin->predecessor->successor1) == twin );
                        assert( (twin->predecessor->successor2) == NULL );
                        /* Fix predecessors */
                        twinsucc->predecessor = pred;
                        twin->predecessor = twinsucc;
                        if ( twinsucc->successor1 != NULL)
                           twinsucc->successor1->predecessor = twin;
                        if ( twinsucc->successor2 != NULL )
                           twinsucc->successor2->predecessor = twin;
                        /* Fix successors */
                        twin->successor1 = twinsucc->successor1;
                        twin->successor2 = twinsucc->successor2;
                        twinsucc->successor1 = twin;
                        twinsucc->successor2 = NULL;
                        pred->successor1 = twinsucc;
                        pred->successor2 = NULL;

                        /* If twinsucc was a leaf, then now twin is a leaf. */
                        if ( twinsucc == permgraphleaves[2 * k + 0] )
                           permgraphleaves[2 * k + 0] = twin;
                        if ( twinsucc == permgraphleaves[2 * k + 1] )
                           permgraphleaves[2 * k + 1] = twin;

                        /* If twinsucc is now child of the root, then the fixing of twinsucc must be applied. */
                        if ( pred->nodetype == SYMRETOPE_ROOT )
                        {
                           SCIP_CALL( enqueueFixing(scip, vars, cons, twinsucc->fixing, permutation, permpows[k],
                              nvars, fq, virtualfixings, infeasible) );
                        }
                     }
                  }
               }
               else
               {
                  /* Remove this whole subtree. */
                  assert( node->nodetype == SYMRETOPE_COND );
                  SCIP_CALL( removeSubtree(scip, node, &permgraphleaves[2 * k]) );
               }
            }
         }

         /* A fixing gets applied. Because the sufficient conditions for completeness are now possibly violated,
          * we need to add this permutation back to the permutation queue, if it's not already there.
          * Note: We do not need to check here explicitly for the conditions, because that is the first step when
          * the permutation is popped from the queue.
          */
         if ( !permsinqueue[k] )
         {
            permsqueue[(*permsqueuesize)++] = k;
            permsinqueue[k] = TRUE;
            assert( *permsqueuesize <= nperms );
         }
      }

      /* End of variable fixing event for fixing (i, fixing). */
   }

   return SCIP_OKAY;
}


/** Perform complete propagation for all individual symresack constraints induced by the chosen permutations in the group
 * When choosing @p basepow equal to 1 and @p support as all entries, the group generated by the permutation of the
 * constraint is considered. Otherwise, when different values are specified, the permutation is restricted to the
 * support and we use the permutation with @p basepow as exponent as generator.
 */
static
SCIP_RETCODE completeFixingsPerPermutation(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_SYMRETOPEGRAPH*  implgraph,          /**< implication graph data structure, initialized empty and for sufficient permutations. */
   SCIP_FIXINGQUEUE*     fixingqueue,        /**< Pointer to the fixing queue */
   int                   basepow,            /**< The exponent of the permutation generating the group that we want to compute the complete set of fixings for */
   int*                  support,            /**< The entries of the permutation that we are interested in. */
   int                   nsupport,           /**< The number of elements in the support */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< Pointer to the virtual fixings structure, or NULL. */
   SCIP_Bool             useproblembounds,   /**< whether or not the bounds of the problem must be used, in addition to the virtual fixings. */
   SCIP_Bool*            checkedentries,     /**< For each variable index, whether their value has been looked up */
   int*                  impactfulentries,   /**< Whether entries appear in any implication tree, or NULL */
   int*                  nimpactfulentries,  /**< Pointer to how many impactful entries there are, already. */
   SCIP_Bool*            entryisimpactful,   /**< Pointer to a boolean array specifying if an entry is impactful. */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_PERMUTATION* permutation;
   int nvars;
   int nperms;
   SCIP_SymretopeGraphNode* permgraphroots;
   SCIP_SymretopeGraphNode** permgraphleaves;
   SCIP_SymretopeGraphNode* permgraphs;
   SCIP_SymretopeGraphNode* root;
   SCIP_SymretopeGraphNode* node;
   SCIP_SymretopeGraphNode* leaf;
   SCIP_SymretopeGraphNode* succ;
   SCIP_SymretopeGraphNode* twin;
   SCIP_SymretopeGraphNode* permgraph;
   int permpow;
   int leafid;
   int fixing;
   int* var1fixes;
   int* var2fixes;
   int var1fix;
   int var2fix;
   int i;
   int i_;
   int j;
   int jj;
   int k;
   SCIP_Bool tightened;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );
   assert( implgraph != NULL );

   *ngen = 0;
   *infeasible = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nvars >= 0 );
   assert( consdata->vars != NULL );
   assert( consdata->permutation != NULL );
   assert( consdata->permutation->order >= 0 );
   nvars = consdata->nvars;

   /* avoid trivial problems */
   if ( nvars < 2 )
      return SCIP_OKAY;

   permutation = consdata->permutation;
   vars = consdata->vars;
   assert( vars != NULL );
   assert( permutation != NULL );
   assert( implgraph->permpows != NULL );

   /* Initialize powers that we want to evaluate */
   if ( support == NULL )
   {
      nperms = consdata->nperms;

      /* Nothing to propagate if there are no permutations. */
      if ( nperms <= 0 )
         return SCIP_OKAY;

      /* We must be able to store (at least) nperms powers in permpows. */
      assert( nperms >= implgraph->maxnperms );

      /* The powers that need to be evaluated. */
      for (k = 0; k < nperms; ++k)
      {
         implgraph->permpows[k] = k + 1;
         assert( implgraph->permpows[k] < permutation->order );
      }
   }
   else
   {
      /* If the support of the cycle is given, then we consider one cycle of size "nsupport".
       * The number of permutations is then the number of elements in the group generated by a cycle of length
       * "nsupport" raised to the power "basepow".
       */
      nperms = nsupport / gcd(nsupport, basepow) - 1;

      /* Nothing to propagate if there are no permutations. */
      if ( nperms <= 0 )
         return SCIP_OKAY;

      /* We must be able to store (at least) nperms powers in permpows. */
      assert( nperms <= implgraph->maxnperms );

      /* The powers that need to be evaluated. */
      for (k = 0; k < nperms; ++k)
      {
         implgraph->permpows[k] = (k + 1) * basepow;
         assert( implgraph->permpows[k] < permutation->order );
      }

      #ifndef NDEBUG
      /* Check that the support is monotonously increasing */
      for (k = 1; k < nsupport; ++k)
      {
         assert( support[k] > support[k - 1] );
      }
      #endif
   }

   /* Get implication tree data structures. */
   permgraphroots = implgraph->permgraphroots;
   permgraphleaves = implgraph->permgraphleaves;
   permgraphs = implgraph->permgraphs;
   assert( permgraphroots != NULL );
   assert( permgraphleaves != NULL );
   assert( permgraphs != NULL );

   /* Create the roots, connected to a single loose end. */
   /* @todo we can use 'clean' arrays for the root too, as in the end it has no successors. */
   for (k = 0; k < nperms; ++k)
   {
      root = &permgraphroots[k];
      root->nodetype = SYMRETOPE_ROOT;
      root->successor1 = NULL;
      root->successor2 = NULL;
      root->predecessor = NULL;
      #ifndef NDEBUG
      root->fixing = -k;   /* For debugging purposes, use the fixing-field to see what the permutation id is */
      /* NOTE: There is a little bit of memory spillage because we do not use "fixing" for the root node. */
      #endif

      /* The root is connected to a single leaf. */
      permgraphleaves[2 * k] = root;
      permgraphleaves[2 * k + 1] = NULL;
   }

   /* Initialize permutation indices. They start at 0. */
   /* Initialize queues: Schedule all permutations that are not the identity. */
   for (k = 0; k < nperms; ++k)
   {
      implgraph->permindices[k] = 0;
      implgraph->permsinqueue[k] = TRUE;
      implgraph->permsqueue[k] = k;
   }
   implgraph->permsqueuesize = nperms;

   /* Variable fixings: To store the fixings of the two variables before extending the leaves. */
   SCIP_CALL( SCIPallocBufferArray(scip, &var1fixes, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &var2fixes, 2) );

   while (implgraph->permsqueuesize > 0)
   {
      /* Pick a permutation that we still need to handle, and remove from queue. */
      k = implgraph->permsqueue[--implgraph->permsqueuesize];
      permpow = implgraph->permpows[k];
      implgraph->permsinqueue[k] = FALSE;

      /* Get the permutation, and the inverse of the permutation. */
      assert( k >= 0 && k < nperms );
      assert( permpow > 0 && permpow < consdata->permutation->order );
      assert( support == NULL || (permpow % nsupport > 0 ) );

      root = &permgraphroots[k];
      permgraph = &permgraphs[2 * nvars * k];

      /* apply events to current permutation until stopping criterion is met */
      while(TRUE)
      {
         SCIP_Bool conditionviolated = FALSE;

         /* Debugging: Now none of the root nodes could be connected to a necessary fixing node in the graphs. */
         #ifndef NDEBUG
         for (j = 0; j < nperms; ++j)
         {
            assert( (&permgraphroots[j])->successor1 == NULL
               || (&permgraphroots[j])->successor1->nodetype == SYMRETOPE_COND );
            assert( (&permgraphroots[j])->successor2 == NULL
               || (&permgraphroots[j])->successor2->nodetype == SYMRETOPE_COND );
            assert(
               ((&permgraphroots[j])->successor1 == NULL && (&permgraphroots[j])->successor2 == NULL)
               ? (
                  (permgraphleaves[2 * j] == (&permgraphroots[j]) || permgraphleaves[2 * j] == NULL)
                  && (permgraphleaves[2 * j + 1] == (&permgraphroots[j]) || permgraphleaves[2 * j + 1] == NULL)
               )
               : TRUE
            );
         }
         #endif

         /* get the variable index @p i from which we start to look for fixings for the current permutation
          *
          * We can stop if one of the following criteria is met:
          *
          * C1: there is no loose end in the implication tree
          * C2: the index lies outside the range of the support we need to consider
          * C3: there is a conditional fixing node on each rooted path to a leaf;
          *     (and at this position, that means that the root has a necessary fixing node as child.)
          *      - i is no 0-fixing, and invperm[i] is no 1-fixing;
          *      - perm[i] > i, invperm[i] > i.
          */

         /* Condition C1: there is no loose end in the implication tree. */
         if ( permgraphleaves[2 * k] == NULL && permgraphleaves[2 * k + 1] == NULL )
            break;

         /* Condition C2: the index lies outside the range of the support we need to consider */
         /* select @p i from the given support or recognize that condition C2 is violated */
         i_ = implgraph->permindices[k];
         if ( support == NULL )
         {
            i = i_;

            if ( i >= nvars )
               conditionviolated = TRUE;
         }
         else
         {
            if ( i_ >= nsupport )
               conditionviolated = TRUE;
            else
            {
               i = support[i_];
               assert( i < nvars );
            }
         }

         if ( conditionviolated )
            break;

         /* If i is a fixed point for this permutation, then we increase the index by one and repeat. */
         j = permGet(permutation, i, -permpow);
         assert( j >= 0 && j < nvars );
         if ( i == j )
         {
            /* This is a fixed point for this permutation. Can continue. */
            ++(implgraph->permindices[k]);
            continue;
         }

         /* For the remainder of the checks, the fixing values of entries i and j are of importance.
          * Hence, their value impacts the behaviour, and they should be marked as impactful. */
         if ( impactfulentries != NULL )
         {
            assert( nimpactfulentries != NULL );
            assert( entryisimpactful != NULL );

            if ( !entryisimpactful[i] )
            {
               impactfulentries[(*nimpactfulentries)++] = i;
               entryisimpactful[i] = TRUE;
               assert( *nimpactfulentries <= nvars );
            }
            if ( !entryisimpactful[j] )
            {
               impactfulentries[(*nimpactfulentries)++] = j;
               entryisimpactful[j] = TRUE;
               assert( *nimpactfulentries <= nvars );
            }
         }

         /* Condition C3:
          * * There is a conditional fixing node on each rooted path to a leaf;
          *   (And at this position, that means that the root has a necessary fixing node as child.)
          * * i is no 0-fixing, and invperm[i] is no 1-fixing;
          * * perm[i] > i, invperm[i] > i.
          * Use that j = invperm[i].
          */
         jj = permGet(permutation, i, permpow);
         if ( jj > i && j > i
            && getVarFixing(vars, i, virtualfixings, useproblembounds, checkedentries) != FIXED0
            && getVarFixing(vars, j, virtualfixings, useproblembounds, checkedentries) != FIXED1
            && (
                  (root->successor1 != NULL && root->successor1->nodetype == SYMRETOPE_COND)
                  || (root->successor2 != NULL && root->successor2->nodetype == SYMRETOPE_COND)
               )
            )
            break;

         /* The sufficient completeness conditions are not met, so apply steps for the index increasing event for k. */

         /* Get the fixings of variables i and j on the possible tree branches. */
         for (leafid = 0; leafid < 2; ++leafid)
         {
            leaf = permgraphleaves[2 * k + leafid];
            if ( leaf == NULL )
            {
               var1fixes[leafid] = -1;
               var2fixes[leafid] = -1;
               continue;
            }
            assert( leaf->successor1 == NULL );
            assert( leaf->successor2 == NULL );

            /* Get the value of var i */
            var1fix = getVarFixing(vars, i, virtualfixings, useproblembounds, checkedentries);
            if ( var1fix == UNFIXED )
            {
               /* We can check whether an internal node exists by checking whether it has a predecessor. */
               if ( (&permgraph[2 * i + leafid])->predecessor != NULL )
               {
                  /* Entry i is fixed on the path after (or possibly before) the branching. */
                  fixing = (&permgraph[2 * i + leafid])->fixing;
                  assert( fixing % nvars == i );

                  /* COMMENT: Use separate function to decouple the encoding of data from their content
                   * var1fix = getFixingValue(fixing, nvars);
                   */
                  var1fix = fixing > nvars ? FIXED1 : FIXED0;
               }
               else if ( (&permgraph[2 * i + (1 - leafid)])->predecessor != NULL )
               {
                  /* Entry i is definitively not fixed on the path before the branching. */
                  fixing = (&permgraph[2 * i + (1 - leafid)])->fixing;
                  assert( fixing % nvars == i );
                  var1fix = fixing > nvars ? FIXED1 : FIXED0;
               }
            }
            var1fixes[leafid] = var1fix;

            /* Get the value of var j */
            var2fix = getVarFixing(vars, j, virtualfixings, useproblembounds, checkedentries);
            if ( var2fix == UNFIXED )
            {
               /* We can check whether a internal node exists by checking whether it has a predecessor. */
               if ( (&permgraph[2 * j + leafid])->predecessor != NULL )
               {
                  /* Entry i is fixed on the path after (or possibly before) the branching. */
                  fixing = (&permgraph[2 * j + leafid])->fixing;
                  assert( fixing % nvars == j );
                  var2fix = fixing > nvars ? FIXED1 : FIXED0;
               }
               else if ( (&permgraph[2 * j + (1 - leafid)])->predecessor != NULL )
               {
                  /* Entry j is definitively not fixed on the path before the branching. */
                  fixing = (&permgraph[2 * j + (1 - leafid)])->fixing;
                  assert( fixing % nvars == j );
                  var2fix = fixing > nvars ? FIXED1 : FIXED0;
               }
            }
            var2fixes[leafid] = var2fix;
         }

         /* Now, extend the leaves. */
         for (leafid = 0; leafid < 2; ++leafid)
         {
            leaf = permgraphleaves[2 * k + leafid];
            if ( leaf == NULL )
               continue;
            assert( leaf->successor1 == NULL );
            assert( leaf->successor2 == NULL );

            var1fix = var1fixes[leafid];
            var2fix = var2fixes[leafid];
            assert( var1fix >= 0 );
            assert( var2fix >= 0 );

            switch (var1fix + FIXEDMAX * var2fix)
            {
               case FIXED0 + FIXEDMAX * FIXED0:   /* (0, 0) */
               case FIXED1 + FIXEDMAX * FIXED1:   /* (1, 1) */
                  /* Nothing to do here. */
                  break;
               case FIXED1 + FIXEDMAX * FIXED0:   /* (1, 0) */
                  /* Remove leaf. */
                  permgraphleaves[2 * k + leafid] = NULL;
                  break;
               case FIXED0 + FIXEDMAX * UNFIXED:  /* (0, _) */
                  /* Locally apply fixing of j to 0. */
                  node = &permgraph[2 * j + leafid];
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );
                  assert( node->nodetype == 0 );
                  assert( node->fixing == 0 );

                  /* Fix links */
                  node->predecessor = leaf;
                  node->nodetype = SYMRETOPE_NECC;
                  node->fixing = j; /* Encoding of fixing i to b is i + nvars * b. This is fixing j to 0. */
                  leaf->successor1 = node;
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  /* node becomes the new leaf. */
                  permgraphleaves[2 * k + leafid] = node;
                  break;
               case UNFIXED + FIXEDMAX * FIXED1:  /* (_, 1) */
                  /* Locally apply fixing of i to 1. */
                  node = &permgraph[2 * i + leafid];
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );
                  assert( node->nodetype == 0 );
                  assert( node->fixing == 0 );

                  /* Fix links */
                  node->predecessor = leaf;
                  node->nodetype = SYMRETOPE_NECC;
                  node->fixing = i + nvars;
                  leaf->successor1 = node;
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  /* node becomes the new leaf. */
                  permgraphleaves[2 * k + leafid] = node;
                  break;
               case FIXED1 + FIXEDMAX * UNFIXED:  /* (1, _) */
                  /* Locally apply conditional fixing of j to 1. */
                  node = &permgraph[2 * j + leafid];
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );
                  assert( node->nodetype == 0 );
                  assert( node->fixing == 0 );

                  /* Fix links */
                  node->predecessor = leaf;
                  node->nodetype = SYMRETOPE_COND;
                  node->fixing = j + nvars;
                  leaf->successor1 = node;
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  /* node becomes the new leaf. */
                  permgraphleaves[2 * k + leafid] = node;
                  break;
               case UNFIXED + FIXEDMAX * FIXED0:  /* (_, 0) */
                  /* Locally apply conditional fixing of i to 0. */
                  node = &permgraph[2 * i + leafid];
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );
                  assert( node->nodetype == 0 );
                  assert( node->fixing == 0 );

                  /* Fix links */
                  node->predecessor = leaf;
                  node->nodetype = SYMRETOPE_COND;
                  node->fixing = i;
                  leaf->successor1 = node;
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  /* node becomes the new leaf. */
                  permgraphleaves[2 * k + leafid] = node;
                  break;
               case FIXED0 + FIXEDMAX * FIXED1:   /* (0, 1) */
                  /* This path is infeasible. Therefore, the path ending in this leaf gets removed to its first
                   * conditional ancestor, and we apply the merging operation. */

                  /* This is no longer going to be a leaf, any more. */
                  permgraphleaves[2 * k + leafid] = NULL;

                  /* Remove all predecessors until a first conditional (or root) node. */
                  while (leaf->nodetype == SYMRETOPE_NECC)
                  {
                     /* Because leaf is a leaf of the tree, it cannot have successors. */
                     assert( leaf->successor1 == NULL );
                     assert( leaf->successor2 == NULL );

                     /* Now, remove leaf. */
                     node = leaf;
                     leaf = leaf->predecessor;

                     /* Because leaf is a necessary fixing node, its predecessor has exactly 1 child. */
                     assert( leaf->successor1 == node );
                     assert( leaf->successor2 == NULL );
                     leaf->successor1 = NULL;

                     /* Empty all the fields of the old leaf. */
                     node->fixing = 0;
                     node->nodetype = 0;
                     node->predecessor = NULL;
                     node->successor1 = NULL;
                     node->successor2 = NULL;
                  }

                  /* Now we are in a conditional or root node. */
                  assert( leaf->nodetype == SYMRETOPE_COND || leaf->nodetype == SYMRETOPE_ROOT );

                  if ( leaf->nodetype == SYMRETOPE_ROOT )
                  {
                     /* We have now found infeasibility! */
                     *infeasible = TRUE;

                     /* Conflict analysis: Infeasibility follows by comparing vars with perm[vars].
                      * Only do so when not peeking. */
                     if ( virtualfixings == NULL && SCIPisConflictAnalysisApplicable(scip) )
                     {
                        SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                        SCIP_CALL( resolveSymretopeConflictVariables(scip, NULL, SCIP_BOUNDTYPE_LOWER, vars, nvars,
                           permutation, permpow, NULL) );
                        SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
                     }

                     goto Cleanup;
                  }
                  else if ( leaf->nodetype == SYMRETOPE_COND )
                  {
                     /* Apply merging operation. */
                     node = leaf->predecessor;

                     /* Get the sibling of "leaf". */
                     twin = node->successor1;
                     if ( twin == leaf )
                        twin = node->successor2;
                     assert( twin != leaf );

                     if ( twin == NULL )
                     {
                        /* leaf has no sibling. Replace leaf by a necessary fixing node with opposing fixing. */
                        leaf->nodetype = SYMRETOPE_NECC;
                        if ( leaf->fixing >= nvars )
                           leaf->fixing = leaf->fixing - nvars;
                        else
                           leaf->fixing = leaf->fixing + nvars;
                        assert( leaf->fixing >= 0 && leaf->fixing < 2 * nvars );
                     }
                     else
                     {
                        /* leaf has sibling. Then we must have the following structure:
                         * "node" is the parent of "leaf" and "cond", each being conditional fixing nodes.
                         * "twin" has ONE child which is a necessary fixing node with the converse fixing of "leaf".
                         * "leaf" had (is removed) ONE child which is a necessary fixing node
                         *   with converse fixing of "twin".
                         *
                         * What we need to do now, is to remove "leaf",
                         * and pull twin->successor1->fixing one to the front.
                         */
                        assert( twin->nodetype == SYMRETOPE_COND );
                        assert( twin->successor1 != NULL );
                        assert( twin->successor1->nodetype == SYMRETOPE_NECC );
                        assert( twin->successor2 == NULL );
                        assert( abs( twin->successor1->fixing - leaf->fixing ) == nvars );

                        /* Remove leaf. */
                        leaf->fixing = 0;
                        leaf->nodetype = 0;
                        leaf->predecessor = NULL;
                        assert( leaf->successor1 == NULL );
                        assert( leaf->successor2 == NULL );

                        /* Overload the "leaf" variable with the successor of twin. Note: this is not not a leaf. */
                        leaf = twin->successor1;
                        /* This can only have one successor. */
                        assert( leaf->successor2 == NULL );

                        /* Make the new configuration: node -> leaf (the new one) -> twin -> successor of leaf. */
                        if ( leaf->successor1 != NULL )
                           leaf->successor1->predecessor = twin;
                        twin->successor1 = leaf->successor1;
                        twin->predecessor = leaf;
                        assert( twin->successor2 == NULL );

                        leaf->predecessor = node;
                        leaf->successor1 = twin;
                        assert( leaf->successor2 == NULL );

                        node->successor1 = leaf;
                        node->successor2 = NULL;

                        /* If the successor of "twin" (now called "leaf") was originally a leaf.
                         * Then now "twin" is a leaf.
                         */
                        if ( permgraphleaves[2 * k + 0] == leaf )
                           permgraphleaves[2 * k + 0] = twin;
                        if ( permgraphleaves[2 * k + 1] == leaf )
                           permgraphleaves[2 * k + 1] = twin;
                     }
                  }

                  break;
               case UNFIXED + FIXEDMAX * UNFIXED: /* (_, _) */
                  /* Now we create two branches. Starting from this leaf, this has two children that are each
                   * a conditional fixing: (i, 0) and (j, 1). These are then linked to neccesary fixings, resp.
                   * (j, 0) and (i, 1), and these are connected to leaves.
                   */

                  /* Check: If this happens we can only be in the situation where we have exactly one leaf (this one) */
                  assert( permgraphleaves[2 * k + (1 - leafid)] == NULL );

                  #ifndef NDEBUG
                  {
                     int dbg;
                     int dbginv;
                     if ( i > 0 )
                     {
                        SCIPdebugMessage("Found (_, _). Checking for permutation %d with power %d at index %d\n", k, permpow, i);
                     }
                     for (dbg = 0; dbg < i; ++dbg)
                     {
                        dbginv = permGet(permutation, dbg, -permpow);
                        /* Ignore fixed points */
                        if ( dbg == dbginv )
                           continue;
                        /* var dbg is fixed, or dbg is somewhere in the tree. */
                        assert(
                           (getVarFixing(vars, dbg, virtualfixings, useproblembounds, checkedentries) != UNFIXED)
                           || ((&permgraph[2 * dbg + leafid])->predecessor != NULL)
                           || ((&permgraph[2 * dbg + (1 - leafid)])->predecessor != NULL)
                        );
                        /* var invperm[dbg] is fixed, or invperm[dbg] is somewhere in the tree. */
                        assert(
                           (getVarFixing(vars, dbginv, virtualfixings, useproblembounds, checkedentries) != UNFIXED)
                           || ((&permgraph[2 * dbginv + leafid])->predecessor != NULL)
                           || ((&permgraph[2 * dbginv + (1 - leafid)])->predecessor != NULL)
                        );
                     }
                  }
                  #endif

                  /* Create two branches.
                   * In one branch we have (i, 0) [COND] -> (j, 0) [NECC] -> leaf0;
                   * In one branch we have (j, 1) [COND] -> (i, 1) [NECC] -> leaf1;
                   */

                  /* First branch. */
                  node = &permgraph[2*i];
                  assert( node->fixing == 0 );
                  assert( node->nodetype == 0 );
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );

                  succ = &permgraph[2*j];
                  assert( succ->fixing == 0 );
                  assert( succ->nodetype == 0 );
                  assert( succ->predecessor == NULL );
                  assert( succ->successor1 == NULL );
                  assert( succ->successor2 == NULL );

                  leaf->successor1 = node;
                  node->predecessor = leaf;
                  node->successor1 = succ;
                  node->nodetype = SYMRETOPE_COND;
                  node->fixing = i; /* Fixing i to 0 */
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  succ->predecessor = node;
                  succ->nodetype = SYMRETOPE_NECC;
                  succ->fixing = j; /* Fixing j to 0 */
                  permgraphleaves[2 * k] = succ;
                  assert( succ->fixing >= 0 && succ->fixing < 2 * nvars );

                  /* Second branch */
                  node = &permgraph[2*j + 1];
                  assert( node->fixing == 0 );
                  assert( node->nodetype == 0 );
                  assert( node->predecessor == NULL );
                  assert( node->successor1 == NULL );
                  assert( node->successor2 == NULL );

                  succ = &permgraph[2*i + 1];
                  assert( succ->fixing == 0 );
                  assert( succ->nodetype == 0 );
                  assert( succ->predecessor == NULL );
                  assert( succ->successor1 == NULL );
                  assert( succ->successor2 == NULL );

                  leaf->successor2 = node;
                  node->predecessor = leaf;
                  node->successor1 = succ;
                  node->nodetype = SYMRETOPE_COND;
                  node->fixing = j + nvars; /* Fixing j to 1 */
                  assert( node->fixing >= 0 && node->fixing < 2 * nvars );

                  succ->predecessor = node;
                  succ->nodetype = SYMRETOPE_NECC;
                  succ->fixing = i + nvars; /* Fixing i to 1 */
                  permgraphleaves[2 * k + 1] = succ;
                  assert( succ->fixing >= 0 && succ->fixing < 2 * nvars );

                  /* We've yet created the other leaf node, which was NULL before.
                   * We should not run this procedure for that node too, so jump out of the outer for-loop.
                   */
                  goto Endofleaves;
               default:
                  /* All situations should be covered above. This should not happen. */
                  assert( FALSE );
            }
         }
         Endofleaves:

         /* Finally, increase i by one in preparation for the next iteration. */
         ++(implgraph->permindices[k]);

         /* If the root has a child that's a necessary fixing node,
          * then we should add that fixing to the set of fixings.
          */
         node = root->successor1;
         assert( node == NULL || node->nodetype != SYMRETOPE_NECC || root->successor2 == NULL );
         if ( node != NULL && node->nodetype == SYMRETOPE_NECC )
         {
            /* A fixing must be applied! */
            SCIP_CALL( enqueueFixing(scip, vars, cons, node->fixing, permutation, implgraph->permpows[k],
                  nvars, fixingqueue, virtualfixings, infeasible) );
            if ( *infeasible )
               goto Cleanup;
         }

         /* Apply variable fixing event if needed. */
         SCIP_CALL( applyFixings(scip, cons, vars, virtualfixings, permutation, permgraphs, permgraphleaves,
               implgraph->permpows,
               nvars, nperms,
               fixingqueue,
               implgraph->permsinqueue, implgraph->permsqueue, &implgraph->permsqueuesize,
               ngen, infeasible, &tightened) );
         if ( *infeasible )
            goto Cleanup;

         /* Debugging: Now none of the root nodes could be connected to a necessary fixing node in the graphs. */
         #ifndef NDEBUG
         for (j = 0; j < nperms; ++j)
         {
            assert( (&permgraphroots[j])->successor1 == NULL
               || (&permgraphroots[j])->successor1->nodetype == SYMRETOPE_COND );
            assert( (&permgraphroots[j])->successor2 == NULL
               || (&permgraphroots[j])->successor2->nodetype == SYMRETOPE_COND );
            assert(
               ((&permgraphroots[j])->successor1 == NULL && (&permgraphroots[j])->successor2 == NULL)
               ? (
                  (permgraphleaves[2 * j] == (&permgraphroots[j]) || permgraphleaves[2 * j] == NULL)
                  && (permgraphleaves[2 * j + 1] == (&permgraphroots[j]) || permgraphleaves[2 * j + 1] == NULL)
               )
               : TRUE
            );
         }
         #endif

         /* Debugging: Check graphs for inconsistencies. */
         #ifndef NDEBUG
         {
            int dk;
            SCIP_SymretopeGraphNode* droot;
            SCIP_SymretopeGraphNode* br1;
            SCIP_SymretopeGraphNode* br2;
            for (dk = 0; dk < nperms; ++dk)
            {
               /* Go to the splitting point */
               droot = &permgraphroots[dk];
               while( droot->successor1 != NULL && droot->successor2 == NULL )
                  droot = droot->successor1;

               assert( !(droot->successor1 == NULL && droot->successor2 != NULL) );
               br1 = droot->successor1;
               br2 = droot->successor2;

               /* Graph does not split */
               if ( br1 == NULL || br2 == NULL )
                  continue;

               assert( br1->successor1 != NULL );
               assert( br2->successor1 != NULL );
               assert( br1->successor2 == NULL );
               assert( br2->successor2 == NULL );
               assert( abs(br1->successor1->fixing - br2->fixing) == nvars );
               assert( abs(br1->fixing - br2->successor1->fixing) == nvars );

               br1 = br1->successor1->successor1;
               br2 = br2->successor1->successor1;
               while (br1 != NULL && br2 != NULL)
               {
                  assert( (br1->fixing % nvars) == (br2->fixing % nvars) );
                  assert( br1->successor2 == NULL );
                  assert( br2->successor2 == NULL );
                  br1 = br1->successor1;
                  br2 = br2->successor1;
               }
               if ( br1 != NULL || br2 != NULL )
               {
                  assert( permgraphleaves[2 * dk] == NULL || permgraphleaves[2 * dk + 1] == NULL );
               }
            }
         }
         #endif
      }
      /* End of index increasing event for permutation k. */
   }

   Cleanup:
   /* Free memory: Variable fixings. */
   SCIPfreeBufferArray(scip, &var2fixes);
   SCIPfreeBufferArray(scip, &var1fixes);

   /* Clear the fixings queue. */
   while( fixingqueue->fixingqueuesize > 0 )
   {
      int fixingvalue;
      int fixingvarid;
      /* Get the fixing encoding, which is i + nvars * b, where b \in {0, 1} */
      fixingvalue = fixingqueue->fixingqueue[--(fixingqueue->fixingqueuesize)];
      assert( fixingvalue >= 0 && fixingvalue < 2 * nvars );
      /* Extract fixing variable index (i) and the associated fixing value (which is b in the description above). */
      fixingvarid = fixingvalue % nvars;
      fixingvalue = fixingvalue >= nvars ? FIXED1 : FIXED0;
      fixingqueue->fixinginqueue[fixingvarid] &= ~fixingvalue;
   }
   #ifndef NDEBUG
   /* Sanity check: The fixinginqueue must indeed be clean. */
   for (i = 0; i < nvars; ++i)
   {
      assert( fixingqueue->fixinginqueue[i] == 0 );
   }
   #endif
   assert( fixingqueue->fixingqueuesize == 0 );

   /* Clear the implication graphs. */
   for (k = 0; k < nperms; ++k)
   {
      root = &permgraphroots[k];

      SCIP_CALL( removeSubtree(scip, root, &permgraphleaves[2 * k]) );
      assert( root->predecessor == NULL );
      assert( root->successor1 == NULL );
      assert( root->successor2 == NULL );
   }
   #ifndef NDEBUG
   /* Sanity check: Make sure this memory is completely nulled by now. */
   for (i = 0; i < 2 * nvars * nperms; ++i)
   {
      node = &permgraphs[i];
      assert( node->fixing == 0 );
      assert( node->nodetype == 0 );
      assert( node->predecessor == NULL );
      assert( node->successor1 == NULL );
      assert( node->successor2 == NULL );
   }
   #endif

   return SCIP_OKAY;
}

/** The propagation function if the generating permutation is monotone and ordered. */
static
SCIP_RETCODE propVariablesMonotoneOrderedHotstart(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< NULL if fixings have to be applied directly. Otherwise a virtual fixings structure */
   SCIP_Bool             useproblembounds,   /**< whether or not the bounds of the problem must be used, in addition to the virtual fixings. */
   SCIP_Bool*            checkedentries,     /**< For each variable index, whether their value has been looked up */
   SCIP_Bool             findcompleteset,    /**< whether the complete set of fixings is sought after, or feasibility only */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen,               /**< pointer to store number of generated bound strengthenings */
   int                   eqpow,              /**< The power of the permutation that generates the strict equality-subgroup (in pseudocode: \mu) */
   int                   startcolour,        /**< The colour to start from. */
   SCIP_SYMRETOPEGRAPH*  implgraph,          /**< Pointer to the implication graph data structure */
   SCIP_FIXINGQUEUE*     fixingqueue,        /**< Pointer to the fixing queue data structure */
   SCIP_CONSDATA*        consdata,           /**< Constraint data of cons */
   SCIP_PERMUTATION*     permutation         /**< The permutation object that generates the group */
)
{
   int c;
   int i;
   int i_;

   /* Since completeFixingsPerPermutation resets ngen, we store the number of fixings in `newngen`.
    * Then, if new fixings are found, we increase `ngen` by this number. */
   int newngen;

   /* Marking for which entries peeking is necessary. */
   int* impactfulentries;
   int nimpactfulentries;
   SCIP_Bool* entryisimpactful;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nvars > 0 );

   assert( permutation != NULL );
   assert( permutation->ismonotone );
   assert( permutation->isordered );
   assert( permutation->maxcyclesize >= 1 );

   /* If we are peeking, then we want to know which entries to peek on. */
   if ( virtualfixings == NULL && findcompleteset )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &impactfulentries, consdata->nvars) );
      nimpactfulentries = 0;
      SCIP_CALL( SCIPallocClearBufferArray(scip, &entryisimpactful, consdata->nvars) );
   }
   else
      impactfulentries = NULL;

   for (c = startcolour; c < permutation->ncycles; ++c)
   {
      int* cycle;
      int cyclen;
      int* subcyclevalues;
      SCIP_Bool unfixedExists;

      /* If the eqpow is the order of the permutation, then only the identity remains and we can stop. */
      if ( eqpow == permutation->order )
         break;

      /* It is not possible that the generating power becomes larger than the order. */
      assert( eqpow < permutation->order );

      /* Get the subcycle and length of subcycle c. */
      cycle = permutation->cycles[c];
      cyclen = permutation->cyclelengths[c];

      /* Skip cycles for which the subcycle is trivial to the power eqpow. */
      if ( eqpow % cyclen == 0 )
      {
         continue;
      }

      /* Determine the complete set of fixings for all symresacks in the group generated by cycle c to the power mu */
      /* Note: The conflict analysis remains correct. This conflict analysis is namely executed on the whole vector
       * rather than the vector restricted to the support. Because (by the inductive argument) the powers that we have
       * selected find equality for the previous cycles. This will also be found by the conflict analysis.
       */
      SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, eqpow, cycle, cyclen, virtualfixings,
         useproblembounds, checkedentries, impactfulentries, &nimpactfulentries, entryisimpactful, infeasible,
         &newngen) );
      *ngen += newngen;

      /* Infeasibility detected. Can stop here. */
      if ( *infeasible )
         break;

      /* Peek-check for additional fixings. */
      if ( virtualfixings == NULL && findcompleteset )
      {
         SCIP_Bool peekinfeasible;
         SCIP_Bool tightened;
         int minunfixedentryinfirsthalfofcycle;
         SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixingspeek;
         int virtualngen;

         SCIP_CALL( allocVirtualFixings(scip, &virtualfixingspeek, consdata->nvars) );

         /* Get minimal unfixed entry in cycle.
          * We have shown that for a cycle (1, ..., n) there exists a vector X with X > perm(X),
          *    (where > is the lexicographic order)
          * Constructively, two variants:
          * 1. If the first unfixed entry "k" is in the first half of the vector:
          *    Set all unfixed entries except for "k" to 0, and "k" to 1, then X > perm(X).
          * 2. If the first unfixed entry "k" is in the second half of the vector:
          *    Set all unfixed entries to 0.
          * Thus, we do not have to peek for these values, saving us some time.
          *
          * Compute the first unfixed entry in the first half of the cycle.
          */
         minunfixedentryinfirsthalfofcycle = -1;
         for (i_ = 0; i_ < cyclen / 2; ++i_)
         {
            i = cycle[i_];
            if ( getVarFixing(consdata->vars, i, virtualfixings, useproblembounds, NULL) == UNFIXED )
            {
               minunfixedentryinfirsthalfofcycle = i;
               break;
            }
         }

         peekinfeasible = FALSE;
         tightened = FALSE;
         // printf("Number of impactful entries: %i / %i \n", nimpactfulentries, consdata->nvars);
         while ( nimpactfulentries > 0 )
         {
            /* Get the entry for which we would like to peek. */
            i = impactfulentries[--nimpactfulentries];

            /* This entry must be in the current cycle. The cycle is monotone, so the minimal entry is known. */
            assert( i >= cycle[0] );
            assert( i <= cycle[cyclen - 1] );

            if ( tightened )
            {
               /* COMMENT: We need to discuss this */
               /* @todo: Why is this needed? It costs an additional O(n^2) but could possibly find additional fixings. */
               /* Compute fixings until the set of fixings is complete for all permutations in the group.
                * Also add new impactfulentries, only these that haven't been found before. */
               /* It is not needed for completeness of the propagation algorithm, but it saves a lot of variable
                * fixings with inferinfo = -1, causing less effective RESPROP calls and leading to more running time. */
               SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, eqpow, cycle, cyclen,
                  virtualfixings, useproblembounds, checkedentries, impactfulentries, &nimpactfulentries,
                  entryisimpactful, infeasible, &newngen) );
               *ngen += newngen;

               /* If infeasibility is found, then we can stop here. */
               if ( *infeasible )
                  break;
            }
            tightened = FALSE;

            /* If entry i is fixed, do not need to peek. */
            if ( getVarFixing(consdata->vars, i, virtualfixings, useproblembounds, checkedentries) != UNFIXED )
               continue;

            /* If i is the first unfixed entry and in the first half of the vector, fixing to 1 is possible.
             * Check: what if we fix to 0?
             */
            if ( i == minunfixedentryinfirsthalfofcycle )
            {
               /* What if variable "i" is 0? */
               /* COMMENT: This causes an O(nvars) overhead in every iteration.
                * If we incorporate this in the functions that fix a variable, we don't have additional overhead.
                */
               if ( virtualfixings == NULL )
                  clearVirtualFixings(virtualfixingspeek);
               else
                  copyVirtualFixings(virtualfixings, virtualfixingspeek);
               assert( getVirtualFixing(virtualfixingspeek, i) == UNFIXED );
               setVirtualFixing(virtualfixingspeek, i, FIXED0);
               SCIP_CALL( propVariablesMonotoneOrderedHotstart(scip, cons, virtualfixingspeek, useproblembounds,
                  checkedentries, FALSE, &peekinfeasible, &virtualngen, eqpow, c, implgraph, fixingqueue, consdata,
                  permutation) );
               if ( peekinfeasible )
               {
                  /* Zero-fixing of "i" is not possible. Fix to 1. */
                  SCIP_CALL( setVarFixing(scip, cons, consdata->vars, i, virtualfixings, FIXED1, infeasible, &tightened, -1) );
                  if ( *infeasible )
                     goto CleanupPeek;
                  if ( tightened )
                     ++(*ngen);

                  continue;
               }
            }
            /* If i is not the first unfixed entry or not in the first half of the vector, fixing to 0 is possible.
             * Check: what if we fix to 1?
             */
            else
            {
               /* What if variable "i" is 1? */
               if ( virtualfixings == NULL )
                  clearVirtualFixings(virtualfixingspeek);
               else
                  copyVirtualFixings(virtualfixings, virtualfixingspeek);
               assert( getVirtualFixing(virtualfixingspeek, i) == UNFIXED );
               setVirtualFixing(virtualfixingspeek, i, FIXED1);

               SCIP_CALL( propVariablesMonotoneOrderedHotstart(scip, cons, virtualfixingspeek, useproblembounds,
                  checkedentries, FALSE, &peekinfeasible, &virtualngen, eqpow, c, implgraph, fixingqueue, consdata,
                  permutation) );
               if ( peekinfeasible )
               {
                  /* One-fixing of "i" is not possible. Fix to 0. */
                  SCIP_CALL( setVarFixing(scip, cons, consdata->vars, i, virtualfixings, FIXED0, infeasible, &tightened, -1) );
                  if ( *infeasible )
                     goto CleanupPeek;
                  if ( tightened )
                     ++(*ngen);

                  continue;
               }
            }
         }

         CleanupPeek:
         freeVirtualFixings(scip, &virtualfixingspeek);
      }

      /* Update eqpow */
      unfixedExists = FALSE;
      SCIP_CALL( SCIPallocBufferArray(scip, &subcyclevalues, cyclen) );

      /* Determine if there are unfixed entries.
       * Additionally, fill "subcyclevalues" with 0- or 1-entries associated to the fixings.
       */
      for (i_ = 0; i_ < cyclen && !unfixedExists; ++i_)
      {
         i = cycle[i_];
         switch (getVarFixing(consdata->vars, i, virtualfixings, useproblembounds, checkedentries))
         {
            case UNFIXED:
               unfixedExists = TRUE;
               break;
            case FIXED0:
               subcyclevalues[i_] = 0;
               break;
            case FIXED1:
               subcyclevalues[i_] = 1;
               break;
            default:
               assert( FALSE );
         }
      }

      /* If an unfixed entry exists,
       * there is a vector such that x \succ \gamma(x) for \gamma being some power of this cycle (!= identity)
       */
      if ( unfixedExists )
         eqpow = lcm(eqpow, cyclen);
      /* If no unfixed entry exists, then everything is fixed. Determine the minimal power k such that
       * \gamma^k(x) = x.
       */
      else
      {
         int k;
         for (k = 1; k < cyclen; ++k)
         {
            /* Check if entry i and i + k of subcyclevalues have the same value. */
            for (i = 0; i < cyclen; ++i)
            {
               if ( subcyclevalues[i] != subcyclevalues[(i + k) % cyclen] )
                  break;
            }
            /* If the loop did break, then i < cyclen and entry i and i + k (mod cyclen) are different entries
             * in the subcycle. So this value of k does not qualify.
             * If the loop did not break, then i = cyclen and the vector and vector shifted by k places is identical
             */
            if ( i == cyclen )
               break;
         }
         eqpow = lcm(eqpow, k);
      }
      SCIPfreeBufferArray(scip, &subcyclevalues);
   }

   /* Cleaning up */
   if ( virtualfixings == NULL && findcompleteset )
   {
      assert( nimpactfulentries == 0 || *infeasible );
      SCIPfreeBufferArray(scip, &impactfulentries);
      SCIPfreeBufferArray(scip, &entryisimpactful);
   }

   return SCIP_OKAY;
}


/** The propagation function if the generating permutation is monotone and ordered. */
static
SCIP_RETCODE propVariablesMonotoneOrdered(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< NULL if fixings have to be applied directly. Otherwise a virtual fixings structure */
   SCIP_Bool             useproblembounds,   /**< whether or not the bounds of the problem must be used, in addition to the virtual fixings. */
   SCIP_Bool*            checkedentries,     /**< For each variable index, whether their value has been looked up */
   SCIP_Bool             findcompleteset,    /**< whether the complete set of fixings is sought after, or feasibility only */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
)
{
   SCIP_SYMRETOPEGRAPH* implgraph;
   SCIP_FIXINGQUEUE* fixingqueue;
   SCIP_CONSDATA* consdata;
   SCIP_PERMUTATION* permutation;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nvars > 0 );

   permutation = consdata->permutation;
   assert( permutation != NULL );
   assert( permutation->ismonotone );
   assert( permutation->isordered );
   assert( permutation->maxcyclesize >= 1 );

   SCIP_CALL( allocSymretopeGraph(scip, &implgraph, consdata->nvars, permutation->maxcyclesize - 1) );
   SCIP_CALL( allocFixingQueue(scip, &fixingqueue, consdata->nvars) );

   SCIP_CALL( propVariablesMonotoneOrderedHotstart(scip, cons, virtualfixings, useproblembounds, checkedentries,
      findcompleteset, infeasible, ngen, 1, 0, implgraph, fixingqueue, consdata, permutation) );

   freeFixingQueue(scip, &fixingqueue);
   freeSymretopeGraph(scip, &implgraph);
   return SCIP_OKAY;
}

/** The propagation function that works for general permutations, but may not always find the complete set of fixings. */
static
SCIP_RETCODE propVariablesStandard(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< structure of virtual fixings, or NULL if fixings are to be applied globally. */
   SCIP_Bool             useproblembounds,   /**< whether or not the bounds of the problem must be used, in addition to the virtual fixings. */
   SCIP_Bool*            checkedentries,     /**< For each variable index, whether their value has been looked up */
   SCIP_Bool             dopeek,             /**< whether the complete set of fixings is sought after, or feasibility only */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_SYMRETOPEGRAPH* implgraph;
   SCIP_FIXINGQUEUE* fixingqueue;
   int* impactfulentries;
   int nimpactfulentries;
   int newngen;
   SCIP_Bool* entryisimpactful;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   #ifdef SCIP_DEBUG
   ++consdata->debugcnt;
   SCIPdebugMsg(scip, "Propagating variables of constraint <%s>; (%d).\n", SCIPconsGetName(cons), consdata->debugcnt);
   #endif

   /* If we are peeking, then we want to know which entries are effective to peek on. */
   if ( virtualfixings == NULL && dopeek )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &impactfulentries, consdata->nvars) );
      nimpactfulentries = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &entryisimpactful, consdata->nvars) );
   }
   else
      impactfulentries = NULL;

   /* Initialize data structure for permutation graph */
   SCIP_CALL( allocSymretopeGraph(scip, &implgraph, consdata->nvars, consdata->nperms) );

   /* And for the fixing queue, that we will recycle for various calls to the propagator */
   SCIP_CALL( allocFixingQueue(scip, &fixingqueue, consdata->nvars) );

   /* First, compute fixings until the set of fixings is complete for all permutations in the group */
   SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, 1, NULL, -1, virtualfixings,
      useproblembounds, checkedentries, impactfulentries, &nimpactfulentries, entryisimpactful, infeasible, &newngen) );
   *ngen += newngen;

   /* If infeasibility is found, then we can stop here. */
   if ( *infeasible )
   {
      freeFixingQueue(scip, &fixingqueue);
      freeSymretopeGraph(scip, &implgraph);

      if ( virtualfixings == NULL && dopeek )
      {
         SCIPfreeBufferArray(scip, &impactfulentries);
         SCIPfreeBufferArray(scip, &entryisimpactful);
      }
      return SCIP_OKAY;
   }

   /* If we do not allow peek, stop here. */
   if ( virtualfixings == NULL && dopeek )
   {
      /* Now peek. For each non-fixed entry, test if feasibility is found if we fix to the opposite value, respectively. */
      int i;
      SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixingspeek;
      int virtualngen;
      SCIP_Bool peekinfeasible;
      SCIP_Bool tightened;

      assert( consdata->vars != NULL );
      assert( consdata->nvars > 0 );
      tightened = FALSE;
      SCIP_CALL( allocVirtualFixings(scip, &virtualfixingspeek, consdata->nvars) );
      // printf("Number of impactful entries: %i / %i \n", nimpactfulentries, consdata->nvars);
      while ( nimpactfulentries > 0 )
      {
         /* Get the entry for which we would like to peek. */
         i = impactfulentries[--nimpactfulentries];
         assert( entryisimpactful[i] );

         if ( tightened )
         {
            /* Compute fixings until the set of fixings is complete for all permutations in the group, again.
             * We add new (not earlier encountered) entries to impactfulentries, if we happen to find them. */
            SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, 1, NULL, -1, virtualfixings,
               useproblembounds, checkedentries, impactfulentries, &nimpactfulentries, entryisimpactful,
               infeasible, &newngen) );
            *ngen += newngen;
            /* If infeasibility is found, then we can stop here. */
            if ( *infeasible )
               goto Cleanup;
         }
         tightened = FALSE;

         if ( getVarFixing(consdata->vars, i, NULL, useproblembounds, NULL) != UNFIXED )
            continue;

         /* What if variable "i" is 0? */
         clearVirtualFixings(virtualfixingspeek);
         setVirtualFixing(virtualfixingspeek, i, FIXED0);
         SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, 1, NULL, -1, virtualfixingspeek,
            useproblembounds, checkedentries, NULL, NULL, NULL, &peekinfeasible, &virtualngen) );
         if ( peekinfeasible )
         {
            /* Zero-fixing of "i" is not possible. Fix to 1. */
            SCIP_CALL( setVarFixing(scip, cons, consdata->vars, i, virtualfixings, FIXED1, infeasible,
               &tightened, -1) );
            if ( *infeasible )
               goto Cleanup;
            if ( tightened )
               ++(*ngen);

            continue;
         }

         /* What if variable "i" is 1? */
         clearVirtualFixings(virtualfixingspeek);
         setVirtualFixing(virtualfixingspeek, i, FIXED1);
         SCIP_CALL( completeFixingsPerPermutation(scip, cons, implgraph, fixingqueue, 1, NULL, -1, virtualfixingspeek,
            useproblembounds, checkedentries, NULL, NULL, NULL, &peekinfeasible, &virtualngen) );
         if ( peekinfeasible )
         {
            /* One-fixing of "i" is not possible. Fix to 0. */
            SCIP_CALL( setVarFixing(scip, cons, consdata->vars, i, virtualfixings, FIXED0, infeasible,
               &tightened, -1) );
            if ( *infeasible )
               goto Cleanup;
            if ( tightened )
               ++(*ngen);

            continue;
         }
      }

      Cleanup:
      freeVirtualFixings(scip, &virtualfixingspeek);

      assert( nimpactfulentries == 0 || *infeasible );
      SCIPfreeBufferArray(scip, &impactfulentries);
      SCIPfreeBufferArray(scip, &entryisimpactful);
   }

   freeFixingQueue(scip, &fixingqueue);
   freeSymretopeGraph(scip, &implgraph);

   return SCIP_OKAY;
}

/** The propagation function for the symretope constraint handler */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings, /**< virtual fixings structure, or NULL if fixings are to be applied globally. */
   SCIP_Bool             useproblembounds,   /**< whether or not the bounds of the problem must be used, in addition to the virtual fixings. */
   SCIP_Bool*            checkedentries,     /**< For each variable index, whether their value has been looked up */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool findcompleteset;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->permutation != NULL );

   conshdlr = SCIPconsGetHdlr(cons);
   assert( conshdlr != NULL );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->permutation != NULL );

   findcompleteset = conshdlrdata->symretopepeek &&
      (SCIPinProbing(scip) ? conshdlrdata->probingpeek : TRUE);

   if ( consdata->permutation->ismonotone && consdata->permutation->isordered )
   {
      // printf("Yep it's monotone and ordered!\n");
      SCIP_CALL( propVariablesMonotoneOrdered(scip, cons, virtualfixings, useproblembounds, checkedentries,
         findcompleteset, infeasible, ngen) );
   }
   else
   {
      // printf("Not monotone and ordered!\n");
      SCIP_CALL( propVariablesStandard(scip, cons, virtualfixings, useproblembounds, checkedentries,
         findcompleteset, infeasible, ngen) );
   }

   return SCIP_OKAY;
}


/** add symresack cover inequality */
static
SCIP_RETCODE addSymresackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variables */
   int*                  coeffs,             /**< coefficient vector of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;
#ifdef SCIP_DEBUG
   SCIP_CONSDATA* consdata;
   char name[SCIP_MAXSTRLEN];
#endif

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars > 0 );
   assert( vars != NULL );
   assert( coeffs != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

#ifdef SCIP_DEBUG
   consdata = SCIPconsGetData(cons);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_cover_%s_%d", SCIPconsGetName(cons), consdata->debugcnt);
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   ++consdata->debugcnt;
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
#endif
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < nvars; ++i)
   {
      if ( coeffs[i] == 1 || coeffs[i] == -1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], (SCIP_Real) coeffs[i]) );
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** Maximize a linear function on a "strict" symresack,
 *  that is a symresack where we do not allow the solution x = gamma(x).
 */
static
SCIP_RETCODE maximizeObjectiveSymresackStrict(
   SCIP*                scip,                /**< SCIP pointer */
   int                  nvars,               /**< number of variables in symresack */
   SCIP_Real*           objective,           /**< the objective vector */
   int*                 perm,                /**< the permutation (without fixed points) as an array */
   int*                 invperm,             /**< the inverse permutation as an array */
   int*                 maxcrit,             /**< pointer to the critical entry where optimality is found at */
   SCIP_Real*           maxsoluval           /**< pointer to store the optimal objective value */
)
{
   /* The maximal objective in every iteration. */
   SCIP_Real tmpobj;
   /* The new value of componentobj when combining two components. */
   SCIP_Real tmpnewcompobj;
   /* helperobj is the sum of all positive objective-sums for all components. */
   SCIP_Real helperobj = 0.0;

   int crit;
   int critinv;
   int i;

   /* For every vertex of degree < 2 we maintain componentends and componentobj. */
   int* componentends;
   SCIP_Real* componentobj;

   assert( scip != NULL );
   assert( nvars > 0 );
   assert( objective != NULL );
   assert( perm != NULL );
   assert( invperm != NULL );
   assert( maxcrit != NULL );
   assert( maxsoluval != NULL );

   /* The current best known critical entry and objective */
   *maxcrit = -1;
   *maxsoluval = -SCIP_DEFAULT_INFINITY;

   SCIP_CALL( SCIPallocBufferArray(scip, &componentends, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &componentobj, nvars) );

   /* Initialization: Every entry is a component in the graph,
    * having the corresponding objective
    */
   for (i = 0; i < nvars; ++i)
   {
      componentends[i] = i;
      componentobj[i] = objective[i];
      if ( SCIPisGT(scip, objective[i], 0.0) )
         helperobj += objective[i];
   }

   /* Iterate over all critical rows, and of the graph maintain the components on the vertices of degree < 2. */
   for (crit = 0; crit < nvars; ++crit)
   {
      critinv = invperm[crit];

      /* This is a fixed point. */
      if ( crit == critinv )
         continue;

      /* If the other end of the component of crit is critinv, then crit cannot be a critical entry. */
      if ( componentends[crit] == critinv )
         continue;

      /* Compute objective for crit as critical entry. Update if it is better than the best found objective */
      tmpobj = helperobj;
      if ( SCIPisLT(scip, componentobj[crit], 0.0) )
         tmpobj += componentobj[crit];
      if ( SCIPisGT(scip, componentobj[critinv], 0.0) )
         tmpobj -= componentobj[critinv];
      if ( SCIPisGT(scip, tmpobj, *maxsoluval) )
      {
         *maxsoluval = tmpobj;
         *maxcrit = crit;
      }

      /* Update helperobj */
      tmpnewcompobj = componentobj[crit] + componentobj[critinv];
      if ( SCIPisGT(scip, componentobj[crit], 0.0) )
         helperobj -= componentobj[crit];
      if ( SCIPisGT(scip, componentobj[critinv], 0.0) )
         helperobj -= componentobj[critinv];
      if ( SCIPisGT(scip, tmpnewcompobj, 0.0) )
         helperobj += tmpnewcompobj;

      /* Update the objective of a component */
      componentobj[componentends[crit]] = tmpnewcompobj;
      componentobj[componentends[critinv]] = tmpnewcompobj;

      /* Connect the endpoints of the newly created path */
      if ( componentends[crit] == crit )
      {
         componentends[crit] = componentends[critinv];
         componentends[componentends[critinv]] = crit;
      }
      else
      {
         componentends[componentends[crit]] = componentends[critinv];
         componentends[componentends[critinv]] = componentends[crit];
      }

      /* Early termination criterion. helperobj is upper bound to tmpobj for every next iteration,
       * so if helperobj <= maxsoluval then we can terminate earlier.
       */
      if ( SCIPisGE(scip, *maxsoluval, helperobj) )
         break;
   }

   /* It is always possible to make the first non-fixed entry critical. */
   assert( *maxcrit >= 0 );

   SCIPfreeBufferArray(scip, &componentobj);
   SCIPfreeBufferArray(scip, &componentends);

   return SCIP_OKAY;
}


/** For a symresack, determine a maximizer for optimizing linear function
 *  over a symresack, where the critical entry is fixed.
 */
static
SCIP_RETCODE maximizeObjectiveSymresackCriticalEntry(
   SCIP*                scip,                /**< SCIP pointer */
   int                  nvars,               /**< number of variables in symresack */
   SCIP_Real*           objective,           /**< the objective vector */
   int*                 perm,                /**< the permutation (without fixed points) as an array */
   int*                 invperm,             /**< the inverse permutation as an array */
   int                  crit,                /**< critical entry where optimality is found at */
   int*                 maxsolu              /**< pointer to the optimal objective array */
)
{
   /* Compute to which components all entries belong. */
   int* entrycomponent;
   SCIP_Real* componentobjective;

   int i;
   int c;

   assert( scip != NULL );
   assert( nvars > 0 );
   assert( objective != NULL );
   assert( perm != NULL );
   assert( invperm != NULL );
   assert( maxsolu != NULL );
   assert( crit >= 0 );
   assert( crit <= nvars );

   SCIP_CALL( SCIPallocBufferArray(scip, &entrycomponent, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &componentobjective, nvars) );

   /* Initially: Everything forms its own component */
   for (i = 0; i < nvars; ++i)
   {
      entrycomponent[i] = i;
      componentobjective[i] = objective[i];
   }
   for (i = 0; i < crit; ++i)
   {
      /* The graph with arcs {i, invperm[i]} if i < c is a collection of paths, cycles and singletons.
       * Label the vertices to the lowest entry in the component,  and store the value of that in this component.
       * Every inner while-loop labels one new vertex per iteration, and a vertex is relabeled exactly once.
       */

      /* Ignore fixed points. */
      if ( i == invperm[i] )
         continue;

      if ( entrycomponent[i] < i )
      {
         /* This entry is already included in a component. */
         continue;
      }

      /* Follow the path forward: Take edges {c, invperm[c]} until c >= crit, or a cycle is found. */
      c = i;
      while( c < crit )
      {
         /* c < crit, so edge {c, invperm[c]} exists. Label invperm[c] as part of component of i */
         c = invperm[c];

         /* Stop if we find a cycle. */
         if ( entrycomponent[c] != c )
            break;

         entrycomponent[c] = i;
         componentobjective[i] += objective[c];
      }

      /* Follow the path backward: Take edges {c, perm[c]} until perm[c] >= crit, or a cycle is found. */
      c = perm[i];
      while( c < crit )
      {
         /* c < crit, so edge {c, invperm[c]} exists. Label c as part of component of i */

         /* Stop if we find a cycle. */
         if ( entrycomponent[c] != c )
            break;

         entrycomponent[c] = i;
         componentobjective[i] += objective[c];
         /* For next iteration: We do another step back */
         c = perm[c];
      }
   }

   /* Now fill the objective vector.
    * If it is a fixed point, it is never part of the constraint.
    * For the component containing crit, set the value to 1.
    * For the component contraining invperm[crit], set the value to 0.
    * For the other components, set the value to 1 if the objective sum is positive.
    * Otherwise to 0.
    */
   for (i = 0; i < nvars; ++i)
   {
      if ( i == invperm[i] )
         maxsolu[i] = 0;
      else if ( entrycomponent[i] == entrycomponent[crit] )
         maxsolu[i] = 1;
      else if ( entrycomponent[i] == entrycomponent[invperm[crit]] )
         maxsolu[i] = 0;
      else if ( SCIPisGT(scip, componentobjective[entrycomponent[i]], 0.0) )
         maxsolu[i] = 1;
      else
         maxsolu[i] = 0;
   }

   SCIPfreeBufferArray(scip, &componentobjective);
   SCIPfreeBufferArray(scip, &entrycomponent);

   return SCIP_OKAY;
}

/** separate symresack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 */
static
SCIP_RETCODE separateSymresackCoversSymretope(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_Real*            vals,               /**< solution values of variables */
   int*                  ngen,               /**< pointer to store the number of separated covers */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real constobjective;
   SCIP_Real* sepaobjective;
   SCIP_Real maxsoluobj;
   int* maxsolu;
   int* invperm;
   int* perm;
   int nvars;
   int maxcrit;
   int i;
   int k;

   *infeasible = FALSE;
   *ngen = 0;

   assert( scip != NULL );
   assert( consdata != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->permutation != NULL );
   assert( consdata->permutation->perm != NULL );
   assert( consdata->permutation->order > 0 );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   nvars = consdata->nvars;

   /* allocate memory for temporary and global solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaobjective, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsolu, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &invperm, nvars) );

   for (k=1; k <= consdata->nperms; ++k)
   {
      SCIP_CALL( getPermArray(consdata->permutation, k, perm, nvars) );
      SCIP_CALL( getPermArray(consdata->permutation, -k, invperm, nvars) );

      #ifndef NDEBUG
      {
      /* Check: Is the permutation correct? */
      int dbg;
      for (dbg = 0; dbg < nvars; ++dbg)
      {
         assert( perm[invperm[dbg]] == dbg );
      }
      }
      #endif

      /* initialize objective */
      constobjective = 1.0; /* constant part of separation objective */
      for (i = 0; i < nvars; ++i)
      {
         if ( i < perm[i] )
            sepaobjective[i] = - vals[i];
         else if ( i > perm[i] )
         {
            sepaobjective[i] = 1.0 - vals[i];
            constobjective += vals[i] - 1.0;
         }
         else
            sepaobjective[i] = 0;
      }

      /* Find critical row of a maximally violated cover */
      SCIP_CALL( maximizeObjectiveSymresackStrict(scip, nvars, sepaobjective, perm, invperm, &maxcrit, &maxsoluobj) );
      assert( maxcrit >= 0 );
      assert( invperm[maxcrit] != maxcrit );
      SCIPdebugMsg(scip, "Critical row %d found; Computing maximally violated cover.\n", maxcrit);
      SCIP_CALL( maximizeObjectiveSymresackCriticalEntry(scip, nvars, sepaobjective, perm, invperm, maxcrit, maxsolu) );

      /* Add constant to maxsoluobj to get the real objective */
      maxsoluobj += constobjective;

      /* Check whether the separation objective is positive, i.e., a violated cover was found. */
      if ( SCIPisEfficacious(scip, maxsoluobj) )
      {
         /* Now add the cut. Reuse array maxsolu as coefficient vector for the constraint. */
         SCIP_Real rhs = -1.0;
         for (i = 0; i < nvars; ++i)
         {
            if ( i < perm[i] )
               maxsolu[i] = -maxsolu[i];
            else if ( i > perm[i] )
            {
               if ( maxsolu[i] == 0 )
                  rhs += 1.0;
               maxsolu[i] = 1 - maxsolu[i];
            }
            else
               maxsolu[i] = 0;
         }

         /* add cover inequality */
         SCIP_CALL( addSymresackInequality(scip, cons, nvars, consdata->vars, maxsolu, rhs, infeasible) );

         if ( ! *infeasible )
            ++(*ngen);
         if ( *infeasible )
            break;
      }
   }

   SCIPfreeBufferArrayNull(scip, &invperm);
   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArrayNull(scip, &maxsolu);
   SCIPfreeBufferArrayNull(scip, &sepaobjective);

   return SCIP_OKAY;
}


/** check whether solution is feasible for symresacks */
static
SCIP_RETCODE checkSymretopeSolution(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constrained for which we check the solution */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_RESULT*          result,             /**< pointer to store whether we detected infeasibility */
   SCIP_Bool             printreason         /**< whether reason for infeasibility should be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_PERMUTATION* permutation;
   SCIP_VAR** vars;
   int nvars;
   int i;
   int j;
   int k;
   SCIP_Real vali;
   SCIP_Real valj;
   int intvali;
   int intvalj;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->permutation != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->permutation != NULL );

   SCIPdebugMsg(scip, "Check method for symretope constraint <%s> (%d rows, %lld perms) ...\n", SCIPconsGetName(cons), consdata->nvars, consdata->permutation->order);

   vars = consdata->vars;
   nvars = consdata->nvars;
   permutation = consdata->permutation;

   /* detect first non-constant pair of variables */
   for (k = 1; k <= consdata->nperms; ++k)
   {
      for (i = 0; i < nvars; ++i)
      {
         j = permGet(permutation, i, -k);
         vali = SCIPgetSolVal(scip, sol, vars[i]);
         assert( SCIPisFeasIntegral(scip, vali) );
         intvali = vali > 0.5 ? 1 : 0;

         valj = SCIPgetSolVal(scip, sol, vars[j]);
         assert( SCIPisFeasIntegral(scip, valj) );
         intvalj = valj > 0.5 ? 1 : 0;

         if ( intvali < intvalj )
         {
            *result = SCIP_INFEASIBLE;
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            if ( printreason )
               SCIPinfoMessage(scip, NULL,
                  "Permutation perm[%d] has first non-constant pair (%d, %d) of variables has pattern (0,1).\n",
                  k, i, j);
            return SCIP_OKAY;
         }
         if ( intvali > intvalj )
            break;
         assert( intvali == intvalj );
      }
   }
   return SCIP_OKAY;
}


/** Upgrade symretope constraints to orbisacks */
static
SCIP_RETCODE orbisackUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            inputvars,          /**< permuted variables array */
   int                   nvars,              /**< size of perm array */
   SCIP_Bool*            upgrade,            /**< whether constraint was upgraded */
   SCIP_Bool             ismodelcons,        /**< whether the symretope is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows = 0;
   int i;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( nvars > 0 );
   assert( inputvars != NULL );
   assert( upgrade != NULL );

   *upgrade = TRUE;

   /* check whether orbisack conshdlr is available */
   conshdlr = SCIPfindConshdlr(scip, "orbisack");
   if ( conshdlr == NULL )
   {
      *upgrade = FALSE;
      SCIPdebugMsg(scip, "Cannot check whether symretope constraint can be upgraded to orbisack constraint. ");
      SCIPdebugMsg(scip, "---> Orbisack constraint handler not found.\n");

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, nvars) );

   /* check whether permutation is a composition of 2-cycles */
   for (i = 0; i < nvars; ++i)
   {
      /* ignore non-binary variables */
      if ( ! SCIPvarIsBinary(inputvars[i]) )
         continue;

      if ( perm[perm[i]] != i )
      {
         *upgrade = FALSE;
         break;
      }

      if ( perm[i] > i )
      {
         vars1[nrows] = inputvars[i];
         vars2[nrows++] = inputvars[perm[i]];

         assert( nrows <= nvars );
      }
   }

   /* if permutation can be upgraded to an orbisack */
   if ( nrows == 0 )
      *upgrade = FALSE;
   else if ( *upgrade )
   {
      SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nrows, FALSE, FALSE, ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


// /** Upgrade symretope constraints to symresack */
// static
// SCIP_RETCODE symresackUpgrade(
//    SCIP*                 scip,               /**< SCIP pointer */
//    SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
//    const char*           name,               /**< name of constraint */
//    int*                  perm,               /**< permutation */
//    SCIP_VAR**            vars,               /**< permuted variables array */
//    int                   nvars,              /**< size of perm array */
//    SCIP_Bool*            upgrade,            /**< whether constraint was upgraded */
//    SCIP_Bool             ismodelcons,        /**< whether the symretope is a model constraint */
//    SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
//                                               *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
//    SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
//                                               *   Usually set to TRUE. */
//    SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
//                                               *   TRUE for model constraints, FALSE for additional, redundant constraints. */
//    SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
//                                               *   TRUE for model constraints, FALSE for additional, redundant constraints. */
//    SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
//                                               *   Usually set to TRUE. */
//    SCIP_Bool             local,              /**< is constraint only valid locally?
//                                               *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
//    SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
//                                               *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
//                                               *   adds coefficients to this constraint. */
//    SCIP_Bool             dynamic,            /**< is constraint subject to aging?
//                                               *   Usually set to FALSE. Set to TRUE for own cuts which
//                                               *   are separated as constraints. */
//    SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
//                                               *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
//    SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
//                                               *   if it may be moved to a more global node?
//                                               *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
//    )
// {
//    SCIP_CONSHDLR* conshdlr;
//    SCIP_CONSHDLR* conshdlrsymretope;
//    SCIP_CONSHDLRDATA* conshdlrsymretopedata;
//    SCIP_Bool* checked;
//    int order;
//    int i;
//    int j;
//    int thiscyclesize;

//    assert( scip != NULL );
//    assert( perm != NULL );
//    assert( nvars > 0 );
//    assert( vars != NULL );
//    assert( upgrade != NULL );

//    *upgrade = FALSE;

//    /* check whether orbisack conshdlr is available */
//    conshdlr = SCIPfindConshdlr(scip, "symresack");
//    if ( conshdlr == NULL )
//    {
//       *upgrade = FALSE;
//       SCIPdebugMsg(scip, "Cannot check whether symretope constraint can be upgraded to symresack constraint. ");
//       SCIPdebugMsg(scip, "---> Symresack constraint handler not found.\n");

//       return SCIP_OKAY;
//    }

//    /* find the symretope constraint handler */
//    conshdlrsymretope = SCIPfindConshdlr(scip, CONSHDLR_NAME);
//    if ( conshdlr == NULL )
//    {
//       SCIPerrorMessage("Symretope constraint handler not found.\n");
//       return SCIP_PLUGINNOTFOUND;
//    }

//    /* We want to upgrade to symresack if the group order is too high. Compute the group order. */
//    order = 1;
//    SCIP_CALL( SCIPallocCleanBufferArray(scip, &checked, nvars) );
//    for (i = 0; i < nvars; ++i)
//    {
//       /* If this index is already processed, don't process. */
//       if ( checked[i] )
//          continue;

//       j = i;
//       thiscyclesize = 0;
//       do
//       {
//          checked[j] = TRUE;
//          j = perm[j];
//          ++thiscyclesize;
//       } while (j != i);
//       assert( j == i );

//       order = lcm(order, thiscyclesize);
//    }

//    for (i = 0; i < nvars; ++i)
//    {
//       assert( checked[i] );
//       checked[i] = FALSE;
//    }
//    SCIPfreeCleanBufferArray(scip, &checked);

//    /* If the order is sufficiently small, stop here and continue with symretope. */
//    conshdlrsymretopedata = SCIPconshdlrGetData(conshdlrsymretope);
//    assert( conshdlrsymretopedata != NULL );
//    if ( order <= conshdlrsymretopedata->maxorder )
//       return SCIP_OKAY;

//    /* Upgrade to symresack. */
//    *upgrade = TRUE;
//    SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars, ismodelcons,
//          initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
//    return SCIP_OKAY;
// }


/** creates a symmetry breaking constraint
 *
 * Depending on the given permutation, either an orbisack or symretope constraint
 * is created.
 */
SCIP_RETCODE SCIPcreateSymbreakConsSymretope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons,        /**< whether the added constraint is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_Bool upgrade = FALSE;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   /* Is it an orbisack? If so, make it an orbisack */
   SCIP_CALL( orbisackUpgrade(scip, cons, name, perm, vars, nvars, &upgrade, ismodelcons,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   if ( upgrade )
      return SCIP_OKAY;

   // /* Is the group size too large to be a symretope constraint? If so, make it a symresack. */
   // SCIP_CALL( symresackUpgrade(scip, cons, name, perm, vars, nvars, &upgrade, ismodelcons,
   //       initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   // if ( upgrade )
   //    return SCIP_OKAY;

   /* The group size is sufficiently small for a symretope constraint. Make it a symretope constraint. */
   SCIP_CALL( SCIPcreateConsSymretope(scip, cons, name, perm, vars, nvars, ismodelcons,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySymretope)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSymretope(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymretope)
{   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata, conshdlr) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeSymretope)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymretope)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* consdata = NULL;
   int nvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* constraint might be empty and not deleted if no presolving took place */
   assert( sourcedata->nvars == 0 || sourcedata->vars != NULL );
   assert( sourcedata->nvars == 0 || sourcedata->permutation != NULL );

   /* create transformed constraint data */
   nvars = sourcedata->nvars;

   /* Create consdata */
   SCIP_CALL( consdataCreate(scip, conshdlr, &consdata, sourcedata->vars, nvars,
      sourcedata->permutation->perm, sourcedata->ismodelcons) );

   /* create transformed constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSymretope)
{
   int c;

   assert( infeasible != NULL );
   *infeasible = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );

      SCIPdebugMsg(scip, "Generating initial symresack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;
   }
   SCIPdebugMsg(scip, "Generated initial symresack cuts.\n");

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSymretope)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* determine maximum number of vars in a symresack constraint */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   conshdlrdata->maxnvars = 0;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* update conshdlrdata if necessary */
      if ( consdata->nvars > conshdlrdata->maxnvars )
         conshdlrdata->maxnvars = consdata->nvars;
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSymretope)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int maxnvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if ( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   maxnvars = conshdlrdata->maxnvars;
   assert( maxnvars > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

   *result = SCIP_DIDNOTFIND;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symretope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      if ( consdata->nvars == 0 )
         continue;

      /* get solution */
      assert( consdata->nvars <= maxnvars );
      SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCoversSymretope(scip, conss[c], consdata, vals, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         goto FREESEPALP;
      }

      if ( ngen > 0 )
      {
         *result = SCIP_SEPARATED;
         if ( !conshdlrdata->sepaallviolperms )
            goto FREESEPALP;
      }
   }

   FREESEPALP:
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolSymretope)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int maxnvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* COMMENT: This function is essentially a copy of CONSSEPALP. I suggest to introduce another function
    * that we call in both CONSSEPALP and CONSEPASOL.
    */

   if ( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   maxnvars = conshdlrdata->maxnvars;
   assert( maxnvars > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

   *result = SCIP_DIDNOTFIND;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      if ( consdata->nvars == 0 )
         continue;

      /* get solution */
      assert( consdata->nvars <= maxnvars );
      SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCoversSymretope(scip, conss[c], consdata, vals, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         goto FREESEPASOL;
      }

      if ( ngen > 0 )
      {
         *result = SCIP_SEPARATED;
         if ( !conshdlrdata->sepaallviolperms )
            goto FREESEPASOL;
      }
   }

 FREESEPASOL:
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions.
 *
 *  To check feasibility, we separate cover inequalities.
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (lp solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   /* COMMENT: also use function to not copy the same functionality in CONSENFOLP and CONSENFORELAX */

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real* vals;
      int maxnvars;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      maxnvars = conshdlrdata->maxnvars;
      assert( maxnvars > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symretope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         /* do not enforce non-model constraints */
         if ( !consdata->ismodelcons )
            continue;

         if ( consdata->nvars == 0 )
            continue;

         /* get solution */
         assert( consdata->nvars <= maxnvars );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCoversSymretope(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         /* SCIPdebugMsg(scip, "Generated symresack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen); */

         if ( ngen > 0 )
         {
            *result = SCIP_SEPARATED;
            if ( !conshdlrdata->sepaallviolperms )
            {
               SCIPfreeBufferArray(scip, &vals);
               return SCIP_OKAY;
            }
         }
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkSymretopeSolution(scip, conss[c], NULL, result, FALSE) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions
 *
 *  To check feasibility, we separate cover inequalities.
 *
 */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSymretope)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (relaxation solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real* vals;
      int maxnvars;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      maxnvars = conshdlrdata->maxnvars;
      assert( maxnvars > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         /* do not enforce non-model constraints */
         if ( !consdata->ismodelcons )
            continue;

         if ( consdata->nvars == 0 )
            continue;

          /* get solution */
         assert( consdata->nvars <= maxnvars );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCoversSymretope(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         if ( ngen > 0 )
         {
            *result = SCIP_SEPARATED;
            if ( !conshdlrdata->sepaallviolperms )
            {
               SCIPfreeBufferArray(scip, &vals);
               return SCIP_OKAY;
            }
         }
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymretope)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* do not check non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkSymretopeSolution(scip, conss[c], sol, result, printreason) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   if ( *result == SCIP_FEASIBLE )
      SCIPdebugMsg(scip, "Solution is feasible.\n");

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   SCIP_Bool success = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Propagation method of symretope constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* Only propagate if there is reason to. */
      if ( !consdata->execprop )
         continue;

      /* Clear the affected entries list */
      for (i = 0; i < consdata->nvars; ++i)
         consdata->affectedentries[i] = FALSE;

      SCIP_CALL( propVariables(scip, conss[c], NULL, TRUE, consdata->affectedentries, &infeasible, &ngen) );

      /* If this subtree is infeasible, cutoff. */
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* If it's not infeasible, do not propagate again until an affected variable is changed. */
      consdata->execprop = FALSE;

      success = success || ( ngen > 0 );

      *result = SCIP_DIDNOTFIND;
   }

   if ( success )
   {
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymretope)
{  /*lint --e{715}*/
   int c;
   SCIP_Bool success = FALSE;
   int oldndelconss;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   oldndelconss = *ndelconss;

   SCIPdebugMsg(scip, "Presolving method of symretope constraint handler. Propagating symretope constraints.\n");
   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_CONSDATA* consdata;
      int ngen = 0;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* avoid trivial problems */
      if ( consdata->nvars == 0 )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         (*ndelconss)++;
      }
      else
      {
         SCIP_CALL( propVariables(scip, conss[c], NULL, TRUE, NULL, &infeasible, &ngen) );
      }

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( ngen > 0 )
      {
         *nfixedvars += ngen;
         success = TRUE;
      }

      *result = SCIP_DIDNOTFIND;
   }

   if ( *ndelconss > oldndelconss ||  success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropSymretope)
{  /*lint --e{715}*/
   int nvars;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_PERMUTATION* permutation;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution method of symretope constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   /* Get instance data */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nvars >= 0 );
   assert( consdata->permutation != NULL );
   vars = consdata->vars;
   nvars = consdata->nvars;
   permutation = consdata->permutation;

   if ( inferinfo < 0 )
   {
      /* If inferinfo < 0, then this is due to peeking.
       * First re-run the propagation algorithm with the converse fixing. This will find infeasibility.
       * While running, store the unfixed variables for which we had to look up the value.
       * Then sparsify: Do we also find infeasibility if some of these are considered UNFIXED?
       */
      SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixings;
      SCIP_SYMRETOPEVIRTUALFIXINGS* virtualfixingsinitial;
      SCIP_Bool* conflictentries;
      SCIP_Bool infeasible;
      int ngen;
      int i;
      int j;
      int infervarid;

      /* Determine the variable index of the fixed variable */
      for (infervarid = 0; infervarid < nvars; ++infervarid)
      {
         if ( vars[infervarid] == infervar )
            break;
      }
      assert( infervarid < nvars );

      /* Initialize virtual fixings.
       * We do not use global bound information, so we store the actual fixing in virtual fixings.
       */
      allocVirtualFixings(scip, &virtualfixings, nvars);
      allocVirtualFixings(scip, &virtualfixingsinitial, nvars);
      for (j = 0; j < nvars; ++j)
      {
         if ( SCIPvarGetUbAtIndex(vars[j], bdchgidx, FALSE) < 0.5 )
            setVirtualFixing(virtualfixings, j, FIXED0);
         else if ( SCIPvarGetLbAtIndex(vars[j], bdchgidx, FALSE) > 0.5 )
            setVirtualFixing(virtualfixings, j, FIXED1);
      }

      /* If the lower bound gets set, then infervar is fixed to 1. (Likewise, if the upper bound is set, fixed to 0)
       * In the following, the goal is to certify infeasibility for the converse fixing.
       * Set converse fixing as virtual fixing.
       */
      assert( getVirtualFixing(virtualfixings, infervarid) == UNFIXED );
      setVirtualFixing(virtualfixings, infervarid, boundtype == SCIP_BOUNDTYPE_LOWER ? FIXED0 : FIXED1);

      /* Copy virtualfixings, as propVariables adds newly found fixings to virtualfixings. */
      copyVirtualFixings(virtualfixings, virtualfixingsinitial);

      /* Store which entries are being looked up. */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &conflictentries, nvars) );

      /* Run the propagation algorithm with the converse fixing applied. This will definitively yield infeasible.
       * In this, do not use global bound information, only the fixings specified in virtualfixings.
       */
      SCIP_CALL( propVariables(scip, cons, virtualfixings, FALSE, conflictentries, &infeasible, &ngen) );
      assert( infeasible );

      /* Now the conflict consists of checkedentries for sure.
       * Make the conflict sparser, by removing virtual fixings from the input as long as it finds feasibility.
       */
      for (i = 0; i < nvars; ++i)
      {
         /* We must certify that infervar gets fixed, so cannot be part of the conflict */
         if ( i == infervarid )
         {
            conflictentries[i] = FALSE;
            continue;
         }

         /* Unfixed entries are not of interest for the conflict. */
         if ( getVirtualFixing(virtualfixingsinitial, i) == UNFIXED )
         {
            conflictentries[i] = FALSE;
            continue;
         }

         if ( conflictentries[i] )
         {
            /* Check what happens if i is not part of the conflict. */
            int entry;

            clearVirtualFixings(virtualfixings);
            for (j = 0; j < virtualfixingsinitial->nvirtualfixings; ++j)
            {
               entry = virtualfixingsinitial->entrystack[j];
               if ( entry == infervarid || (entry != i && conflictentries[entry]) )
                  setVirtualFixing(virtualfixings, entry, virtualfixingsinitial->entrylookup[entry]);
            }

            /* Run propagator virtually to determine if it is feasible or not. */
            SCIP_CALL( propVariables(scip, cons, virtualfixings, FALSE, NULL, &infeasible, &ngen) );

            /* If entry i is not necessary for certifying infeasibility, remove from conflict. */
            if ( infeasible )
               conflictentries[i] = FALSE;
         }
      }

      /* add the conflict */
      for (j = 0; j < nvars; ++j)
      {
         /* We must certify that infervar gets fixed, so cannot be part of the conflict */
         if ( j == infervarid )
            continue;

         /* Add to the conflict if it is looked up by the code, and fixed to 0 or 1. */
         if ( conflictentries[j] )
         {
            if ( SCIPvarGetUbAtIndex(vars[j], bdchgidx, FALSE) < 0.5 )
            {
               SCIP_CALL( SCIPaddConflictUb(scip, vars[j], bdchgidx) );
            }
            else if ( SCIPvarGetLbAtIndex(vars[j], bdchgidx, FALSE) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictLb(scip, vars[j], bdchgidx) );
            }
         }
      }

      SCIPfreeBufferArray(scip, &conflictentries);
      freeVirtualFixings(scip, &virtualfixingsinitial);
      freeVirtualFixings(scip, &virtualfixings);

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
   else
   {
      /* If inferinfo >= 0, then this is found without peeking. Inferinfo encodes the permutation at which it fails. */
      /* Get the power of the permutation for which this fixing is found. This is stored in "inferinfo" */
      /* Now list the conflict that makes sure that the bound "boundtype" must be tightened for variable "infervar" */
      SCIP_CALL( resolveSymretopeConflictVariables(scip, infervar, boundtype, vars, nvars, permutation, inferinfo, bdchgidx) );

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
}


/** lock variables
 *
 *  We assume we have only one global (void) constraint and lock all binary variables
 *  which do not correspond to fixed points of the permutation.
 *
 * - Symresack constraints may get violated if the variables with a negative coefficient
 *   in the FD inequality are rounded down, we therefor call
 *   SCIPaddVarLocksType(..., nlockspos, nlocksneg).
 * - Symresack constraints may get violated if the variables with a positive coefficient
 *   in the FD inequality are rounded up, we therefor call
 *   SCIPaddVarLocksType(..., nlocksneg, nlockspo ).
 * - Symretope constraints can be violated if a variable is rounded up or down.
 */
static
SCIP_DECL_CONSLOCK(consLockSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_PERMUTATION* permutation;
   int c;
   int i;
   int* cycle;
   int cyclen;
   int cycmax;
   int cycmin;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking method for symretope constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );

   vars = consdata->vars;

   assert( consdata->permutation != NULL );
   permutation = consdata->permutation;
   assert( permutation->cycles != NULL );
   assert( permutation->cyclelengths != NULL );

   for (c = 0; c < permutation->ncycles; ++c)
   {
      cycle = permutation->cycles[c];
      cyclen = permutation->cyclelengths[c];
      assert( cyclen > 0 );

      /* If the cycle is trivial, then the value does not matter at all for this constraint. */
      if ( cyclen == 1 )
         continue;

      /* Determine minimal and maximal entries of the cycles */
      cycmin = cycle[0];
      cycmax = cycle[0];
      for (i = 1; i < cyclen; ++i)
      {
         assert( cycmax >= cycmin );
         if ( cycle[i] > cycmax )
            cycmax = cycle[i];
         else if ( cycle[i] < cycmin )
            cycmin = cycle[i];
      }

      for (i = 0; i < cyclen; ++i)
      {
         /* Minimal entry in orbit, then rounding down could cause infeasibilities. */
         if ( cycle[i] == cycmin )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, vars[cycle[i]], locktype, nlockspos, nlocksneg) );
         }
         /* Maximal entry in orbit, then rounding up could cause infeasibilities. */
         else if ( cycle[i] == cycmax )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, vars[cycle[i]], locktype, nlocksneg, nlockspos) );
         }
         /* Other entries in orbit, then rounding either way could cause infeasibilities. */
         else
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, vars[cycle[i]], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySymrestope)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for symresack constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->permutation != NULL );
   assert( sourcedata->permutation->perm != NULL );
   assert( sourcedata->nvars > 0 );

   conshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert( conshdlrdata != NULL );

   /* do not copy non-model constraints */
   if ( !sourcedata->ismodelcons && !conshdlrdata->forceconscopy )
   {
      *valid = FALSE;

      return SCIP_OKAY;
   }

   sourcevars = sourcedata->vars;
   nvars = sourcedata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for (i = 0; i < nvars && *valid; ++i)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i], &(vars[i]), varmap, consmap, global, valid) );
      assert( !(*valid) || vars[i] != NULL );
   }

   /* only create the target constraint, if all variables could be copied */
   if ( *valid )
   {
      /* create copied constraint */
      if ( name == NULL )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsSymretope(scip, cons, name, sourcedata->permutation->perm, vars, nvars,
            sourcedata->ismodelcons, initial, separate, enforce, check, propagate, local,
            modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSymretope)
{  /*lint --e{715}*/
   const char* s;
   char* endptr;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int* perm;
   int val;
   int nvars = 0;
   int cnt = 0;
   int nfoundpermidx = 0;
   int maxnvars = 128;

   assert( success != NULL );

   *success = TRUE;
   s = str;

   /* skip white space */
   while ( *s != '\0' && isspace((unsigned char)*s) )
      ++s;

   if ( strncmp(s, "symretope(", 10) != 0 )
   {
      SCIPerrorMessage("Syntax error - expected \"symretope(\", but got '%s'", s);
      *success = FALSE;
      return SCIP_OKAY;
   }
   s += 10;

   /* loop through string */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, maxnvars) );

   do
   {
      if ( cnt > 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected two arrays, but got more\n");
         *success = FALSE;

         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &vars);
      }

      /* skip whitespace and ',' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
         ++s;

      /* if we could not find starting indicator of array */
      if ( *s != '[' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '[' to start new array\n");
         *success = FALSE;

         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &vars);
      }
      ++s;

      /* read array, cnt = 0: variables; cnt = 1: permutation*/
      if ( cnt == 0 )
      {
         do
         {
            /* skip whitespace and ',' */
            while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
               ++s;

            /* parse variable name */
            SCIP_CALL( SCIPparseVarName(scip, s, &var, &endptr) );
            if ( var == NULL )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
               *success = FALSE;
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
            s = endptr;
            assert( s != NULL );

            vars[nvars++] = var;

            if ( nvars >= maxnvars )
            {
               int newsize;

               newsize = SCIPcalcMemGrowSize(scip, nvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &vars, newsize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &perm, newsize) );
               maxnvars = newsize;
            }
         }
         while ( *s != ']' );
      }
      else
      {
         do
         {
            /* skip whitespace and ',' */
            while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' ) )
               ++s;

            /* parse integer value */
            if ( ! SCIPstrToIntValue(s, &val, &endptr) )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "could not extract int from string '%s'\n", str);
               *success = FALSE;
               SCIPfreeBufferArray(scip, &perm);
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
            s = endptr;
            assert( s != NULL );

            perm[nfoundpermidx++] = val;

            if ( nfoundpermidx > nvars )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "permutation is longer than vars array\n");
               *success = FALSE;
               SCIPfreeBufferArray(scip, &perm);
               SCIPfreeBufferArray(scip, &vars);
               return SCIP_OKAY;
            }
         }
         while ( *s != ']' );
      }
      ++s;
      ++cnt;
   }
   while ( *s != ')' );

   if ( nfoundpermidx == nvars )
   {
      /* Do NOT add symretope as a model constraint. */
      SCIP_CALL( SCIPcreateConsBasicSymretope(scip, cons, name, perm, vars, nvars, FALSE) );
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
         "Length of permutation is not equal to number of given variables.\n");
      *success = FALSE;
   }

   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "Printing method for symretope constraint handler\n");

   /* we do not have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->permutation != NULL );
   assert( consdata->permutation->perm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   perm = consdata->permutation->perm;

   SCIPinfoMessage(scip, file, "symretope([");
   SCIP_CALL( SCIPwriteVarName(scip, file, vars[0], TRUE) );

   for (i = 1; i < nvars; ++i)
   {
      SCIPinfoMessage(scip, file, ",");
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[i], TRUE) );
   }
   SCIPinfoMessage(scip, file, "],[%d", perm[0]);
   for (i = 1; i < nvars; ++i)
      SCIPinfoMessage(scip, file, ",%d", perm[i]);
   SCIPinfoMessage(scip, file, "])");

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;

      for (i = 0; i < consdata->nvars; ++i)
         vars[cnt++] = consdata->vars[i];
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSymretope)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( nvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/** creates the handler for symresack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSymretope(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSymretope, consEnfopsSymretope, consCheckSymretope, consLockSymretope,
         conshdlrdata) );
   assert( conshdlr != NULL );

   /* include event handler */
   conshdlrdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr,
      EVENTHDLR_SYMRETOPE_NAME, EVENTHDLR_SYMRETOPE_DESC, eventExec, NULL) );
   assert( conshdlrdata->eventhdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySymretope, consCopySymrestope) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSymretope) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymretope) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymretope) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSymretope) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSymretope) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSymretope) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymretope, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSymretope) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSymretope, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSymretope) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSymretope, consSepasolSymretope, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymretope) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSymretope) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSymretope) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forceconscopy",
         "Whether symresack constraints should be forced to be copied to sub SCIPs.",
         &conshdlrdata->forceconscopy, TRUE, DEFAULT_FORCECONSCOPY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/peek",
         "Whether additional symretope fixings should be determined by testing feasibility by testing unfixed entries.",
         &conshdlrdata->symretopepeek, TRUE, DEFAULT_SYMRETOPEPEEK, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxgrouporder",
         "Maximal group order for symretope constraint before restricting the number of considered permutations.",
         &conshdlrdata->maxorder, TRUE, DEFAULT_SYMRETOPEMAXORDER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxgroupordernvars",
         "Maximal value of group order multiplied with group support  before restricting number of permutations.",
         &conshdlrdata->maxordernvars, TRUE, DEFAULT_SYMRETOPEMAXORDERNVARS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/sepaallviolperms",
         "Whether a separating inequality should be added only for one violated symresack (FALSE) or for all violating symresacks (TRUE)",
         &conshdlrdata->sepaallviolperms, TRUE, DEFAULT_SEPAALLVIOLPERMS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/probingpeek",
         "Whether peeking should be done during probing.",
         &conshdlrdata->probingpeek, TRUE, DEFAULT_PROBINGPEEK, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a symretope constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSymretope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation generating the symretope group */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons,        /**< whether the symretope is a model constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( nvars > 0 );

   /* find the symretope constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Symretope constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, conshdlr, &consdata, vars, nvars, perm, ismodelcons) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures a symretope constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values
 *
 *  In a presolving step, we remove all fixed points and cycles that act on non-binary variables of the permutation
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsBasicSymretope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation generating the symretope group */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
   SCIP_Bool             ismodelcons         /**< whether the symretope is a model constraint */
   )
{
   SCIP_CALL( SCIPcreateConsSymretope(scip, cons, name, perm, vars, nvars, ismodelcons,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
