/* Copyright 2004,2007,2008,2010,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : parser.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a static mapper.                **/
/**                This module is the strategy lexical and **/
/**                syntactic analyzer.                     **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     02 may 1996     **/
/**                # Version 3.2  : from : 07 oct 1996     **/
/**                                 to     19 oct 1996     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     10 sep 2001     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to     02 feb 2004     **/
/**                # Version 5.0  : from : 20 feb 2008     **/
/**                                 to     20 feb 2008     **/
/**                # Version 5.1  : from : 22 oct 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 01 jun 2012     **/
/**                                 to     30 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define PARSER

#include "module.h"
#include "common.h"

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

#include "parser.h"
#include "parser_yy.h"

/*
**  The static and global variables.
*/

static StratTab             stratdummytab = { NULL, NULL, NULL }; /* Dummy strategy table for the dummy empty object       */
Strat                       stratdummy = { &stratdummytab, STRATNODEEMPTY }; /* Dummy empty object for offset computations */

/**************************/
/*                        */
/* The strategy routines. */
/*                        */
/**************************/

/* This routine parses the given strategy
** string and builds the corresponding
** strategy tree.
** It returns:
** - !NULL  : pointer to the strategy.
** - NULL   : on error.
*/

Strat *
stratInit (
const StratTab * const      strattab,             /*+ Pointer to strategy parsing table +*/
const char * const          string)               /*+ Strategy string to parse          +*/
{
#ifdef SCOTCH_DEBUG_PARSER1
  if ((strattab == NULL) || (string == NULL)) {
    errorPrint ("stratInit: invalid parameter");
    return     (NULL);
  }
#endif /* SCOTCH_DEBUG_PARSER1 */

  return (stratParserParse (strattab, string));   /* Parse strategy string */
}

/* This routine frees a strategy structure.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
stratExit (
Strat * const               strat)
{
  StratParamTab *   paratab;                      /* Table of method parameters                  */
  byte *            paraofft;                     /* Offset of parameter within method structure */
  unsigned int      i;
  int               o;

  if (strat == NULL)                              /* If node does not exist, abort */
    return (0);

  o = 0;                                          /* Assume everything will be all right */
  switch (strat->type) {                          /* Recursively free sub-strategies     */
    case STRATNODECONCAT :
      o  = stratExit (strat->data.concat.strat[0]);
      o |= stratExit (strat->data.concat.strat[1]);
      break;
    case STRATNODECOND :
      o  = stratTestExit (strat->data.cond.test);
      o |= stratExit     (strat->data.cond.strat[0]);
      if (strat->data.cond.strat[1] != NULL)
        o |= stratExit (strat->data.cond.strat[1]);
      break;
    case STRATNODESELECT :
      o  = stratExit (strat->data.select.strat[0]);
      o |= stratExit (strat->data.select.strat[1]);
      break;
    case STRATNODEEMPTY :                         /* Empty strategy node         */
      if (strat == &stratdummy)                   /* If node is empty dummy node */
        return (0);                               /* Return without freeing it   */
      break;
    case STRATNODEMETHOD :                        /* Method strategy node       */
      paratab = strat->tabl->paratab;             /* Free the method parameters */
      for (i = 0; paratab[i].name != NULL; i ++) {
        if ((paratab[i].meth == strat->data.method.meth) && /* For all parameters of that method    */
            (paratab[i].type == STRATPARAMSTRAT)) { /* Which are non-deprecated strategy parameters */
          paraofft = (byte *) &strat->data.method.data + /* Compute parameter offset within method  */
                      (paratab[i].dataofft -
                       paratab[i].database);
          o |= stratExit (*((Strat **) paraofft)); /* Perform recursion */
        }
      }
      break;
#ifdef SCOTCH_DEBUG_PARSER2
    default :
      errorPrint ("stratExit: invalid strategy node");
      o = 1;
      break;
#endif /* SCOTCH_DEBUG_PARSER2 */
  }

  memFree (strat);                                /* Free strategy structure itself */
  return  (o);                                    /* Return output code             */
}

/* This routine displays the given
** strategy structure.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
stratSave (
const Strat * const         strat,
FILE * const                stream)
{
  unsigned int      paraflag;                     /* Flag set if method has parameters           */
  StratParamTab *   paratab;                      /* Pointer to method parameter table           */
  byte *            paraofft;                     /* Offset of parameter within method structure */
  unsigned int      i;
  int               o;

  o = 0;
  switch (strat->type) {                          /* Recursively view sub-strategies */
    case STRATNODECOND :
      if ((fprintf (stream, "(/(") == EOF) ||
          (stratTestSave (strat->data.cond.test, stream) != 0) ||
          (fprintf (stream, ")?(") == EOF) ||
          (stratSave (strat->data.cond.strat[0], stream) != 0))
        o = 1;
      if ((o == 0) && (strat->data.cond.strat[1] != NULL)) {
        if ((fprintf (stream, "):(") == EOF) ||
            (stratSave (strat->data.cond.strat[1], stream) != 0))
          o = 1;
      }
      if (o == 0)
        o = (fprintf (stream, ");)") == EOF);
      break;
    case STRATNODECONCAT :
      if ((stratSave (strat->data.concat.strat[0], stream) != 0) ||
          (stratSave (strat->data.concat.strat[1], stream) != 0))
        o = 1;
      break;
    case STRATNODESELECT :
      if ((fprintf (stream, "(") == EOF) ||
          (stratSave (strat->data.select.strat[0], stream) != 0) ||
          (fprintf (stream, "|") == EOF) ||
          (stratSave (strat->data.select.strat[1], stream) != 0) ||
          (fprintf (stream, ")") == EOF))
        o = 1;
    case STRATNODEEMPTY :
      break;
    case STRATNODEMETHOD :
      if (fprintf (stream, "%s", strat->tabl->methtab[strat->data.method.meth].name) == EOF) { /* Print method name */
        o = 1;
        break;
      }
      paraflag = 0;                               /* No method parameters seen yet */
      paratab  = strat->tabl->paratab;
      for (i = 0; paratab[i].name != NULL; i ++) {
        if ((paratab[i].meth == strat->data.method.meth) && /* For all parameters of that method  */
            ((paratab[i].type & STRATPARAMDEPRECATED) == 0)) { /* Which are not deprecated        */
          paraofft = (byte*) &strat->data.method.data + /* Compute parameter offset within method */
                     (paratab[i].dataofft -
                      paratab[i].database);
          if (fprintf (stream, "%c%s=",           /* Open or continue parameter list */
                       ((paraflag ++ == 0) ? '{' : ','),
                       paratab[i].name) == EOF) {
            o = 1;
            break;
          }
          switch (paratab[i].type) {              /* Print parameter value         */
            case STRATPARAMCASE :                 /* Case value                    */
              o = (fprintf (stream, "%c",         /* Print corresponding character */
                            ((char *) paratab[i].datasltr)[*((unsigned int *) paraofft)]) == EOF);
              break;
            case STRATPARAMINT :                  /* Integer value */
              o = (fprintf (stream, INTSTRING, *((INT *) paraofft)) == EOF);
              break;
            case STRATPARAMDOUBLE :               /* Double value */
              o = (fprintf (stream, "%g", *((double *) paraofft)) == EOF);
              break;
            case STRATPARAMSTRAT :                /* Strategy                      */
              o = stratSave (*((Strat **) paraofft), stream); /* Perform recursion */
              break;
            case STRATPARAMSTRING :               /* String value */
              o = (fprintf (stream, "%s", (char *) paraofft) == EOF);
              break;
#ifdef SCOTCH_DEBUG_PARSER2
            default :
              errorPrint ("stratSave: invalid parameter type");
              return     (1);
#endif /* SCOTCH_DEBUG_PARSER2 */
          }
        }
        if (o != 0)                               /* If an error has occured */
          break;                                  /* Abort the loop          */
      }
      if ((o == 0) && (paraflag != 0))            /* If there is a parameter list */
        o |= (fprintf (stream, "}") == EOF);      /* Close it                     */
      break;
#ifdef SCOTCH_DEBUG_PARSER2
    default :
      errorPrint ("stratSave: invalid strategy node");
      return     (1);
#endif /* SCOTCH_DEBUG_PARSER2 */
  }
  if (o != 0) {
    errorPrint ("stratSave: bad output");
  }

  return (o);
}

/*****************************************/
/*                                       */
/* These routines handle strategy tests. */
/*                                       */
/*****************************************/

/* This routine evaluates the
** given condition.
** It returns:
** - 0   : on success; eval updated.
** - !0  : on error.
*/

int
stratTestEval (
const StratTest * restrict const  test,
StratTest * restrict const        eval,           /*+ Place where to return final value                      +*/
const void * restrict const       data)           /*+ Pointer to data structure where to read variables from +*/
{
  StratTest         val[2];                       /* Temporary evaluation variables */
  StratTestType     sign;                         /* Sign of comparison             */
  int               o;

#ifdef SCOTCH_DEBUG_PARSER1
  if ((test == NULL) || (eval == NULL) || (data == NULL)) {
    errorPrint ("stratTestEval: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_PARSER1 */

  o = 0;                                          /* Assume no error */
  switch (test->typetest) {
    case STRATTESTNOT :                           /* Not operator */
      o = stratTestEval (test->data.test[0], eval, data);
#ifdef SCOTCH_DEBUG_PARSER2
      if ((o == 0) && (eval->typenode != STRATPARAMLOG)) {
        errorPrint ("stratTestEval: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_PARSER2 */
      eval->data.val.vallog = 1 - eval->data.val.vallog;
      break;
    case STRATTESTAND :                           /* And operator */
      o = stratTestEval (test->data.test[0], eval, data);
#ifdef SCOTCH_DEBUG_PARSER2
      if ((o == 0) && (eval->typenode != STRATPARAMLOG)) {
        errorPrint ("stratTestEval: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_PARSER2 */
      if ((o == 0) && (eval->data.val.vallog == 1)) {
        o = stratTestEval (test->data.test[1], eval, data);
#ifdef SCOTCH_DEBUG_PARSER2
        if (eval->typenode != STRATPARAMLOG) {
          errorPrint ("stratTestEval: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_PARSER2 */
      }
      break;
    case STRATTESTOR :                            /* Or operator */
      o = stratTestEval (test->data.test[0], eval, data);
#ifdef SCOTCH_DEBUG_PARSER2
      if ((o == 0) && (eval->typenode != STRATPARAMLOG)) {
        errorPrint ("stratTestEval: internal error (4)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_PARSER2 */
      if ((o == 0) && (eval->data.val.vallog == 0)) {
        o = stratTestEval (test->data.test[1], eval, data);
#ifdef SCOTCH_DEBUG_PARSER2
        if (eval->typenode != STRATPARAMLOG) {
          errorPrint ("stratTestEval: internal error (5)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_PARSER2 */
      }
      break;
    case STRATTESTLT :                            /* Less-than operator    */
    case STRATTESTEQ :                            /* Equal-to operator     */
    case STRATTESTGT :                            /* Greater-than operator */
      o  = stratTestEval (test->data.test[0], &val[0], data);
      o |= stratTestEval (test->data.test[1], &val[1], data);
      o |= stratTestEvalCast (&val[0], &val[1]);
      if (o != 0)
        break;
      sign = STRATTESTNBR;                        /* In case of error */
      switch (val[0].typenode) {
        case STRATPARAMDOUBLE :
          sign = (val[0].data.val.valdbl < val[1].data.val.valdbl) ? STRATTESTLT : ((val[0].data.val.valdbl > val[1].data.val.valdbl) ? STRATTESTGT : STRATTESTEQ);
          break;
        case STRATPARAMINT :
          sign = (val[0].data.val.valint < val[1].data.val.valint) ? STRATTESTLT : ((val[0].data.val.valint > val[1].data.val.valint) ? STRATTESTGT : STRATTESTEQ);
          break;
#ifdef SCOTCH_DEBUG_PARSER2
        default :
          errorPrint ("stratTestEval: internal error (6)");
          o = 1;
          break;
#endif /* SCOTCH_DEBUG_PARSER2 */
      }
      eval->typenode        = STRATPARAMLOG;      /* Build test result */
      eval->data.val.vallog = (sign == test->typetest);
      break;
    case STRATTESTADD :                           /* Addition operator */
      o  = stratTestEval (test->data.test[0], &val[0], data);
      o |= stratTestEval (test->data.test[1], &val[1], data);
      o |= stratTestEvalCast (&val[0], &val[1]);
      if (o != 0)
        break;
      if (val[0].typenode == STRATPARAMDOUBLE)
        eval->data.val.valdbl = val[0].data.val.valdbl + val[1].data.val.valdbl;
      else
        eval->data.val.valint = val[0].data.val.valint + val[1].data.val.valint;
      eval->typenode = val[0].typenode;
      break;
    case STRATTESTSUB :                           /* Subtraction operator */
      o  = stratTestEval (test->data.test[0], &val[0], data);
      o |= stratTestEval (test->data.test[1], &val[1], data);
      o |= stratTestEvalCast (&val[0], &val[1]);
      if (o != 0)
        break;
      if (val[0].typenode == STRATPARAMDOUBLE)
        eval->data.val.valdbl = val[0].data.val.valdbl - val[1].data.val.valdbl;
      else
        eval->data.val.valint = val[0].data.val.valint - val[1].data.val.valint;
      eval->typenode = val[0].typenode;
      break;
    case STRATTESTMUL :                           /* Multiplication operator */
      o  = stratTestEval (test->data.test[0], &val[0], data);
      o |= stratTestEval (test->data.test[1], &val[1], data);
      o |= stratTestEvalCast (&val[0], &val[1]);
      if (o != 0)
        break;
      if (val[0].typenode == STRATPARAMDOUBLE)
        eval->data.val.valdbl = val[0].data.val.valdbl * val[1].data.val.valdbl;
      else
        eval->data.val.valint = val[0].data.val.valint * val[1].data.val.valint;
      eval->typenode = val[0].typenode;
      break;
    case STRATTESTMOD :                           /* Modulus operator */
      o  = stratTestEval (test->data.test[0], &val[0], data);
      o |= stratTestEval (test->data.test[1], &val[1], data);
      o |= stratTestEvalCast (&val[0], &val[1]);
      if (o != 0)
        break;
      if (val[0].typenode == STRATPARAMDOUBLE)
        eval->data.val.valdbl = fmod (val[0].data.val.valdbl, val[1].data.val.valdbl);
      else
        eval->data.val.valint = val[0].data.val.valint % val[1].data.val.valint;
      eval->typenode = val[0].typenode;
      break;
    case STRATTESTVAL :                           /* Constant value */
      *eval = *test;                              /* Copy value     */
      break;
    case STRATTESTVAR :                           /* Variable */
      switch (test->typenode) {
        case STRATPARAMDOUBLE :
          eval->data.val.valdbl = *((double *) ((byte *) data + test->data.var.datadisp));
          break;
        case STRATPARAMINT :
          eval->data.val.valint = *((INT *) ((byte *) data + test->data.var.datadisp));
          break;
#ifdef SCOTCH_DEBUG_PARSER1
        default :
          errorPrint ("stratTestEval: internal error (7)");
          o = 1;
          break;
#endif /* SCOTCH_DEBUG_PARSER1 */
      }
      eval->typenode = test->typenode;
      break;
#ifdef SCOTCH_DEBUG_PARSER1
    default :
      errorPrint ("stratTestEval: invalid condition type (%u)", test->typetest);
      o = 1;
      break;
#endif /* SCOTCH_DEBUG_PARSER1 */
  }
  eval->typetest = STRATTESTVAL;

  return (o);
}

/* This routine casts the type of one
** of the two input values so as to
** get the same type for both values.
** It returns:
** - VOID  : in all cases;
*/

static
int
stratTestEvalCast (
StratTest * const           test0,
StratTest * const           test1)
{
#ifdef SCOTCH_DEBUG_PARSER2
  if (((test0->typenode != STRATPARAMINT) && (test0->typenode != STRATPARAMDOUBLE)) ||
      ((test1->typenode != STRATPARAMINT) && (test1->typenode != STRATPARAMDOUBLE))) {
    errorPrint ("stratTestEvalCast: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_PARSER2 */

  if (test0->typenode != test1->typenode) {       /* If value types differ */
    if (test0->typenode == STRATPARAMDOUBLE) {
      test1->typenode        = STRATPARAMDOUBLE;
      test1->data.val.valdbl = (double) test1->data.val.valint;
    }
    else {
      test0->typenode        = STRATPARAMDOUBLE;
      test0->data.val.valdbl = (double) test0->data.val.valint;
    }
  }

  return (0);
}

/* This routine fres the given
** strategy condition.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
stratTestExit (
StratTest * const           test)
{
  int               o;                            /* Output condition flag */

#ifdef SCOTCH_DEBUG_PARSER1
  if (test == NULL) {
    errorPrint ("stratTestExit: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_PARSER1 */

  o = 0;
  switch (test->typetest) {
    case STRATTESTNOT :                           /* Not operator */
      o = stratTestExit (test->data.test[0]);     /* Free the son */
      break;
    case STRATTESTAND :                           /* And operator            */
    case STRATTESTOR  :                           /* Or operator             */
    case STRATTESTLT  :                           /* Less-than operator      */
    case STRATTESTEQ  :                           /* Equal-to operator       */
    case STRATTESTGT  :                           /* Greater-than operator   */
    case STRATTESTMOD :                           /* Modulus operator        */
    case STRATTESTMUL :                           /* Multiplication operator */
    case STRATTESTADD :                           /* Addition operator       */
    case STRATTESTSUB :                           /* Subtraction operator    */
      o  = stratTestExit (test->data.test[0]);    /* Free the sons           */
      o |= stratTestExit (test->data.test[1]);
      break;
    case STRATTESTVAL :                           /* Constant value */
    case STRATTESTVAR :                           /* Variable       */
      break;
#ifdef SCOTCH_DEBUG_PARSER1
    default :
      errorPrint ("stratTestExit: invalid condition type (%u)", test->typetest);
      o = 1;
      break;
#endif /* SCOTCH_DEBUG_PARSER1 */
  }

  memFree (test);                                 /* Free the structure */
  return  (o);
}

/* This routine displays the
** given strategy condition.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static char                 strattestsaveop[STRATTESTNBR] = "|&!=><+-*%##";
static char *               strattestsavepa[2][2] = { { "(", ")" }, { "", "" } };

int
stratTestSave (
const StratTest * const     test,
FILE * const                stream)
{
  int               i;
  int               o;

#ifdef SCOTCH_DEBUG_PARSER1
  if ((test == NULL) || (stream == NULL)) {
    errorPrint ("stratTestSave: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_PARSER1 */

  o = 0;                                          /* Assume no error */
  switch (test->typetest) {
    case STRATTESTNOT :                           /* Not operator */
      if ((fprintf (stream, "!(") == EOF) ||
          (stratTestSave (test->data.test[0], stream) != 0) ||
          (fprintf (stream, ")") == EOF))
        o = 1;
      break;
    case STRATTESTAND :                           /* And operator            */
    case STRATTESTOR :                            /* Or operator             */
    case STRATTESTEQ :                            /* Equal-to operator       */
    case STRATTESTGT :                            /* Greater-than operator   */
    case STRATTESTLT :                            /* Less-than operator      */
    case STRATTESTADD :                           /* Addition operator       */
    case STRATTESTSUB :                           /* Subtraction operator    */
    case STRATTESTMUL :                           /* Multiplication operator */
    case STRATTESTMOD :                           /* Modulus operator        */
      i = (test->data.test[0]->typetest < test->typetest) ? 1 : 0;
      fprintf (stream, "%s", strattestsavepa[i][0]);
      o = stratTestSave (test->data.test[0], stream);
      fprintf (stream, "%s", strattestsavepa[i][1]);
      if (o == 0) {
        fprintf (stream, "%c", strattestsaveop[test->typetest]);
        i = (test->data.test[1]->typetest < test->typetest) ? 1 : 0;
        fprintf (stream, "%s", strattestsavepa[i][0]);
        stratTestSave (test->data.test[1], stream);
        fprintf (stream, "%s", strattestsavepa[i][1]);
      }
      break;
    case STRATTESTVAL :                           /* Constant value */
      switch (test->typenode) {
        case STRATPARAMDOUBLE :
          o = (fprintf (stream, "%lf", test->data.val.valdbl) == EOF);
          break;
        case STRATPARAMINT :
          o = (fprintf (stream, INTSTRING, (INT) test->data.val.valint) == EOF);
          break;
#ifdef SCOTCH_DEBUG_PARSER2
        default :
          errorPrint ("stratTestSave: invalid value type");
          o = 1;
#endif /* SCOTCH_DEBUG_PARSER2 */
      }
      break;
    case STRATTESTVAR :                           /* Variable */
      for (i = 0; test->data.var.datatab->condtab[i].name != NULL; i ++) {
        if ((test->data.var.datatab->condtab[i].dataofft -
             test->data.var.datatab->condtab[i].database) == test->data.var.datadisp)
          break;
      }
      if (test->data.var.datatab->condtab[i].name == NULL) {
        errorPrint ("stratTestSave: invalid variable displacement");
        return     (1);
      }
      o = (fprintf (stream, "%s", test->data.var.datatab->condtab[i].name) == EOF);
      break;
#ifdef SCOTCH_DEBUG_PARSER2
    default :
      errorPrint ("stratTestSave: invalid condition type (%u)", test->typetest);
      o = 1;
#endif /* SCOTCH_DEBUG_PARSER2 */
  }

  return (o);
}
