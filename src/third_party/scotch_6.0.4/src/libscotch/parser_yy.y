%{
/* Copyright 2004,2007,2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : parser_yy.y                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the syntactic parser     **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     13 jun 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to     27 feb 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to     11 jun 2004     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to     24 jul 2011     **/
/**                # Version 6.0  : from : 30 sep 2014     **/
/**                                 to     30 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define PARSER_YY

#include "module.h"
#include "common.h"

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

#include "parser.h"
#include "parser_ll.h"
#include "parser_yy.h"

/* #define SCOTCH_DEBUG_PARSER3 */
#ifdef SCOTCH_DEBUG_PARSER3
extern int                  yydebug;
#define YYDEBUG                     1
#endif /* SCOTCH_DEBUG_PARSER3 */

/*
**  The static and global definitions.
**  See also at the end of this file.
*/

static const StratTab *     parserstrattab;       /* Pointer to parsing tables          */
static Strat *              parserstratcurr = NULL; /* Pointer to current strategy node */
static StratParamTab *      parserparamcurr = NULL; /* Pointer to current parameter     */

extern unsigned int         parsermethtokentab[]; /* Pre-definition for stupid compilers */

%}

%union {
  char                      CASEVAL;              /* Case value          */
  StratTest *               TEST;                 /* Test type           */
  StratTestType             TESTOP;               /* Relational type     */
  double                    DOUBLE;               /* Double-precision    */
  INT                       INTEGER;              /* Integer             */
  char                      STRING[PARSERSTRINGLEN]; /* Character string */
  struct {
    const StratTab *        tabl;                 /* Current tables    */
    Strat *                 strat;                /* Current method    */
    StratParamTab *         param;                /* Current parameter */
  } SAVE;                                         /* Parameter type    */
  Strat *                   STRAT;                /* Strategy tree     */
}

%token               METHODNAME
%token               PARAMNAME
%token               VALCASE       VALDOUBLE     VALINT        VALSTRING
%token               VALSTRAT      VALPARAM      VALTEST

%type <TEST>         TEST          TESTOR        TESTAND       TESTNOT
%type <TEST>         TESTREL       TESTEXPR1     TESTEXPR2     TESTEXPR3
%type <TEST>         TESTEXPR4     TESTVAL       TESTVAR
%type <TESTOP>       TESTRELOP     TESTEXPR1OP   TESTEXPR2OP   TESTEXPR3OP
%type <CASEVAL>      VALCASE
%type <DOUBLE>       VALDOUBLE     VALSDOUBLE
%type <INTEGER>      VALINT        VALSINT
%type <STRING>       VALSTRING
%type <STRAT>        STRATCONCAT   STRATTEST     STRATTESTELSE STRATEMPTY
%type <STRAT>        STRATGROUP    STRATMETHOD   STRATSELECT
%type <STRING>       METHODNAME    PARAMNAME

%start STRAT

%%

/*
**  These rules define the strategy grammar.
*/

STRAT         : STRATSELECT
              {
                parserstratcurr = ($1);           /* Save pointer to root of tree */
              }
              ;

STRATSELECT   : STRATSELECT '|' STRATEMPTY
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (2)");
                  stratExit  ($1);
                  stratExit  ($3);
                  YYABORT;
                }

                strat->tabl                 = parserstrattab;
                strat->type                 = STRATNODESELECT;
                strat->data.select.strat[0] = ($1);
                strat->data.select.strat[1] = ($3);

                ($$) = strat;
              }
              | STRATEMPTY
              ;

STRATEMPTY    : STRATCONCAT
              |
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (3)");
                  YYABORT;
                }

                strat->tabl = parserstrattab;
                strat->type = STRATNODEEMPTY;

                ($$) = strat;
              }
              ;

STRATCONCAT   : STRATCONCAT STRATTEST
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (4)");
                  stratExit  ($1);
                  stratExit  ($2);
                  YYABORT;
                }

                strat->tabl                 = parserstrattab;
                strat->type                 = STRATNODECONCAT;
                strat->data.concat.strat[0] = ($1);
                strat->data.concat.strat[1] = ($2);

                ($$) = strat;
              }
              | STRATTEST
              ;

STRATTEST     :
              {
                stratParserSelect (VALTEST);      /* Parse parameter tokens */
              }
                '/' TEST
              {
                stratParserSelect (VALSTRAT);     /* Parse strategy tokens */
              }
                '?' STRATSELECT STRATTESTELSE ';'
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint  ("stratParserParse: out of memory (1)");
                  stratExit ($6);
                  if (($7) != NULL)
                    stratExit ($7);
                  stratTestExit ($3);
                  YYABORT;
                }

                strat->tabl               = parserstrattab;
                strat->type               = STRATNODECOND;
                strat->data.cond.test     = ($3);
                strat->data.cond.strat[0] = ($6);
                strat->data.cond.strat[1] = ($7);

                ($$) = strat;
              }
              | STRATGROUP
              ;

STRATTESTELSE : ':' STRATSELECT
              {
                ($$) = ($2);
              }
              |
              {
                ($$) = NULL;
              }
              ;

STRATGROUP    : '(' STRATSELECT ')'
              {
                ($$) = ($2);
              }
              | STRATMETHOD
              ;

STRATMETHOD   : METHODNAME
              {
                Strat *           strat;
                int               meth;
                int               methlen;
                StratMethodTab *  methtab;
                int               i, j;

                meth    =
                methlen = 0;                      /* No method recognized yet   */
                methtab = parserstrattab->methtab; /* Point to the method table */
                for (i = 0; methtab[i].name != NULL; i ++) {
                  if ((strncasecmp (($1),         /* Find longest matching code name */
                       methtab[i].name,
                       j = strlen (methtab[i].name)) == 0) &&
                      (j > methlen)) {
                    meth    = methtab[i].meth;
                    methlen = j;
                  }
                }
                if (methlen == 0) {               /* If method name not known */
                  errorPrint ("stratParserParse: invalid method name \"%s\", before \"%s\"",
                              ($1), stratParserRemain ());
                  YYABORT;
                }
                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (5)");
                  YYABORT;
                }

                strat->tabl             = parserstrattab;
                strat->type             = STRATNODEMETHOD;
                strat->data.method.meth = meth;   /* Set method type         */
                if (methtab[meth].data != NULL)   /* If default values exist */
                  memcpy (&strat->data.method.data, /* Set values to default */
                          methtab[meth].data,
                          sizeof (StratNodeMethodData));

                parserstratcurr = strat;          /* Structure available for parameter processing */
              }
                METHODPARAM
              {
                StratParamTab *   paratab;
                int               i;

                paratab = parserstrattab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].name != NULL; i ++) {
                  if ((paratab[i].meth == parserstratcurr->data.method.meth) && /* If a strategy parameter found for this method */
                      (paratab[i].type == STRATPARAMSTRAT)) {
                    if (*((Strat **) ((byte *) &parserstratcurr->data.method.data + /* And this parameter has not been set */
                        (paratab[i].dataofft - paratab[i].database))) == NULL)
                      errorPrintW ("stratParserParse: strategy parameter \"%s\" of method \"%s\" not set, before \"%s\"",
                                   paratab[i].name, parserstrattab->methtab[parserstratcurr->data.method.meth].name, stratParserRemain ());
                  }
                }

                ($$) = parserstratcurr;           /* Return current structure */
                parserstratcurr = NULL;           /* No current structure     */
              }
              ;

METHODPARAM   :
              {
                stratParserSelect (VALPARAM);     /* Parse parameter tokens */
              }
                '{' PARAMLIST
              {
                stratParserSelect (VALSTRAT);     /* Parse strategy tokens */
              }
                '}'
              |                                   /* No parameters at all */
              ;

PARAMLIST     : PARAMLIST ',' PARAMPARAM
              | PARAMPARAM
              ;

PARAMPARAM    : PARAMNAME
              {
                int               para;
                int               paralen;
                StratParamTab *   paratab;
                int               i, j;

                para    =
                paralen = 0;                      /* No parameter recognized yet   */
                paratab = parserstrattab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].name != NULL; i ++) {
                  if ((paratab[i].meth == parserstratcurr->data.method.meth) &&
                      (strncasecmp (($1),         /* Find longest matching parameter name */
                                    paratab[i].name,
                                    j = strlen (paratab[i].name)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid method parameter name \"%s\", before \"%s\"",
                              ($1), stratParserRemain ());
                  YYABORT;
                }

                ($<SAVE>$).tabl = parserstrattab; /* Save current strategy tables */
                parserparamcurr = &paratab[para]; /* Save current parameter value */
                stratParserSelect (parsermethtokentab[parserparamcurr->type & ~STRATPARAMDEPRECATED]); /* Get non-deprecated type */
                if (parserparamcurr->type == STRATPARAMSTRAT) /* If parameter is a strategy           */
                  parserstrattab = (StratTab *) parserparamcurr->datasltr; /* Use new strategy tables */
              }
                '=' PARAMVAL
              {
                stratParserSelect (VALPARAM);     /* Go-on reading parameters        */
                parserstrattab = ($<SAVE>2).tabl; /* Restore current strategy tables */
              }
              ;

PARAMVAL      : VALCASE
              {
                char              c;              /* Character read             */
                char *            p;              /* Pointer to selector string */
                int               i;              /* Index in selector string   */

                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
                  c = ($1);                       /* First, use char as is                                   */
                  for (p = (char *) parserparamcurr->datasltr, i = 0;
                       (*p != '\0') && (*p != c);
                       p ++, i ++) ;
                  if (*p == '\0') {               /* Char was not found         */
                    c = tolower (c);              /* Convert char to lower case */
                    for (p = (char *) parserparamcurr->datasltr, i = 0;
                         (*p != '\0') && (*p != c);
                         p ++, i ++) ;
                    if (*p == '\0') {
                      errorPrint ("stratParserParse: invalid method parameter switch \"%s=%c\", before \"%s\"",
                                  parserparamcurr->name, ($1), stratParserRemain ());
                      YYABORT;
                    }
                  }

#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (int)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (1)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((int *) ((byte *) &parserstratcurr->data.method.data +
                             (parserparamcurr->dataofft -
                              parserparamcurr->database))) = i;
                }
              }
              | VALSDOUBLE
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (double)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (2)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((double *) ((byte *) &parserstratcurr->data.method.data +
                                (parserparamcurr->dataofft -
                                 parserparamcurr->database))) = ($1);
                }
              }
              | VALSINT
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (INT)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (3)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((INT *) ((byte *) &parserstratcurr->data.method.data +
                             (parserparamcurr->dataofft -
                              parserparamcurr->database))) = (INT) ($1);
                }
              }
              | VALSTRING
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + strlen ($1) + 1) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (4)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  strcpy ((char *) ((byte *) &parserstratcurr->data.method.data +
                                    (parserparamcurr->dataofft -
                                     parserparamcurr->database)),
                          ($1));
                }
              }
              |
              {
                ($<SAVE>$).strat = parserstratcurr;
                ($<SAVE>$).param = parserparamcurr;
                parserstratcurr  = NULL;
                parserparamcurr  = NULL;
              }
                STRATSELECT
              {
                parserstratcurr = ($<SAVE>1).strat; /* Restore current method    */
                parserparamcurr = ($<SAVE>1).param; /* Restore current parameter */

                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (Strat *)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (5)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((Strat **) ((byte *) &parserstratcurr->data.method.data +
                                (parserparamcurr->dataofft -
                                 parserparamcurr->database))) = ($2);
                }
              }
              | error
              {
                errorPrint ("stratParserParse: invalid value for parameter \"%s\" of method \"%s\", before \"%s\"",
                            parserparamcurr->name, parserstratcurr->tabl->methtab[parserstratcurr->data.method.meth].name, stratParserRemain ());
                YYABORT;
              }
              ;

TEST          : TESTOR
              ;

TESTOR        : TESTOR '|' TESTAND
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (6)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }

                test->typetest     = STRATTESTOR;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              | TESTAND
              ;

TESTAND       : TESTAND '&' TESTNOT
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (7)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }

                test->typetest     = STRATTESTAND;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              | TESTNOT
              ;

TESTNOT       : '!' TESTNOT
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (8)");
                  stratTestExit ($2);
                  YYABORT;
                }

                test->typetest     = STRATTESTNOT;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ($2);

                ($$) = test;
              }
              | '(' TESTOR ')'
              {
                ($$) = ($2);
              }
              | TESTREL
              ;

TESTREL       : TESTEXPR1 TESTRELOP TESTEXPR1
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (9)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                test->typetest     = ($2);
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              ;

TESTRELOP     : '<'
              {
                ($$) = STRATTESTLT;
              }
              | '='
              {
                ($$) = STRATTESTEQ;
              }
              | '>'
              {
                ($$) = STRATTESTGT;
              }
              ;

TESTEXPR1     : TESTEXPR1 TESTEXPR1OP TESTEXPR2
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (10)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                test->typetest     = ($2);
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              | TESTEXPR2
              ;

TESTEXPR1OP   : '+'
              {
                ($$) = STRATTESTADD;
              }
              | '-'
              {
                ($$) = STRATTESTSUB;
              }
              ;

TESTEXPR2     : TESTEXPR2 TESTEXPR2OP TESTEXPR3
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  stratTestExit ($1);
                  stratTestExit ($3);
                  errorPrint    ("stratParserParse: out of memory (11)");
                  YYABORT;
                }
                test->typetest     = ($2);
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              | TESTEXPR3
              ;

TESTEXPR2OP   : '*'
              {
                ($$) = STRATTESTMUL;
              }
              ;

TESTEXPR3     : TESTEXPR3 TESTEXPR3OP TESTEXPR4
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (12)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                test->typetest     = ($2);
                test->data.test[0] = ($1);
                test->data.test[1] = ($3);

                ($$) = test;
              }
              | TESTEXPR4
              ;

TESTEXPR3OP   : '%'
              {
                ($$) = STRATTESTMOD;
              }
              ;

TESTEXPR4     : '(' TESTEXPR1 ')'
              {
                ($$) = ($2);
              }
              | TESTVAL
              | TESTVAR
              ;

TESTVAL       : VALSDOUBLE
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (13)");
                  YYABORT;
                }

                test->typetest        = STRATTESTVAL;
                test->typenode        = STRATPARAMDOUBLE;
                test->data.val.valdbl = ($1);

                ($$) = test;
              }
              | VALSINT
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (14)");
                  YYABORT;
                }

                test->typetest        = STRATTESTVAL;
                test->typenode        = STRATPARAMINT;
                test->data.val.valint = ($1);

                ($$) = test;
              }
              ;

TESTVAR       : PARAMNAME
              {
                StratTest *       test;
                StratParamTab *   condtab;
                int               para;
                int               paralen;
                int               i, j;

                para    =
                paralen = 0;                      /* No parameter recognized yet */
                condtab = parserstrattab->condtab; /* Point to parameter table   */
                for (i = 0; condtab[i].name != NULL; i ++) {
                  if ((strncasecmp (($1),         /* Find longest matching parameter name */
                                    condtab[i].name,
                                    j = strlen (condtab[i].name)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid graph parameter name \"%s\", before \"%s\"",
                              ($1), stratParserRemain ());
                  YYABORT;
                }

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (15)");
                  YYABORT;
                }

                test->typetest          = STRATTESTVAR;
                test->typenode          = condtab[para].type;
                test->data.var.datatab  = parserstrattab;
                test->data.var.datadisp = condtab[para].dataofft -
                                          condtab[para].database;

                ($$) = test;
              }
              ;

VALSDOUBLE    : TESTEXPR1OP VALDOUBLE
              {
                ($$) = (($1) == STRATTESTSUB) ? - ($2) : ($2);
              }
              | VALDOUBLE
              ;

VALSINT       : TESTEXPR1OP VALINT
              {
                ($$) = (($1) == STRATTESTSUB) ? - ($2) : ($2);
              }
              | VALINT
              ;

%%

/*
**  The static and global definitions (bis).
**  These are put at the end of the file because
**  the token values that they use are not yet
**  defined in the first section of the file.
*/

unsigned int                parsermethtokentab[] = { /* Table for parameter/token type conversion */
                              VALCASE,
                              VALDOUBLE,
                              VALINT,
                              -1,                 /* No logical parameters */
                              VALSTRAT,
                              VALSTRING,
                              -1                  /* One more value to detect array overflow */
                            };

/************************************/
/*                                  */
/* These routines drive the parser. */
/*                                  */
/************************************/

/* This routine is the entry point for
** the strategy parser.
** It returns:
** - !NULL  : pointer to the strategy.
** - NULL   : on error.
*/

Strat *
stratParserParse (
const StratTab * const      strattab,             /*+ Pointer to parsing tables +*/
const char * const          string)               /*+ Strategy string to parse  +*/
{
  yyclearin;                                      /* Reset the parser state */

#ifdef SCOTCH_DEBUG_PARSER3
  yydebug = 1;                                    /* Set debugging if needed */
#endif /* SCOTCH_DEBUG_PARSER3 */

  stratParserInit (string);                       /* Initialize the lexical parser           */
  parserstrattab  = strattab;                     /* Point to the parsing tables             */
  parserstratcurr = NULL;                         /* Clear up the temporary strategy pointer */

  if (stratParserParse2 () != 0) {                /* Parse the strategy string */
    if (parserstratcurr != NULL)
      stratExit (parserstratcurr);
    return (NULL);
  }

  return (parserstratcurr);                       /* Return strategy pointer */
}

/* This routine displays the parser error message.
** It returns:
** - 1  : in all cases.
*/

static
int
stratParserError (
const char * const          errstr)
{
  errorPrint ("stratParserParse: invalid strategy string, before \"%s\"", stratParserRemain ());
  return     (1);
}
