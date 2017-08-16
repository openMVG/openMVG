%{
/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : parser_ll.l                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the lexical parser       **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     23 aug 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to     05 jun 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     11 sep 2001     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to     23 dec 2001     **/
/**                # Version 5.1  : from : 09 jun 2009     **/
/**                                 to     24 jul 2011     **/
/**                                                        **/
/**   NOTES      : # In order for flex to read its input   **/
/**                  with getc() instead of fread, we set  **/
/**                  YY_ALWAYS_INTERACTIVE to 1. This may  **/
/**                  not always work with future releases. **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define PARSER_LL

#include "module.h"
#include "common.h"

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

#include "parser.h"
#include "parser_ll.h"
#include "parser_yy.h"
#include "parser_ly.h"                            /*+ Definitions produced by yacc +*/

/* Assume no interactive parsing. */

#ifdef X_OSDOS                                    /* Available only with MKS LEX */
#ifdef YY_INTERACTIVE
#undef YY_INTERACTIVE
#endif /* YY_INTERACTIVE */
#define YY_INTERACTIVE              0
#endif /* X_OSDOS */

#ifdef FLEX_SCANNER
#define YY_ALWAYS_INTERACTIVE       1             /* Set the parser as interactive and read one char at a time */
#define YY_INPUT(buf,result,max_size) { int c = stratParserInput (); result = (c == 0) ? YY_NULL : ((buf)[0] = c, 1); }
#else /* FLEX_SCANNER */
#undef  getc                                      /* Redirect I/O functions */
#define getc                        yygetc
#undef  yygetc
#define yygetc(stream)              stratParserInput ()
#endif /* FLEX_SCANNER */

#define YY_NO_UNPUT                               /* No prototype for yyunput as not defined     */
#define YY_SKIP_YYWRAP                            /* No prototype for yywrap as defined as macro */
#define yywrap()                    (1)           /* Always return end-of-file on end-of-string  */

/*
**  The static variables.
*/

static const char *         stratparserstringptr; /* Pointer to the string to parse */

%}

IDENT                       [A-Za-z][0-9A-Za-z]*
INTEGER                     [0-9]+
FLOAT                       [0-9]+(\.[0-9]+)?([Ee][-+]?[0-9]+)?
STRING                      \"[^\"]*\"

%s lstrat
%s lparam
%s lparamcase
%s lparamdouble
%s lparamint
%s lparamstring
%s ltest

%%
<lstrat>[0-9A-Za-z]         {
                              strncpy (yylval.STRING, yytext, PARSERSTRINGLEN);
                              yylval.STRING[PARSERSTRINGLEN - 1] = '\0';
                              return (METHODNAME);
                            }

<lparam>{IDENT}             {
                              strncpy (yylval.STRING, yytext, PARSERSTRINGLEN);
                              yylval.STRING[PARSERSTRINGLEN - 1] = '\0';
                              return (PARAMNAME);
                            }

<lparamcase>[0-9A-Za-z]     {
                              yylval.CASEVAL = yytext[0];
                              return (VALCASE);
                            }

<lparamdouble>({FLOAT}|{INTEGER}) {
                              yylval.DOUBLE = atof (yytext);
                              return (VALDOUBLE);
                            }

<lparamint>{INTEGER}        {
                              yylval.INTEGER = (INT) atol (yytext);
                              return (VALINT);
                            }
<lparamint>{FLOAT}          {
                              yylval.INTEGER = (INT) atof (yytext); /* FLOAT is put after so that INTEGER can be matched */
                              return (VALINT);
                            }

<lparamstring>{STRING}      {
                              yytext[yyleng - 1] = '\0'; /* Remove the heading and trailing \" */
                              strncpy (yylval.STRING, yytext + 1, PARSERSTRINGLEN);
                              yylval.STRING[PARSERSTRINGLEN - 1] = '\0';
                              return (VALSTRING);
                            }

<ltest>{INTEGER}            {
                              yylval.INTEGER = (INT) atol (yytext);
                              return (VALINT);
                            }
<ltest>{FLOAT}              {
                              yylval.DOUBLE = atof (yytext);
                              return (VALDOUBLE);
                            }
<ltest>{IDENT}              {
                              strncpy (yylval.STRING, yytext, PARSERSTRINGLEN);
                              yylval.STRING[PARSERSTRINGLEN - 1] = '\0';
                              return (PARAMNAME);
                            }

[ \t\n]*                    ;
.                           return (yytext[0]);

%%

/*******************************************/
/*                                         */
/* These routines handle data input to the */
/* lexical analyzer.                       */
/*                                         */
/*******************************************/

/* This routine initializes the
** lexical analyzer.
** It returns:
** - VOID  : in all cases.
*/

void
stratParserInit (
const char * const          string)               /*+ Strategy string to parse +*/
{
#ifdef FLEX_SCANNER
  yyrestart (yyin);                               /* (Re-)initialize the parser */
#endif /* FLEX_SCANNER */

  stratParserSelect (VALSTRAT);                   /* Begin with a strategy        */
  stratparserstringptr = string;                  /* Point to beginning of string */
}

/* This routine reads a single character
** from the input string.
** It returns:
** - 0   : if end of string reached.
** - !0  : character from string.
*/

static
int
stratParserInput ()
{
  if (*stratparserstringptr == '\0')              /* If end-of-string reached  */
    return (0);                                   /* Return end-of-file token  */
  else                                            /* Else return the character */
    return ((int) (unsigned char) *stratparserstringptr ++);
}

/* This routine returns the pointer to the
** remaining part of the string.
*/

const char *
stratParserRemain ()
{
  return (stratparserstringptr);
}

/* This routine selects the sub-parser
** to parse the input.
** It returns:
** - VOID  : in all cases.
*/

void
stratParserSelect (
unsigned int      type)
{
  switch (type) {
    case VALCASE :
      BEGIN lparamcase;
      break;
    case VALDOUBLE :
      BEGIN lparamdouble;
      break;
    case VALINT :
      BEGIN lparamint;
      break;
    case VALSTRING :
      BEGIN lparamstring;
      break;
    case VALPARAM :
      BEGIN lparam;
      break;
    case VALSTRAT :
      BEGIN lstrat;
      break;
    case VALTEST :
      BEGIN ltest;
      break;
  }
}
