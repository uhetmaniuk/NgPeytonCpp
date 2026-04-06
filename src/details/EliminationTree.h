#ifndef NGPEYTONCPP_DETAILS_ELIMINATIONTREE_H
#define NGPEYTONCPP_DETAILS_ELIMINATIONTREE_H

#include "details/Utilities.h"

namespace NgPeytonCpp { namespace details { namespace f2c {

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***************     EPOST2 ..... ETREE POSTORDERING #2  *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF THE */
/*       ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE */
/*       CORRESPONDING PARENT AND COLCNT VECTORS ARE ALSO MODIFIED TO */
/*       REFLECT THE REORDERING. */

/*   INPUT PARAMETERS: */
/*       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT */
/*                           IS NEQNS). */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHR VECTOR. */

/*   UPDATED PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR. */
/*       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR. */

/*   OUTPUT PARAMETERS: */
/*       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING. */

/*   WORKING PARAMETERS: */
/*       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE */
/*                           TREE. */

/* *********************************************************************** */

template <typename Index>
void epost2(
  Index root, const Index* fson, Index* brothr, Index* invpos, Index* parent,
  Index* colcnt, Index* stack) {
  /* Local variables */
  Index num, node, itop, ndpar, nunode;

  /* Parameter adjustments */
  --stack;
  --colcnt;
  --parent;
  --invpos;
  --brothr;
  --fson;

  /* Function Body */
  num = 0;
  itop = 0;
  node = root;
  /*       ------------------------------------------------------------- */
  /*       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES */
  /*       ALONG THE TRAVERSAL INTO THE STACK. */
  /*       ------------------------------------------------------------- */
L100:
  ++itop;
  stack[itop] = node;
  node = fson[node];
  if (node > 0) {
    goto L100;
  }
  /*           ---------------------------------------------------------- */
  /*           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT. */
  /*           ---------------------------------------------------------- */
L200:
  if (itop <= 0) {
    goto L300;
  }
  node = stack[itop];
  --itop;
  ++num;
  invpos[node] = num;
  /*               ---------------------------------------------------- */
  /*               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE. */
  /*               ---------------------------------------------------- */
  node = brothr[node];
  if (node <= 0) {
    goto L200;
  }
  goto L100;

L300:
  /*       ------------------------------------------------------------ */
  /*       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR */
  /*       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR. */
  /*       ------------------------------------------------------------ */
  for (node = 1; node <= num; ++node) {
    nunode = invpos[node];
    ndpar = parent[node];
    if (ndpar > 0) {
      ndpar = invpos[ndpar];
    }
    brothr[nunode] = ndpar;
    /* L400: */
  }

  for (nunode = 1; nunode <= num; ++nunode) {
    parent[nunode] = brothr[nunode];
    /* L500: */
  }

  /*       ---------------------------------------------- */
  /*       PERMUTE COLCNT(*) TO REFLECT THE NEW ORDERING. */
  /*       ---------------------------------------------- */
  for (node = 1; node <= num; ++node) {
    nunode = invpos[node];
    stack[nunode] = colcnt[node];
    /* L600: */
  }

  for (node = 1; node <= num; ++node) {
    colcnt[node] = stack[node];
    /* L700: */
  }
} /* epost2 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BTREE2 ..... BINARY TREE REPRESENTATION OF ETREE     ******* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       TO DETERMINE A BINARY TREE REPRESENTATION OF THE ELIMINATION */
/*       TREE, FOR WHICH EVERY "LAST CHILD" HAS THE MAXIMUM POSSIBLE */
/*       COLUMN NONZERO COUNT IN THE FACTOR.  THE RETURNED REPRESENTATION */
/*       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT OF */
/*       THE BINARY TREE IS ALWAYS NEQNS. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */
/*                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF */
/*                           THE ROOTS. */
/*       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR. */

/*   OUTPUT PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */

/*   WORKING PARAMETERS: */
/*       LSON            -   LAST SON VECTOR. */

/* *********************************************************************** */

template <typename Index>
void btree2(
  Index neqns, const Index* parent, const Index* colcnt, Index* fson,
  Index* brothr, Index* lson) {
  /* Local variables */
  Index node, ndpar, lroot, ndlson;

  if (neqns <= 0) {
    return;
  }

  for (node = 0; node < neqns; ++node) {
    fson[node] = 0;
    brothr[node] = 0;
    lson[node] = 0;
  }

  lroot = neqns;
  /*       ------------------------------------------------------------ */
  /*       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING. */
  /*       ------------------------------------------------------------ */
  if (neqns <= 1) {
    return;
  }

  /* Parameter adjustments */
  --lson;
  --brothr;
  --fson;
  --colcnt;
  --parent;

  for (node = neqns - 1; node >= 1; --node) {
    ndpar = parent[node];
    if (ndpar <= 0 || ndpar == node) {
      /*               ------------------------------------------------- */
      /*               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST. */
      /*               SET NODE TO BE ONE OF THE ROOTS OF THE TREES. */
      /*               ------------------------------------------------- */
      brothr[lroot] = node;
      lroot = node;
    } else {
      /*               ------------------------------------------- */
      /*               OTHERWISE, BECOMES FIRST SON OF ITS PARENT. */
      /*               ------------------------------------------- */
      ndlson = lson[ndpar];
      if (ndlson != 0) {
        if (colcnt[node] >= colcnt[ndlson]) {
          brothr[node] = fson[ndpar];
          fson[ndpar] = node;
        } else {
          brothr[ndlson] = node;
          lson[ndpar] = node;
        }
      } else {
        fson[ndpar] = node;
        lson[ndpar] = node;
      }
    }
    /* L300: */
  }
  brothr[lroot] = 0;
} /* btree2 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **********     CHORDR ..... CHILD REORDERING                *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       REARRANGE THE CHILDREN OF EACH VERTEX SO THAT THE LAST ONE */
/*       MAXIMIZES (AMONG THE CHILDREN) THE NUMBER OF NONZEROS IN THE */
/*       CORRESPONDING COLUMN OF L.  ALSO DETERMINE AN NEW POSTORDERING */
/*       BASED ON THE STRUCTURE OF THE MODIFIED ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE NEW */
/*                           POSTORDERING. */
/*       COLCNT          -   COLUMN COUNTS IN L UNDER INITIAL ORDERING; */
/*                           MODIFIED TO REFLECT THE NEW ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BTREE2, EPOST2, INVINV. */

/* *********************************************************************** */

template <typename Index>
void chordr(
  Index neqns, const Index* xadj, const Index* adjncy, Index* perm, Index* invp,
  Index* colcnt, Index* parent, Index* fson, Index* brothr, Index* invpos) {
  (void)xadj;
  (void)adjncy;
  /*       ---------------------------------------------------------- */
  /*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE, */
  /*       SO THAT EACH "LAST CHILD" MAXIMIZES AMONG ITS SIBLINGS THE */
  /*       NUMBER OF NONZEROS IN THE CORRESPONDING COLUMNS OF L. */
  /*       ---------------------------------------------------------- */
  btree2(neqns, parent, colcnt, fson, brothr, invpos);

  /*       ---------------------------------------------------- */
  /*       POSTORDER THE ELIMINATION TREE (USING THE NEW BINARY */
  /*       REPRESENTATION. */
  /*       ---------------------------------------------------- */
  epost2(neqns, fson, brothr, invpos, parent, colcnt, perm);

  /*       -------------------------------------------------------- */
  /*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
  /*       -------------------------------------------------------- */
  invinv(neqns, invp, invpos, perm);
} /* chordr */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***************     ETPOST ..... ETREE POSTORDERING     *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (SEPT 17, 1986) */

/*   PURPOSE: */
/*       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF */
/*       THE ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE */
/*       CORRESPONDING PARENT VECTOR IS ALSO MODIFIED TO REFLECT */
/*       THE REORDERING. */

/*   INPUT PARAMETERS: */
/*       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT */
/*                           IS NEQNS). */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHR VECTOR. */

/*   UPDATED PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR. */

/*   OUTPUT PARAMETERS: */
/*       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING. */

/*   WORKING PARAMETERS: */
/*       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE */
/*                           TREE. */

/* *********************************************************************** */

template <typename Index>
void etpost(
  Index root, const Index* fson, Index* brothr, Index* invpos, Index* parent,
  Index* stack) {
  /* Local variables */
  Index num, node, itop, ndpar, nunode;

  /* Parameter adjustments */
  --invpos;

  /* Function Body */
  num = 0;
  itop = 0;
  node = root;
  /*       ------------------------------------------------------------- */
  /*       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES */
  /*       ALONG THE TRAVERSAL INTO THE STACK. */
  /*       ------------------------------------------------------------- */
L100:
  ++itop;
  stack[itop - 1] = node;
  node = fson[node - 1];
  if (node > 0) {
    goto L100;
  }
  /*           ---------------------------------------------------------- */
  /*           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT. */
  /*           ---------------------------------------------------------- */
L200:
  if (itop <= 0) {
    goto L300;
  }
  node = stack[itop - 1];
  --itop;
  ++num;
  invpos[node] = num;
  /*               ---------------------------------------------------- */
  /*               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE. */
  /*               ---------------------------------------------------- */
  node = brothr[node - 1];
  if (node <= 0) {
    goto L200;
  }
  goto L100;

L300:
  /*       ------------------------------------------------------------ */
  /*       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR */
  /*       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR. */
  /*       ------------------------------------------------------------ */
  for (node = 0; node < num; ++node) {
    nunode = invpos[node + 1];
    ndpar = parent[node];
    if (ndpar > 0) {
      ndpar = invpos[ndpar];
    }
    brothr[nunode - 1] = ndpar;
  }

  for (nunode = 0; nunode < num; ++nunode) {
    parent[nunode] = brothr[nunode];
  }
} /* etpost */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BETREE ..... BINARY TREE REPRESENTATION OF ETREE     ******* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE BINARY TREE REPRESENTATION OF THE ELIMINATION */
/*       TREE GIVEN BY THE PARENT VECTOR.  THE RETURNED REPRESENTATION */
/*       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT */
/*       OF THE BINARY TREE IS ALWAYS NEQNS. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */
/*                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF */
/*                           THE ROOTS. */

/*   OUTPUT PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */

/* *********************************************************************** */

template <typename Index>
void betree(Index neqns, const Index* parent, Index* fson, Index* brothr) {
  /* Local variables */
  Index node, ndpar, lroot;

  /* Function Body */
  if (neqns <= 0) {
    return;
  }

  for (node = 0; node < neqns; ++node) {
    fson[node] = 0;
    brothr[node] = 0;
  }

  lroot = neqns;
  /*       ------------------------------------------------------------ */
  /*       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING. */
  /*       ------------------------------------------------------------ */
  if (neqns <= 1) {
    return;
  }

  /* Parameter adjustments */
  --brothr;
  --fson;
  --parent;

  for (node = neqns - 1; node >= 1; --node) {
    ndpar = parent[node];
    if (ndpar <= 0 || ndpar == node) {
      /*               ------------------------------------------------- */
      /*               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST. */
      /*               SET NODE TO BE ONE OF THE ROOTS OF THE TREES. */
      /*               ------------------------------------------------- */
      brothr[lroot] = node;
      lroot = node;
    } else {
      /*               ------------------------------------------- */
      /*               OTHERWISE, BECOMES FIRST SON OF ITS PARENT. */
      /*               ------------------------------------------- */
      brothr[node] = fson[ndpar];
      fson[ndpar] = node;
    }
    /* L300: */
  }
  brothr[lroot] = 0;
} /* betree */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************     ETREE ..... ELIMINATION TREE     ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND */
/*       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */
/*       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */

/*   WORKING PARAMETERS: */
/*       ANCSTR          -   THE ANCESTOR VECTOR. */

/* *********************************************************************** */

template <typename Index = int64_t>
void etree(
  Index neqns, Index* xadj, Index* adjncy, Index* perm, Index* invp,
  Index* parent, Index* ancstr) {
  /* Local variables */
  Index i, j, nbr, node, next, jstop, jstrt;

  if (neqns <= 0) {
    return;
  }

  /* Parameter adjustments */
  --ancstr;
  --parent;
  --invp;
  --perm;
  --adjncy;
  --xadj;

  for (i = 1; i <= neqns; ++i) {
    parent[i] = 0;
    ancstr[i] = 0;
    node = perm[i];

    jstrt = xadj[node];
    jstop = xadj[node + 1] - 1;
    if (jstrt <= jstop) {
      for (j = jstrt; j <= jstop; ++j) {
        nbr = adjncy[j];
        nbr = invp[nbr];
        if (nbr < i) {
          /*                       ------------------------------------------- */
          /*                       FOR EACH NBR, FIND THE ROOT OF ITS CURRENT */
          /*                       ELIMINATION TREE.  PERFORM PATH COMPRESSION */
          /*                       AS THE SUBTREE IS TRAVERSED. */
          /*                       ------------------------------------------- */
L100:
          if (ancstr[nbr] == i) {
            continue;
          }
          if (ancstr[nbr] > 0) {
            next = ancstr[nbr];
            ancstr[nbr] = i;
            nbr = next;
            goto L100;
          }
          /*  -------------------------------------------- */
          /*    NOW, NBR IS THE ROOT OF THE SUBTREE. */
          /*    MAKE I THE PARENT NODE OF THIS ROOT. */
          /*  ------------------------------------------- */
          parent[nbr] = i;
          ancstr[nbr] = i;
        }
      }
    }
  }
} /* etree */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **********     ETORDR ..... ELIMINATION TREE REORDERING     *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE AN EQUIVALENT REORDERING BASED ON THE STRUCTURE OF */
/*       THE ELIMINATION TREE.  A POSTORDERING OF THE GIVEN ELIMINATION */
/*       TREE IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE EQUIVALENT */
/*                           ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BETREE, ETPOST, ETREE , INVINV. */

/* *********************************************************************** */

template <typename Index = int64_t>
void etordr(
  Index neqns, Index* xadj, Index* adjncy, Index* perm, Index* invp,
  Index* parent, Index* fson, Index* brothr, Index* invpos) {
  /*       ----------------------------- */
  /*       COMPUTE THE ELIMINATION TREE. */
  /*       ----------------------------- */
  etree(neqns, xadj, adjncy, perm, invp, parent, invpos);

  /*       -------------------------------------------------------- */
  /*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE. */
  /*       -------------------------------------------------------- */
  betree(neqns, parent, fson, brothr);

  /*       ------------------------------- */
  /*       POSTORDER THE ELIMINATION TREE. */
  /*       ------------------------------- */
  etpost(neqns, fson, brothr, invpos, parent, perm);

  /*       -------------------------------------------------------- */
  /*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
  /*       -------------------------------------------------------- */
  invinv(neqns, invp, invpos, perm);
} /* etordr */

}}}  // namespace NgPeytonCpp::details::f2c

#endif  // NGPEYTONCPP_DETAILS_ELIMINATIONTREE_H
