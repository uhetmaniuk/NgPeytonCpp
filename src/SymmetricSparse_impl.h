#ifndef CPPSELINV_SYMMETRICSPARSE_IMPL_H
#define CPPSELINV_SYMMETRICSPARSE_IMPL_H

#include <cmath>
#include <iostream>
#include <memory>
#include <sys/time.h>
#include <vector>

#include "SymmetricSparse.h"


#ifdef METIS
extern void METIS_EdgeND(int *n, int *xadj, int *adj, int *numflag,
                         int *options, int *perm, int *invp);

extern void METIS_NodeND(int *n, int *xadj, int *adj, int *numflag,
                         int *options, int *perm, int *invp);
#endif

namespace NgPeytonCpp {

namespace details {

        double gtimer() {
            timeval tp;
            struct timezone tz;
            gettimeofday(&tp, &tz);
            return 1000.0 * tp.tv_sec + tp.tv_usec / 1000.0;
        } /* gtimer */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ************     ASSMB .... INDEXED ASSEMBLY OPERATION     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS AN INDEXED ASSEMBLY (I.E., SCATTER-ADD) */
/*       OPERATION, ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE */
/*       CHOLESKY CODES. */

/*   INPUT PARAMETERS: */
/*       M               -   NUMBER OF ROWS IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN Y. */
/*       Y               -   BLOCK UPDATE TO BE INCORPORATED INTO FACTOR */
/*                           STORAGE. */
/*       RELIND          -   RELATIVE INDICES FOR MAPPING THE UPDATES */
/*                           ONTO THE TARGET COLUMNS. */
/*       XLNZ            -   POINTERS TO THE START OF EACH COLUMN IN THE */
/*                           TARGET MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   CONTAINS COLUMNS MODIFIED BY THE UPDATE */
/*                           MATRIX. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void assmb(Index m, Index q, Scalar *y, Index *relind, Index *xlnz,
                   Scalar *lnz, Index lda) {
            /* Local variables */
            Index ir, il1, iy1, icol, ycol, lbot1, yoff1;

            yoff1 = 0;
            for (icol = 1; icol <= q; ++icol) {
                ycol = lda - relind[icol - 1];
                lbot1 = xlnz[ycol] - 1;
                for (ir = icol; ir <= m; ++ir) {
                    il1 = lbot1 - relind[ir - 1];
                    iy1 = yoff1 + ir;
                    lnz[il1 - 1] += y[iy1 - 1];
                    y[iy1 - 1] = static_cast<Scalar>(0);
                }
                yoff1 = iy1 - icol;
            }

        } /* assmb */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***********     INVINV ..... CONCATENATION OF TWO INVP     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO PERFORM THE MAPPING OF */
/*           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP */
/*       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION */
/*       VECTOR PERM IS ALSO COMPUTED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR. */

/*   UPDATED PARAMETERS: */
/*       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON */
/*                           OUTPUT, IT CONTAINS THE NEW INVERSE */
/*                           PERMUTATION. */

/*   OUTPUT PARAMETER: */
/*       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS */
/*                           INVP2). */

/* *********************************************************************** */

        template<typename Index>
        void invinv(Index neqns, Index *invp, const Index *invp2, Index *perm) {
            /* Local variables */
            Index i, node, interm;

            for (i = 0; i < neqns; ++i) {
                interm = invp[i];
                invp[i] = invp2[interm - 1];
            }

            for (i = 0; i < neqns; ++i) {
                node = invp[i] - 1;
                perm[node] = i + 1;
            }

        } /* invinv */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************    FSUP2  ..... FIND SUPERNODES #2   ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO */
/*       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE */
/*       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE */
/*       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS. */

/*   ASSUMPTIONS: */
/*       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT */
/*       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   OUTPUT PARAMETERS: */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           SUPERNODE PARTITIONING. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMEBER 22, 1994. */

/* *********************************************************************** */

        template<typename Index>
        void fsup2(Index neqns, Index nsuper, const Index *snode,
                   Index *xsuper) {
            Index kcol, ksup, lstsup;
            /*       ------------------------------------------------- */
            /*       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*). */
            /*       ------------------------------------------------- */
            lstsup = nsuper + 1;
            for (kcol = neqns; kcol >= 1; --kcol) {
                ksup = snode[kcol - 1];
                if (ksup != lstsup) {
                    xsuper[lstsup - 1] = kcol + 1;
                }
                lstsup = ksup;
            }
            xsuper[0] = 1;
        } /* fsup2 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************    FSUP1 ..... FIND SUPERNODES #1    ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE FIRST OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT RETURNS ONLY THE NUMBER OF */
/*       SUPERNODES NSUPER AND THE SUPERNODE MEMBERSHIP VECTOR SNODE(*), */
/*       WHICH IS OF LENGTH NEQNS.  THE VECTORS OF LENGTH NSUPER ARE */
/*       COMPUTED SUBSEQUENTLY BY THE COMPANION ROUTINE FSUP2. */

/*   METHOD AND ASSUMPTIONS: */
/*       THIS ROUTINE USES THE ELIMINATION TREE AND THE FACTOR COLUMN */
/*       COUNTS TO COMPUTE THE SUPERNODE PARTITION; IT ALSO ASSUMES A */
/*       POSTORDERING OF THE ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           FACTOR COLUMN COUNTS: I.E., THE NUMBER OF */
/*                           NONZERO ENTRIES IN EACH COLUMN OF L */
/*                           (INCLUDING THE DIAGONAL ENTRY). */

/*   OUTPUT PARAMETERS: */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMBER 11, 1994. */

/* *********************************************************************** */

        template<typename Index>
        void fsup1(Index neqns, const Index *etpar, const Index *colcnt,
                   Index &nofsub, Index &nsuper, Index *snode) {
            /*       -------------------------------------------- */
            /*       COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION. */
            /*       -------------------------------------------- */
            nsuper = 1;
            snode[0] = 1;
            nofsub = colcnt[0];
            for (Index kcol = 2; kcol <= neqns; ++kcol) {
                if (etpar[kcol - 2] == kcol) {
                    if (colcnt[kcol - 2] == colcnt[kcol - 1] + 1) {
                        snode[kcol - 1] = nsuper;
                        continue;
                    }
                }
                ++nsuper;
                snode[kcol - 1] = nsuper;
                nofsub += colcnt[kcol - 1];
            }
        } /* fsup1 */

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

        template<typename Index>
        void epost2(Index root, const Index *fson, Index *brothr,
                   Index *invpos,
                   Index *parent, Index *colcnt, Index *stack) {

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

        template<typename Index>
        void btree2(Index neqns, const Index *parent, const Index *colcnt, Index *fson,
                   Index *brothr, Index *lson) {
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

        template<typename Index>
        void chordr(Index neqns, const Index *xadj, const Index *adjncy,
                    Index *perm, Index *invp,
                    Index *colcnt, Index *parent, Index *fson, Index *brothr,
                    Index *invpos) {

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

        template<typename Index>
        void etpost(Index root, const Index *fson, Index *brothr,
                    Index *invpos,
                   Index *parent, Index *stack) {

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

        template<typename Index>
        void betree(Index neqns, const Index *parent, Index *fson, Index *brothr) {

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

/*********************************************************************/
/*********************************************************************/

/*     FLO2HO */

/*        This routine converts the lower triangular part of symmetric */
/*        sparse matrix to its full representation (i.e. both */
/*        upper and lower triangular parts are stored. The diagonal */
/*        of the matrix is INCLUDED in newvals in this version. */

/*     Arguments: */

/*     N        INTEGER (INPUT) */
/*              Dimension of the matrix. */

/*     COLPTR   INTEGER array of size N+1. (INPUT) */
/*              The column pointers of the lower triangular input matrix. */

/*     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT) */
/*              The row indices of the lower triangular input matrix. */

/*     NZVALS   COMPLEX*16 array of size COLPTR(N+1)-1. (INPUT) */
/*              The nonzero entries in the lower triangular part of the */
/*              matrix column-wise compressed. */

/*     NEWPTR   INTEGER array of size N+1. (OUTPUT) */
/*              The column pointers of the converted matrix. */

/*     NEWIND   INTEGER array of size H_NNZ. (OUTPUT) */
/*              The row indices of the converted matrix. */

/*     NEWVALS  COMPLEX*16 array of size H_NNZ (OUTPUT) */
/*              The nonzero entries of the upper and lower triangular */
/*              part of the matrix column-wise compressed. */
/* The diagonal of the matrix is not included. */

/*     IP       INTEGER array of size N. (WORK) */
/*              Work array for storing pointers. */

/*     === work array to accumulate the non-zero */
/*         counts for each row === */

        template<typename Scalar, typename Index>
        void flo2ho(Index n,
                    const Index *colptr, const Index *rowind, const Scalar *nzvals,
                    Index *newptr, Index *newind, Scalar *newvals, Index *ip) {

            /* Local variables */
            Index i, j, iend, inew, irow, ibegin;

            for (i = 0; i < n; ++i)
                ip[i] = 0;

            /* Parameter adjustments */
            --ip;

            --newptr;
            --newind;
            --newvals;

            /*     === nonzero count for each row of the */
            /*         lower triangular part */
            /*         (including the diagonal) === */

            for (j = 0; j < n; ++j) {
                for (i = colptr[j]; i < colptr[j + 1]; ++i) {
                    irow = rowind[i - 1];
                    ++ip[irow];
                }
            }

            /*     === nonzero count for each column */
            /*         (excluding the diagonal)      === */

            for (j = 0; j < n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                if (iend > ibegin) {
                    ip[j + 1] += iend - ibegin;
                }
            }

            /*     === compute pointers to the beginning of each column === */

            newptr[1] = 1;
            for (i = 1; i <= n; ++i) {
                newptr[i + 1] = newptr[i] + ip[i];
            }

            for (i = 1; i <= n; ++i)
                ip[i] = newptr[i];

            /*     === copy the upper triangular part === */
            /*            (excluding the diagonal) */

            for (j = 0; j < n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                if (ibegin < iend) {
                    for (i = ibegin + 1; i <= iend; ++i) {
                        irow = rowind[i - 1];
                        newind[ip[irow]] = j + 1;
                        newvals[ip[irow]] = nzvals[i - 1];
                        ++ip[irow];
                    }
                }
            }

            /*     === copy the lower triangular part === */
            /*            (including the diagonal) */
            for (j = 0; j < n; ++j) {
                inew = ip[j + 1];
                for (i = colptr[j]; i < colptr[j + 1]; ++i) {
                    newind[inew] = rowind[i - 1];
                    newvals[inew] = nzvals[i - 1];
                    ++inew;
                }
            }

        } /* flo2ho */

/*********************************************************************/
/*********************************************************************/

/*     Purpose: */

/*        This routine converts the lower triangular representation */
/*        of symmetric sparse matrix pattern to its full representation */
/*        (i.e. both upper and lower triangular indices are stored. */
/*        Diagonal entries are EXCLUDED from NEWPTR and NEWIND. */

/*     Arguments: */

/*     N        INTEGER (INPUT) */
/*              Dimension of the matrix. */

/*     HNNZ     INTEGER (OUTPUT) */
/*              The number of nonzeros in the entire (both upper and lower */
/*              triangular part) of the matrix. */

/*     COLPTR   INTEGER array of size N+1. (INPUT) */
/*              The column pointers of the lower triangular input matrix. */

/*     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT) */
/*              The row indices of the lower triangular input matrix. */

/*     NEWPTR   INTEGER array of size N+1. (OUTPUT) */
/*              The column pointers of the converted matrix. */

/*     NEWIND   INTEGER array of size H_NNZ. (OUTPUT) */
/*              The row indices of the converted matrix. */

/*     IP       INTEGER array of size N. (WORK) */
/*              Work array for pointers. */

        template<typename Index>
        void ilo2ho(Index n, Index &hnnz, const Index *colptr, const Index *rowind,
                    Index *newptr, Index *newind, Index *ip) {

            /* Local variables */
            Index i, j, nnz, iend, inew, irow, ibegin;

            nnz = colptr[n] - 1;

            /*     === work array to accumulate the non-zero */
            /*         counts for each row === */

            for (i = 0; i < n; ++i) {
                ip[i] = 0;
            }

            /* Parameter adjustments */
            --ip;
            --newptr;
            --colptr;
            --rowind;
            --newind;

            /*     === nonzero count for each row of the */
            /*         lower triangular part === */

            for (j = 1; j <= n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                for (i = ibegin; i <= iend; ++i) {
                    irow = rowind[i];
                    ++ip[irow];
                }
            }

            /*     === nonzero count for each column === */

            for (j = 1; j <= n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                if (iend >= ibegin) {
                    ip[j] = ip[j] + iend - ibegin - 1;
                }
            }

            /*     === compute pointers to the beginning of each column === */

            newptr[1] = 1;
            for (i = 1; i <= n; ++i) {
                newptr[i + 1] = newptr[i] + ip[i];
            }
            hnnz = newptr[n + 1] - 1;

            for (i = 1; i <= n; ++i) {
                ip[i] = newptr[i];
            }

            /*     === copy the upper triangular part === */

            for (j = 1; j <= n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                    for (i = ibegin + 1; i <= iend; ++i) {
                        irow = rowind[i];
                        newind[ip[irow]] = j;
                        ++ip[irow];
                    }
            }

            /*     === copy the lower triangular part === */

            for (j = 1; j <= n; ++j) {
                ibegin = colptr[j];
                iend = colptr[j + 1] - 1;
                inew = ip[j];
                    for (i = ibegin + 1; i <= iend; ++i) {
                        newind[inew] = rowind[i];
                        ++inew;
                    }
            }
        } /* ilo2ho */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE */
/*       REQUIRED BY BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */

/*       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE. */
/*       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY */
/*       TIGHT. */

/*       ---------------------------------------- */
/*       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR */
/*       NEEDED BY BLKFCT. */
/*       ---------------------------------------- */
/* *********************************************************************** */

        template<typename Index>
        void fntsiz(Index nsuper, const Index *xsuper, const Index *snode,
                    const Index *xlindx, const Index *lindx,
                    Index &tmpsiz) {

            /* Local variables */
            Index ii, iend, clen, ksup, bound, ncols, width, tsize, ibegin, length,
                    cursup, nxtsup;

            /* Parameter adjustments */
            --lindx;
            --xlindx;
            --snode;
            --xsuper;

            /* Function Body */
            tmpsiz = 0;
            for (ksup = nsuper; ksup >= 1; --ksup) {
                ncols = xsuper[ksup + 1] - xsuper[ksup];
                ibegin = xlindx[ksup] + ncols;
                iend = xlindx[ksup + 1] - 1;
                length = iend - ibegin + 1;
                bound = length * (length + 1) / 2;
                if (bound > tmpsiz) {
                    cursup = snode[lindx[ibegin]];
                    clen = xlindx[cursup + 1] - xlindx[cursup];
                    width = 0;
                    for (ii = ibegin; ii <= iend; ++ii) {
                        nxtsup = snode[lindx[ii]];
                        if (nxtsup == cursup) {
                            ++width;
                            if (ii == iend) {
                                if (clen > length) {
                                    tsize = length * width - (width - 1) * width / 2;
                                    tmpsiz = std::max(tsize, tmpsiz);
                                }
                            }
                        } else {
                            if (clen > length) {
                                tsize = length * width - (width - 1) * width / 2;
                                tmpsiz = std::max(tsize, tmpsiz);
                            }
                            length -= width;
                            bound = length * (length + 1) / 2;
                            if (bound <= tmpsiz) {
                                break;
                            }
                            width = 1;
                            cursup = nxtsup;
                            clen = xlindx[cursup + 1] - xlindx[cursup];
                        }
                    }
                }
            }
        } /* fntsiz */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     ***** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES */
/*       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE */
/*       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           POINTERS IN THE SUPERNODE INDICES. */
/*       CACHSZ          -   CACHE SIZE IN KILO BYTES. */
/*                           IF THERE IS NO CACHE, SET CACHSZ = 0. */

/*   OUTPUT PARAMETERS: */
/*       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE */
/*                           FINE PARTITIONING. */

/* *********************************************************************** */

        template<typename Index>
        void fnsplt(Index neqns, Index nsuper,
                    const Index *xsuper, const Index *xlindx,
                   Index cachsz, Index *split) {
            /* Local variables */
            Index used, ksup, cache, ncols, width, height, curcol, fstcol, lstcol,
                    nxtblk;

            /*       -------------------------------------------- */
            /*       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE. */
            /*       -------------------------------------------- */
            if (cachsz <= 0) {
                cache = 2000000000;
            } else {
                constexpr double tmpf = 1024.0 / 8 * 9.;
                cache = static_cast<Index>(static_cast<double>(cachsz) * tmpf);
            }

            /*       --------------- */
            /*       INITIALIZATION. */
            /*       --------------- */
            for (Index kcol = 0; kcol < neqns; ++kcol) {
                split[kcol] = 0;
            }

            /* Parameter adjustments */
            --split;
            --xlindx;
            --xsuper;

            /*       --------------------------- */
            /*       FOR EACH SUPERNODE KSUP ... */
            /*       --------------------------- */
            for (ksup = 1; ksup <= nsuper; ++ksup) {
                /*           ----------------------- */
                /*           ... GET SUPERNODE INFO. */
                /*           ----------------------- */
                height = xlindx[ksup + 1] - xlindx[ksup];
                fstcol = xsuper[ksup];
                lstcol = xsuper[ksup + 1] - 1;
                width = lstcol - fstcol + 1;
                nxtblk = fstcol;
                /*           -------------------------------------- */
                /*           ... UNTIL ALL COLUMNS OF THE SUPERNODE */
                /*               HAVE BEEN PROCESSED ... */
                /*           -------------------------------------- */
                curcol = fstcol - 1;
                L200:
                /*               ------------------------------------------- */
                /*               ... PLACE THE FIRST COLUMN(S) IN THE CACHE. */
                /*               ------------------------------------------- */
                ++curcol;
                if (curcol < lstcol) {
                    ++curcol;
                    ncols = 2;
                    used = height * 3 - 1;
                    height += -2;
                } else {
                    ncols = 1;
                    used = height << 1;
                    --height;
                }

                /*               -------------------------------------- */
                /*               ... WHILE THE CACHE IS NOT FILLED AND */
                /*                   THERE ARE COLUMNS OF THE SUPERNODE */
                /*                   REMAINING TO BE PROCESSED ... */
                /*               -------------------------------------- */
                while (used + height < cache && curcol < lstcol) {
                    /*                   -------------------------------- */
                    /*                   ... ADD ANOTHER COLUMN TO CACHE. */
                    /*                   -------------------------------- */
                    ++curcol;
                    ++ncols;
                    used += height;
                    --height;
                }
                /*               ------------------------------------- */
                /*               ... RECORD THE NUMBER OF COLUMNS THAT */
                /*                   FILLED THE CACHE. */
                /*               ------------------------------------- */
                split[nxtblk] = ncols;
                ++nxtblk;
                /*               -------------------------- */
                /*               ... GO PROCESS NEXT BLOCK. */
                /*               -------------------------- */
                if (curcol < lstcol) {
                    goto L200;
                }
            }
        } /* fnsplt */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************     FCNTHN  ..... FIND NONZERO COUNTS    *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN */
/*       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM. */

/*       TECHNIQUES: */
/*       1) SUPERNODE DETECTION. */
/*       2) PATH HALVING. */
/*       3) NO UNION BY RANK. */

/*   NOTES: */
/*       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE. */
/*       2) IT ASSUMES NO DIAGONAL ENTRIES IN THE ADJACENCY STRUCTURE, */
/*          I.E., NO SELF-LOOPS. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                           THE ADJACENCY STRUCTURE, EXCLUDING THE */
/*                           DIAGONAL ENTRIES. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH ROW OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                           THE DIAGONAL ENTRIES. */

/*   WORK PARAMETERS: */
/*       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE */
/*                           DISJOINT SETS (I.E., SUBTREES). */
/*       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS LEAF OF EACH ROW SUBTREE. */
/*       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL */
/*                           (DISTANCE FROM THE ROOT). */
/*       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS */
/*                           USED TO COMPUTE COLUMN COUNTS. */
/*       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT. */
/*       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           NUMBER OF CHILDREN. */
/*       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE. */

/*   FIRST CREATED ON    APRIL 12, 1990. */
/*   LAST UPDATED ON     JANUARY 12, 1995. */
/*   LAST UPDATED ON     OCTOBER 20, 1997. */

/* *********************************************************************** */

        template<typename Index>
        void fcnthn(Index neqns, const Index *xadj, const Index *adjncy, const Index *perm,
                   const Index *invp, const Index *etpar,
                   Index *rowcnt, Index *colcnt, Index &nlnz,
                   Index *set, Index *prvlf, Index *level, Index *weight, Index *fdesc,
                   Index *nchild, Index *prvnbr) {

            /* Local variables */
            Index j, k, lca, temp, xsup, last1, last2, lflag, pleaf, hinbr, jstop,
                    jstrt, ifdesc, oldnbr, parent, lownbr;

            /* *********************************************************************** */

            /*       -------------------------------------------------- */
            /*       COMPUTE LEVEL(*), FDESC(*), NCHILD(*). */
            /*       INITIALIZE XSUP, ROWCNT(*), COLCNT(*), */
            /*                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*). */
            /*       -------------------------------------------------- */
            /* Parameter adjustments */
            --prvnbr;
            --prvlf;
            --set;
            --colcnt;
            --rowcnt;
            --etpar;
            --invp;
            --perm;
            --adjncy;
            --xadj;

            level[0] = 0;
            for (k = neqns; k >= 1; --k) {
                rowcnt[k] = 1;
                colcnt[k] = 0;
                set[k] = k;
                prvlf[k] = 0;
                level[k] = level[etpar[k]] + 1;
                weight[k] = 1;
                fdesc[k] = k;
                nchild[k] = 0;
                prvnbr[k] = 0;
            }
            nchild[0] = 0;
            fdesc[0] = 0;
            for (k = 1; k <= neqns; ++k) {
                parent = etpar[k];
                weight[parent] = 0;
                ++nchild[parent];
                ifdesc = fdesc[k];
                if (ifdesc < fdesc[parent]) {
                    fdesc[parent] = ifdesc;
                }
            }
            /*       ------------------------------------ */
            /*       FOR EACH ``LOW NEIGHBOR'' LOWNBR ... */
            /*       ------------------------------------ */
            for (lownbr = 1; lownbr <= neqns; ++lownbr) {
                lflag = 0;
                ifdesc = fdesc[lownbr];
                oldnbr = perm[lownbr];
                jstrt = xadj[oldnbr];
                jstop = xadj[oldnbr + 1] - 1;
                /*           -------------------------------------------- */
                /*           ISOLATED VERTEX IS LEAF IN OWN LEAF SUBTREE. */
                /*           -------------------------------------------- */
                if (jstrt > jstop) {
                    lflag = 1;
                }
                /*           ----------------------------------------------- */
                /*           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ... */
                /*           ----------------------------------------------- */
                for (j = jstrt; j <= jstop; ++j) {
                    hinbr = invp[adjncy[j]];
                    if (hinbr > lownbr) {
                        if (ifdesc > prvnbr[hinbr]) {
                            /*                       ------------------------- */
                            /*                       INCREMENT WEIGHT(LOWNBR). */
                            /*                       ------------------------- */
                            ++weight[lownbr];
                            pleaf = prvlf[hinbr];
                            /*                       ----------------------------------------- */
                            /*                       IF HINBR HAS NO PREVIOUS ``LOW NEIGHBOR'' */
                            /*                       THEN ... */
                            /*                       ----------------------------------------- */
                            if (pleaf == 0) {
                                /*                  ---------------------------------------- */
                                /*     ... ACCUMULATE LOWNBR-->HINBR PATH * LENGTH */
                                /*                               IN ROWCNT(HINBR). */
                                /*                           -----------------------------------------
             */
                                rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[hinbr];
                            } else {
                                /*                           -----------------------------------------
             */
                                /*                           ... OTHERWISE, LCA <-- FIND(PLEAF),
             * WHICH */
                                /*                               IS THE LEAST COMMON ANCESTOR OF
             * PLEAF */
                                /*                               AND LOWNBR. */
                                /*                               (PATH HALVING.) */
                                /*                           -----------------------------------------
             */
                                last1 = pleaf;
                                last2 = set[last1];
                                lca = set[last2];
                                L300:
                                if (lca != last2) {
                                    set[last1] = lca;
                                    last1 = lca;
                                    last2 = set[last1];
                                    lca = set[last2];
                                    goto L300;
                                }
                                /*                           -------------------------------------
             */
                                /*                           ACCUMULATE PLEAF-->LCA PATH LENGTH IN
             */
                                /*                           ROWCNT(HINBR). */
                                /*                           DECREMENT WEIGHT(LCA). */
                                /*                           -------------------------------------
             */
                                rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[lca];
                                --weight[lca];
                            }
                            /*                       ----------------------------------------------
           */
                            /*                       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF
           * HINBR. */
                            /*                       ----------------------------------------------
           */
                            prvlf[hinbr] = lownbr;
                            lflag = 1;
                        }
                        /*                   --------------------------------------------------
         */
                        /*                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.
         */
                        /*                   --------------------------------------------------
         */
                        prvnbr[hinbr] = lownbr;
                    }
                    /* L500: */
                }
                /*           ---------------------------------------------------- */
                /*           DECREMENT WEIGHT ( PARENT(LOWNBR) ). */
                /*           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP). */
                /*           ---------------------------------------------------- */
                parent = etpar[lownbr];
                --weight[parent];
                if (lflag == 1 || nchild[lownbr] >= 2) {
                    xsup = lownbr;
                }
                set[xsup] = parent;
                /* L600: */
            }
            /*       --------------------------------------------------------- */
            /*       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS. */
            /*       --------------------------------------------------------- */
            nlnz = 0;
            for (k = 1; k <= neqns; ++k) {
                temp = colcnt[k] + weight[k];
                colcnt[k] = temp;
                nlnz += temp;
                parent = etpar[k];
                if (parent != 0) {
                    colcnt[parent] += temp;
                }
            }
        } /* fcnthn */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******         IGATHR .... INTEGER GATHER OPERATION      ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A STANDARD INTEGER GATHER */
/*               OPERATION. */

/*     INPUT PARAMETERS - */
/*        KLEN   - LENGTH OF THE LIST OF GLOBAL INDICES. */
/*        LINDX  - LIST OF GLOBAL INDICES. */
/*        INDMAP - INDEXED BY GLOBAL INDICES, IT CONTAINS THE */
/*                 REQUIRED RELATIVE INDICES. */

/*     OUTPUT PARAMETERS - */
/*        RELIND - LIST RELATIVE INDICES. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void igathr(Index klen, const Index *lindx, const Index *indmap,
                    Index *relind) {
            /* Parameter adjustments */
            --indmap;
            /* Function Body */
            for (Index i = 0; i < klen; ++i) {
                relind[i] = indmap[lindx[i]];
            }
        } /* igathr */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */

/*     ------------------------------------------------------ */
/*     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ... */
/*     ------------------------------------------------------ */

        template<typename Scalar, typename Index>
        void inpnv(Index neqns, Index *xadjf, Index *adjf, Scalar *anzf,
                    Index *perm, Index *invp, Index nsuper, Index *xsuper,
                    Index *xlindx, Index *lindx, Index *xlnz, Scalar *lnz, Index iwsiz,
                    Index *offset, Index &iflag) {
            /* System generated locals */
            Index i__4;

            /* Local variables */
            Index i, j, ii, jlen, oldj, last, jsuper;

            /* Function Body */
            iflag = 0;
            if (iwsiz < neqns) {
                iflag = -1;
                std::cerr << "\n";
                std::cerr << "*** INTEGER WORK SPACE = " << iwsiz << "\n";
                std::cerr << "*** IS SMALLER THAN REQUIRED = " << neqns << "\n";
                std::cerr << "\n";
                return;
            }

            /* Parameter adjustments */
            --offset;
            --lnz;
            --xlnz;
            --lindx;
            --xlindx;
            --xsuper;
            --invp;
            --perm;
            --anzf;
            --adjf;
            --xadjf;

            for (jsuper = 1; jsuper <= nsuper; ++jsuper) {

                /*           ---------------------------------------- */
                /*           FOR EACH SUPERNODE, DO THE FOLLOWING ... */
                /*           ---------------------------------------- */

                /*           ----------------------------------------------- */
                /*           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT. */
                /*           ----------------------------------------------- */
                jlen = xlindx[jsuper + 1] - xlindx[jsuper];
                for (ii = xlindx[jsuper]; ii < xlindx[jsuper + 1]; ++ii) {
                    i = lindx[ii];
                    --jlen;
                    offset[i] = jlen;
                }

                for (j = xsuper[jsuper]; j < xsuper[jsuper + 1]; ++j) {
                    /*               ----------------------------------------- */
                    /*               FOR EACH COLUMN IN THE CURRENT SUPERNODE, */
                    /*               FIRST INITIALIZE THE DATA STRUCTURE. */
                    /*               ----------------------------------------- */
                    for (ii = xlnz[j]; ii < xlnz[j + 1]; ++ii) {
                        lnz[ii] = 0;
                    }

                    /*               ----------------------------------- */
                    /*               NEXT INPUT THE INDIVIDUAL NONZEROS. */
                    /*               ----------------------------------- */
                    oldj = perm[j];
                    last = xlnz[j + 1] - 1;
                    for (ii = xadjf[oldj]; ii < xadjf[oldj + 1]; ++ii) {
                        i = invp[adjf[ii]];
                        if (i >= j) {
                            i__4 = last - offset[i];
                            lnz[i__4] = anzf[ii];
                        }
                    }
                }
            }
        } /* inpnv */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******         LDINDX .... LOAD INDEX VECTOR             ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE COMPUTES THE SECOND INDEX VECTOR */
/*               USED TO IMPLEMENT THE DOUBLY-INDIRECT SAXPY-LIKE */
/*               LOOPS THAT ALLOW US TO ACCUMULATE UPDATE */
/*               COLUMNS DIRECTLY INTO FACTOR STORAGE. */

/*     INPUT PARAMETERS - */
/*        JLEN   - LENGTH OF THE FIRST COLUMN OF THE SUPERNODE, */
/*                 INCLUDING THE DIAGONAL ENTRY. */
/*        LINDX  - THE OFF-DIAGONAL ROW INDICES OF THE SUPERNODE, */
/*                 I.E., THE ROW INDICES OF THE NONZERO ENTRIES */
/*                 LYING BELOW THE DIAGONAL ENTRY OF THE FIRST */
/*                 COLUMN OF THE SUPERNODE. */

/*     OUTPUT PARAMETERS - */
/*        INDMAP - THIS INDEX VECTOR MAPS EVERY GLOBAL ROW INDEX */
/*                 OF NONZERO ENTRIES IN THE FIRST COLUMN OF THE */
/*                 SUPERNODE TO ITS POSITION IN THE INDEX LIST */
/*                 RELATIVE TO THE LAST INDEX IN THE LIST.  MORE */
/*                 PRECISELY, IT GIVES THE DISTANCE OF EACH INDEX */
/*                 FROM THE LAST INDEX IN THE LIST. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void ldindx(Index jlen, Index *lindx, Index *indmap) {
            /* Local variables */
            Index j, jsub, curlen;

            /* Parameter adjustments */
            --indmap;

            curlen = jlen;
            for (j = 0; j < jlen; ++j) {
                jsub = lindx[j];
                --curlen;
                indmap[jsub] = curlen;
            }
        } /* ldindx */

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

        template<typename Index = int64_t>
        void etree(Index neqns, Index *xadj, Index *adjncy, Index *perm, Index *invp,
                   Index *parent, Index *ancstr) {
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

        template<typename Index = int64_t>
        void etordr(Index neqns, Index *xadj, Index *adjncy, Index *perm, Index *invp,
                    Index *parent, Index *fson, Index *brothr, Index *invpos) {

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

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDUPD */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* *****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     ************* */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES */
/*        AFTER A MULTIPLE ELIMINATION STEP. */

/*     INPUT PARAMETERS - */
/*        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED */
/*                 NODES (I.E., NEWLY FORMED ELEMENTS). */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - ADJACENCY STRUCTURE. */
/*        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
/*        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) */
/*                 INTEGER. */

/*     UPDATED PARAMETERS - */
/*        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE. */
/*        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE. */
/*        QSIZE  - SIZE OF SUPERNODE. */
/*        LLIST  - WORKING LINKED LIST. */
/*        MARKER - MARKER VECTOR FOR DEGREE UPDATE. */
/*        TAG    - TAG VALUE. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        int mmdupd(Index ehead, Index neqns, const Index *xadj, Index *adjncy, Index delta,
                   Index &mdeg, Index *dhead, Index *dforw, Index *dbakw, Index *qsize,
                   Index *llist, Index *marker, Index maxint, Index &tag) {
            /* Local variables */
            Index i, j, iq2, deg, deg0, node, mtag, link, mdeg0, enode, fnode, nabor,
                    elmnt, istop, jstop, q2head, istrt, jstrt, qxhead;

            /* Parameter adjustments */
            --marker;
            --llist;
            --qsize;
            --dbakw;
            --dforw;
            --dhead;
            --adjncy;
            --xadj;

            /* Function Body */
            mdeg0 = mdeg + delta;
            elmnt = ehead;
            L100:
            /*            ------------------------------------------------------- */
            /*            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING. */
            /*            (RESET TAG VALUE IF NECESSARY.) */
            /*            ------------------------------------------------------- */
            if (elmnt <= 0) {
                return 0;
            }
            mtag = tag + mdeg0;
            if (mtag < maxint) {
                goto L300;
            }
            tag = 1;
            for (i = 1; i <= neqns; ++i) {
                if (marker[i] < maxint) {
                    marker[i] = 0;
                }
            }
            mtag = tag + mdeg0;
            L300:
            /*            --------------------------------------------- */
            /*            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED */
            /*            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN */
            /*            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE */
            /*            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0, */
            /*            NUMBER OF NODES IN THIS ELEMENT. */
            /*            --------------------------------------------- */
            q2head = 0;
            qxhead = 0;
            deg0 = 0;
            link = elmnt;
            L400:
            istrt = xadj[link];
            istop = xadj[link + 1] - 1;
            for (i = istrt; i <= istop; ++i) {
                enode = adjncy[i];
                link = -enode;
                if (enode < 0) {
                    goto L400;
                } else if (enode == 0) {
                    goto L800;
                } else {
                    goto L500;
                }

                L500:
                if (qsize[enode] == 0) {
                    goto L700;
                }
                deg0 += qsize[enode];
                marker[enode] = mtag;
                /*                        ---------------------------------- */
                /*                        IF ENODE REQUIRES A DEGREE UPDATE, */
                /*                        THEN DO THE FOLLOWING. */
                /*                        ---------------------------------- */
                if (dbakw[enode] != 0) {
                    goto L700;
                }
                /*                            --------------------------------------- */
                /*                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS. */
                /*                            --------------------------------------- */
                if (dforw[enode] == 2) {
                    goto L600;
                }
                llist[enode] = qxhead;
                qxhead = enode;
                goto L700;
                L600:
                llist[enode] = q2head;
                q2head = enode;
                L700:;
            }
            L800:
            /*            -------------------------------------------- */
            /*            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING. */
            /*            -------------------------------------------- */
            enode = q2head;
            iq2 = 1;
            L900:
            if (enode <= 0) {
                goto L1500;
            }
            if (dbakw[enode] != 0) {
                goto L2200;
            }
            ++(tag);
            deg = deg0;
            /*                    ------------------------------------------ */
            /*                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR. */
            /*                    ------------------------------------------ */
            istrt = xadj[enode];
            nabor = adjncy[istrt];
            if (nabor == elmnt) {
                nabor = adjncy[istrt + 1];
            }
            /*                    ------------------------------------------------ */
            /*                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT. */
            /*                    ------------------------------------------------ */
            link = nabor;
            if (dforw[nabor] < 0) {
                goto L1000;
            }
            deg += qsize[nabor];
            goto L2100;
            L1000:
            /*                        -------------------------------------------- */
            /*                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT, */
            /*                        DO THE FOLLOWING. */
            /*                        -------------------------------------------- */
            istrt = xadj[link];
            istop = xadj[link + 1] - 1;
            for (i = istrt; i <= istop; ++i) {
                node = adjncy[i];
                link = -node;
                if (node == enode) {
                    goto L1400;
                }
                if (node < 0) {
                    goto L1000;
                } else if (node == 0) {
                    goto L2100;
                } else {
                    goto L1100;
                }

                L1100:
                if (qsize[node] == 0) {
                    goto L1400;
                }
                if (marker[node] >= tag) {
                    goto L1200;
                }
                /*                                ------------------------------------- */
                /*                                CASE WHEN NODE IS NOT YET CONSIDERED. */
                /*                                ------------------------------------- */
                marker[node] = tag;
                deg += qsize[node];
                goto L1400;
                L1200:
                /*                            ---------------------------------------- */
                /*                            CASE WHEN NODE IS INDISTINGUISHABLE FROM */
                /*                            ENODE.  MERGE THEM INTO A NEW SUPERNODE. */
                /*                            ---------------------------------------- */
                if (dbakw[node] != 0) {
                    goto L1400;
                }
                if (dforw[node] != 2) {
                    goto L1300;
                }
                qsize[enode] += qsize[node];
                qsize[node] = 0;
                marker[node] = maxint;
                dforw[node] = -enode;
                dbakw[node] = -(maxint);
                goto L1400;
                L1300:
                /*                            -------------------------------------- */
                /*                            CASE WHEN NODE IS OUTMATCHED BY ENODE. */
                /*                            -------------------------------------- */
                if (dbakw[node] == 0) {
                    dbakw[node] = -(maxint);
                }
                L1400:;
            }
            goto L2100;
            L1500:
            /*                ------------------------------------------------ */
            /*                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING. */
            /*                ------------------------------------------------ */
            enode = qxhead;
            iq2 = 0;
            L1600:
            if (enode <= 0) {
                goto L2300;
            }
            if (dbakw[enode] != 0) {
                goto L2200;
            }
            ++(tag);
            deg = deg0;
            /*                        --------------------------------- */
            /*                        FOR EACH UNMARKED NABOR OF ENODE, */
            /*                        DO THE FOLLOWING. */
            /*                        --------------------------------- */
            istrt = xadj[enode];
            istop = xadj[enode + 1] - 1;
            for (i = istrt; i <= istop; ++i) {
                nabor = adjncy[i];
                if (nabor == 0) {
                    goto L2100;
                }
                if (marker[nabor] >= tag) {
                    goto L2000;
                }
                marker[nabor] = tag;
                link = nabor;
                /*                                ------------------------------ */
                /*                                IF UNELIMINATED, INCLUDE IT IN */
                /*                                DEG COUNT. */
                /*                                ------------------------------ */
                if (dforw[nabor] < 0) {
                    goto L1700;
                }
                deg += qsize[nabor];
                goto L2000;
                L1700:
                /*                                    ------------------------------- */
                /*                                    IF ELIMINATED, INCLUDE UNMARKED */
                /*                                    NODES IN THIS ELEMENT INTO THE */
                /*                                    DEGREE COUNT. */
                /*                                    ------------------------------- */
                jstrt = xadj[link];
                jstop = xadj[link + 1] - 1;
                for (j = jstrt; j <= jstop; ++j) {
                    node = adjncy[j];
                    link = -node;
                    if (node < 0) {
                        goto L1700;
                    } else if (node == 0) {
                        goto L2000;
                    } else {
                        goto L1800;
                    }

                    L1800:
                    if (marker[node] >= tag) {
                        continue;
                    }
                    marker[node] = tag;
                    deg += qsize[node];
                }
                L2000:;
            }
            L2100:
            /*                    ------------------------------------------- */
            /*                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE */
            /*                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY. */
            /*                    ------------------------------------------- */
            deg = deg - qsize[enode] + 1;
            fnode = dhead[deg];
            dforw[enode] = fnode;
            dbakw[enode] = -deg;
            if (fnode > 0) {
                dbakw[fnode] = enode;
            }
            dhead[deg] = enode;
            if (deg < mdeg) {
                mdeg = deg;
            }
            L2200:
            /*                    ---------------------------------- */
            /*                    GET NEXT ENODE IN CURRENT ELEMENT. */
            /*                    ---------------------------------- */
            enode = llist[enode];
            if (iq2 == 1) {
                goto L900;
            }
            goto L1600;
            L2300:
            /*            ----------------------------- */
            /*            GET NEXT ELEMENT IN THE LIST. */
            /*            ----------------------------- */
            tag = mtag;
            elmnt = llist[elmnt];
            goto L100;

        } /* mmdupd */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDINT */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* ***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE */
/*        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE */
/*        ALGORITHM. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - ADJACENCY STRUCTURE. */

/*     OUTPUT PARAMETERS - */
/*        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE. */
/*        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE). */
/*        LLIST  - LINKED LIST. */
/*        MARKER - MARKER VECTOR. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void mmdint(Index neqns, const Index *xadj, Index *dhead, Index *dforw,
                   Index *dbakw, Index *qsize, Index *llist, Index *marker) {
            /* Local variables */
            Index ndeg, node, fnode;

            /* Parameter adjustments */
            --marker;
            --llist;
            --qsize;
            --dhead;

            /* Function Body */
            for (node = 1; node <= neqns; ++node) {
                dhead[node] = 0;
                qsize[node] = 1;
                marker[node] = 0;
                llist[node] = 0;
                /* L100: */
            }
            /*        ------------------------------------------ */
            /*        INITIALIZE THE DEGREE DOUBLY LINKED LISTS. */
            /*        ------------------------------------------ */
            for (node = 0; node < neqns; ++node) {
                ndeg = xadj[node + 1] - xadj[node] + 1;
                fnode = dhead[ndeg];
                dforw[node] = fnode;
                dhead[ndeg] = node + 1;
                if (fnode > 0) {
                    dbakw[fnode] = node + 1;
                }
                dbakw[node] = -ndeg;
            }

        } /* mmdint */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDELM */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* **     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF */
/*        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH */
/*        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO */
/*        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE */
/*        ELIMINATION GRAPH. */

/*     INPUT PARAMETERS - */
/*        MDNODE - NODE OF MINIMUM DEGREE. */
/*        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT) */
/*                 INTEGER. */
/*        TAG    - TAG VALUE. */

/*     UPDATED PARAMETERS - */
/*        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE. */
/*        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE. */
/*        QSIZE  - SIZE OF SUPERNODE. */
/*        MARKER - MARKER VECTOR. */
/*        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        int mmdelm(Index mdnode, const Index *xadj, Index *adjncy, Index *dhead,
                   Index *dforw, Index *dbakw, Index *qsize, Index *llist,
                   Index *marker, Index maxint, Index tag) {
            /* Local variables */
            Index i, j, npv, node, link, rloc, rlmt, nabor, rnode, elmnt, xqnbr, istop,
                    jstop, istrt, jstrt, nxnode, pvnode, nqnbrs;

            /*        ----------------------------------------------- */
            /*        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE. */
            /*        ----------------------------------------------- */
            /* Parameter adjustments */
            --marker;
            --llist;
            --qsize;
            --dbakw;
            --dforw;
            --dhead;
            --adjncy;
            --xadj;

            /* Function Body */
            marker[mdnode] = tag;
            istrt = xadj[mdnode];
            istop = xadj[mdnode + 1] - 1;
            /*        ------------------------------------------------------- */
            /*        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED */
            /*        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION */
            /*        FOR THE NEXT REACHABLE NODE. */
            /*        ------------------------------------------------------- */
            elmnt = 0;
            rloc = istrt;
            rlmt = istop;
            for (i = istrt; i <= istop; ++i) {
                nabor = adjncy[i];
                if (nabor == 0) {
                    goto L300;
                }
                if (marker[nabor] >= tag) {
                    goto L200;
                }
                marker[nabor] = tag;
                if (dforw[nabor] < 0) {
                    goto L100;
                }
                adjncy[rloc] = nabor;
                ++rloc;
                goto L200;
                L100:
                llist[nabor] = elmnt;
                elmnt = nabor;
                L200:;
            }
            L300:
            /*            ----------------------------------------------------- */
            /*            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS. */
            /*            ----------------------------------------------------- */
            if (elmnt <= 0) {
                goto L1000;
            }
            adjncy[rlmt] = -elmnt;
            link = elmnt;
            L400:
            jstrt = xadj[link];
            jstop = xadj[link + 1] - 1;
            for (j = jstrt; j <= jstop; ++j) {
                node = adjncy[j];
                link = -node;
                if (node < 0) {
                    goto L400;
                } else if (node == 0) {
                    goto L900;
                } else {
                    goto L500;
                }
                L500:
                if (marker[node] >= tag || dforw[node] < 0) {
                    goto L800;
                }
                marker[node] = tag;
                /*                            --------------------------------- */
                /*                            USE STORAGE FROM ELIMINATED NODES */
                /*                            IF NECESSARY. */
                /*                            --------------------------------- */
                L600:
                if (rloc < rlmt) {
                    goto L700;
                }
                link = -adjncy[rlmt];
                rloc = xadj[link];
                rlmt = xadj[link + 1] - 1;
                goto L600;
                L700:
                adjncy[rloc] = node;
                ++rloc;
                L800:;
            }
            L900:
            elmnt = llist[elmnt];
            goto L300;
            L1000:
            if (rloc <= rlmt) {
                adjncy[rloc] = 0;
            }
            /*        -------------------------------------------------------- */
            /*        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ... */
            /*        -------------------------------------------------------- */
            link = mdnode;
            L1100:
            istrt = xadj[link];
            istop = xadj[link + 1] - 1;
            for (i = istrt; i <= istop; ++i) {
                rnode = adjncy[i];
                link = -rnode;
                if (rnode < 0) {
                    goto L1100;
                } else if (rnode == 0) {
                    goto L1800;
                } else {
                    goto L1200;
                }
                L1200:
                /*                -------------------------------------------- */
                /*                IF RNODE IS IN THE DEGREE LIST STRUCTURE ... */
                /*                -------------------------------------------- */
                pvnode = dbakw[rnode];
                if (pvnode == 0 || pvnode == -(maxint)) {
                    goto L1300;
                }
                /*                    ------------------------------------- */
                /*                    THEN REMOVE RNODE FROM THE STRUCTURE. */
                /*                    ------------------------------------- */
                nxnode = dforw[rnode];
                if (nxnode > 0) {
                    dbakw[nxnode] = pvnode;
                }
                if (pvnode > 0) {
                    dforw[pvnode] = nxnode;
                }
                npv = -pvnode;
                if (pvnode < 0) {
                    dhead[npv] = nxnode;
                }
                L1300:
                /*                ---------------------------------------- */
                /*                PURGE INACTIVE QUOTIENT NABORS OF RNODE. */
                /*                ---------------------------------------- */
                jstrt = xadj[rnode];
                jstop = xadj[rnode + 1] - 1;
                xqnbr = jstrt;
                for (j = jstrt; j <= jstop; ++j) {
                    nabor = adjncy[j];
                    if (nabor == 0) {
                        goto L1500;
                    }
                    if (marker[nabor] >= tag) {
                        goto L1400;
                    }
                    adjncy[xqnbr] = nabor;
                    ++xqnbr;
                    L1400:;
                }
                L1500:
                /*                ---------------------------------------- */
                /*                IF NO ACTIVE NABOR AFTER THE PURGING ... */
                /*                ---------------------------------------- */
                nqnbrs = xqnbr - jstrt;
                if (nqnbrs > 0) {
                    goto L1600;
                }
                /*                    ----------------------------- */
                /*                    THEN MERGE RNODE WITH MDNODE. */
                /*                    ----------------------------- */
                qsize[mdnode] += qsize[rnode];
                qsize[rnode] = 0;
                marker[rnode] = maxint;
                dforw[rnode] = -(mdnode);
                dbakw[rnode] = -(maxint);
                goto L1700;
                L1600:
                /*                -------------------------------------- */
                /*                ELSE FLAG RNODE FOR DEGREE UPDATE, AND */
                /*                ADD MDNODE AS A NABOR OF RNODE. */
                /*                -------------------------------------- */
                dforw[rnode] = nqnbrs + 1;
                dbakw[rnode] = 0;
                adjncy[xqnbr] = mdnode;
                ++xqnbr;
                if (xqnbr <= jstop) {
                    adjncy[xqnbr] = 0;
                }

                L1700:;
            }
            L1800:
            return 0;

        } /* mmdelm */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDNUM */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* *****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     ************* */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN */
/*        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION */
/*        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE */
/*        MINIMUM DEGREE ORDERING ALGORITHM. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION. */

/*     UPDATED PARAMETERS - */
/*        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT, */
/*                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED */
/*                 INTO THE NODE -INVP(NODE); OTHERWISE, */
/*                 -INVP(NODE) IS ITS INVERSE LABELLING. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE PERMUTATION VECTOR. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void mmdnum(Index neqns, Index *perm, Index *invp, Index *qsize) {
            /* Local variables */
            Index num, node, root, nextf, father, nqsize;

            /* Parameter adjustments */
            --qsize;
            --invp;
            --perm;

            /* Function Body */
            for (node = 1; node <= neqns; ++node) {
                nqsize = qsize[node];
                if (nqsize <= 0) {
                    perm[node] = invp[node];
                }
                if (nqsize > 0) {
                    perm[node] = -invp[node];
                }
                /* L100: */
            }
            /*        ------------------------------------------------------ */
            /*        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING. */
            /*        ------------------------------------------------------ */
            for (node = 1; node <= neqns; ++node) {
                if (perm[node] > 0) {
                    goto L500;
                }
                /*                ----------------------------------------- */
                /*                TRACE THE MERGED TREE UNTIL ONE WHICH HAS */
                /*                NOT BEEN MERGED, CALL IT ROOT. */
                /*                ----------------------------------------- */
                father = node;
                L200:
                if (perm[father] > 0) {
                    goto L300;
                }
                father = -perm[father];
                goto L200;
                L300:
                /*                ----------------------- */
                /*                NUMBER NODE AFTER ROOT. */
                /*                ----------------------- */
                root = father;
                num = perm[root] + 1;
                invp[node] = -num;
                perm[root] = num;
                /*                ------------------------ */
                /*                SHORTEN THE MERGED TREE. */
                /*                ------------------------ */
                father = node;
                L400:
                nextf = -perm[father];
                if (nextf <= 0) {
                    goto L500;
                }
                perm[father] = -root;
                father = nextf;
                goto L400;
                L500:;
            }
            /*        ---------------------- */
            /*        READY TO COMPUTE PERM. */
            /*        ---------------------- */
            for (node = 1; node <= neqns; ++node) {
                num = -invp[node];
                invp[node] = num;
                perm[num] = node;
            }

        } /* mmdnum */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = GENMMD */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* ****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE */
/*        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION */
/*        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE */
/*        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS */
/*        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM */
/*        EXTERNAL DEGREE. */
/*        --------------------------------------------- */
/*        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE */
/*        DESTROYED. */
/*        --------------------------------------------- */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE. */
/*        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
/*        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER */
/*                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING */
/*                 NODES. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE MINIMUM DEGREE ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        NNZL   - NUMBER OF NONZEROS IN LOWER TRIAGULAR FACTOR. */
/*        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE */
/*                 SCHEME. */
/*        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING */
/*                 THE DIAGONAL ENTRY (DEGREE+1). */
/*        NSUPER - NUMBER OF SUPERNODES. */
/*        XSUPER - FIRST COLUMN OF EACH SUPERNODE. */
/*        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN. */

/*     WORKING PARAMETERS - */
/*        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS. */
/*        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK. */
/*        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK. */
/*        QSIZE  - VECTOR FOR SIZE OF SUPERNODES. */
/*        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS. */
/*        MARKER - A TEMPORARY MARKER VECTOR. */

/*     PROGRAM SUBROUTINES - */
/*        MMDELM, MMDINT, MMDNUM, MMDUPD. */

/* *********************************************************************** */

        template<typename Index>
        void genmmd(Index neqns, const Index *xadj, Index *adjncy, Index *invp, Index *perm,
                   Index delta, Index *dhead, Index *qsize, Index *llist,
                   Index *marker, Index maxint, Index &nnzl, Index &nofsub,
                   Index *colcnt, Index &nsuper, Index *xsuper, Index *snode) {
            /* Local variables */
            Index i, cc, tag, num, mdeg, kcol, ehead, mdlmt, mdnode;
            Index fstcol;
            Index nextmd, lstcol, ksuper;

            if (neqns <= 0) {
                return;
            }

            /* Parameter adjustments */
            --snode;
            --xsuper;
            --colcnt;
            --marker;
            --llist;
            --qsize;
            --dhead;
            --perm;
            --invp;
            --adjncy;
            --xadj;

            /*        ------------------------------------------------ */
            /*        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM. */
            /*        ------------------------------------------------ */
            nsuper = 0;
            mmdint(neqns, &xadj[1], &dhead[1], &invp[1], &perm[1], &qsize[1],
                   &llist[1], &marker[1]);

            /*        ---------------------------------------------- */
            /*        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1. */
            /*        ---------------------------------------------- */
            num = 1;

            /*        ----------------------------- */
            /*        ELIMINATE ALL ISOLATED NODES. */
            /*        ----------------------------- */
            nextmd = dhead[1];
            L100:
            if (nextmd <= 0) {
                goto L200;
            }
            mdnode = nextmd;
            nextmd = invp[mdnode];
            marker[mdnode] = maxint;
            invp[mdnode] = -num;
            ++(nsuper);
            xsuper[nsuper] = num;
            colcnt[num] = 1;
            ++num;
            goto L100;

            L200:
            /*        ---------------------------------------- */
            /*        SEARCH FOR NODE OF THE MINIMUM DEGREE. */
            /*        MDEG IS THE CURRENT MINIMUM DEGREE; */
            /*        TAG IS USED TO FACILITATE MARKING NODES. */
            /*        ---------------------------------------- */
            if (num > neqns) {
                goto L1000;
            }
            tag = 1;
            dhead[1] = 0;
            mdeg = 2;
            L300:
            if (dhead[mdeg] > 0) {
                goto L400;
            }
            ++mdeg;
            goto L300;
            L400:
            /*            ------------------------------------------------- */
            /*            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS */
            /*            WHEN A DEGREE UPDATE IS TO BE PERFORMED. */
            /*            ------------------------------------------------- */
            mdlmt = mdeg + delta;
            ehead = 0;

            L500:
            mdnode = dhead[mdeg];
            if (mdnode > 0) {
                goto L600;
            }
            ++mdeg;
            if (mdeg > mdlmt) {
                goto L900;
            }
            goto L500;
            L600:
            /*                ---------------------------------------- */
            /*                REMOVE MDNODE FROM THE DEGREE STRUCTURE. */
            /*                ---------------------------------------- */
            nextmd = invp[mdnode];
            dhead[mdeg] = nextmd;
            if (nextmd > 0) {
                perm[nextmd] = -mdeg;
            }
            ++(nsuper);
            xsuper[nsuper] = num;
            colcnt[num] = mdeg + qsize[mdnode] - 1;
            invp[mdnode] = -num;
            if (num + qsize[mdnode] > neqns) {
                goto L1000;
            }
            /*                ---------------------------------------------- */
            /*                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH */
            /*                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY. */
            /*                ---------------------------------------------- */
            ++tag;
            if (tag < maxint) {
                goto L800;
            }
            tag = 1;
            for (i = 1; i <= neqns; ++i) {
                if (marker[i] < maxint) {
                    marker[i] = 0;
                }
            }
            L800:
            mmdelm(mdnode, &xadj[1], &adjncy[1], &dhead[1], &invp[1], &perm[1],
                   &qsize[1], &llist[1], &marker[1], maxint, tag);
            num += qsize[mdnode];
            llist[mdnode] = ehead;
            ehead = mdnode;
            if (delta >= 0) {
                goto L500;
            }
            L900:
            /*            ------------------------------------------- */
            /*            UPDATE DEGREES OF THE NODES INVOLVED IN THE */
            /*            MINIMUM DEGREE NODES ELIMINATION. */
            /*            ------------------------------------------- */
            if (num > neqns) {
                goto L1000;
            }
            mmdupd(ehead, neqns, &xadj[1], &adjncy[1], delta, mdeg, &dhead[1], &invp[1],
                   &perm[1], &qsize[1], &llist[1], &marker[1], maxint, tag);
            goto L300;

            L1000:
            xsuper[nsuper + 1] = neqns + 1;
            nnzl = 0;
            nofsub = 0;
            for (ksuper = 1; ksuper <= nsuper; ++ksuper) {
                fstcol = xsuper[ksuper];
                lstcol = xsuper[ksuper + 1] - 1;
                cc = colcnt[fstcol];
                nofsub += cc;
                for (kcol = fstcol; kcol <= lstcol; ++kcol) {
                    snode[kcol] = ksuper;
                    colcnt[kcol] = cc;
                    nnzl += cc;
                    --cc;
                }
            }
            mmdnum(neqns, &perm[1], &invp[1], &qsize[1]);

        } /* genmmd */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     ORDMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE CALLS LIU'S MULTIPLE MINIMUM DEGREE */
/*               ROUTINE. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE. */
/*        IWSIZ  - SIZE OF INTEGER WORKING STORAGE. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE MINIMUM DEGREE ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        NNZL   - NUMBER OF NONZEROS IN LOWER TRIANGULAR FACTOR. */
/*        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE */
/*                 SCHEME. */
/*        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING */
/*                 THE DIAGONAL ENTRY (DEGREE+1). */
/*        NSUPER - NUMBER OF SUPERNODES. */
/*        XSUPER - FIRST COLUMN OF EACH SUPERNODE. */
/*        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN. */
/*        SFIFLG - SFIFLG=.F. MEANS SKIP SYMBOLIC FACTORIZATION */
/*                 INITIALIZATION (SFINIT), SFIFLG=.T. MEANS EXECUTE */
/*                 SFINIT. */
/*        IFLAG  - ERROR FLAG. */
/*                   0: SUCCESSFUL ORDERING */
/*                  -1: INSUFFICIENT WORKING STORAGE */
/*                      [IWORK(*)]. */

/*     WORKING PARAMETERS - */
/*        IWORK  - INTEGER WORKSPACE OF LENGTH 4*NEQNS. */

/* *********************************************************************** */

        template<typename Index>
        void ordmmd(Index neqns, const Index *xadj, Index *adjncy, Index *invp,
                    Index *perm, Index iwsiz, Index *iwork, Index &nnzl, Index &nofsub,
                    Index *colcnt, Index &nsuper, Index *xsuper, Index *snode,
                    bool &sfiflg, Index &iflag) {
            /* Local variables */
            Index delta, maxint;

            iflag = 0;
            if (iwsiz < neqns << 2) {
                iflag = -1;
                std::cerr << "\n";
                std::cerr << "*** INTEGER WORK SPACE = " << iwsiz << "\n";
                std::cerr << "*** IS SMALLER THAN REQUIRED = " << neqns << "\n";
                std::cerr << "\n";
                return;
            }

            /*       DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
            /*       MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER */
            /*                (ANY SMALLER ESTIMATE WILL DO) FOR MARKING */
            /*                NODES. */
            delta = 0;
            maxint = 32767;
            genmmd(neqns, xadj, adjncy, invp, perm, delta, iwork,
                   &iwork[neqns], &iwork[(neqns << 1)], &iwork[neqns * 3],
                   maxint, nnzl, nofsub, colcnt, nsuper, xsuper, snode);
            sfiflg = false;
        } /* ordmmd */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING */
/*       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */
/*       CACHSZ          -   CACHE SIZE (IN KBYTES). */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void bfinit(Index neqns, Index nsuper, Index *xsuper, Index *snode,
                    Index *xlindx, Index *lindx, Index cachsz, Index &tmpsiz,
                    Index *split) {

            /*       --------------------------------------------------- */
            /*       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT. */
            /*       --------------------------------------------------- */
            fntsiz(nsuper, xsuper, snode, xlindx, lindx, tmpsiz);

            /*       ------------------------------- */
            /*       PARTITION SUPERNODES FOR CACHE. */
            /*       ------------------------------- */
            fnsplt(neqns, nsuper, xsuper, xlindx, cachsz, split);

        } /* bfinit */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFC2 ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A */
/*       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS */
/*       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE */
/*       ASSOCIATED WITH THE FACTOR MATRIX L. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                               1 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE */
/*                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST. */
/*       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                           CURRENT LINKED LIST OF MERGED INDICES (THE */
/*                           "REACH" SET). */
/*       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES */
/*                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S */
/*                           INDEX SET. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        int symfc2(Index neqns, Index adjlen, Index *xadj, Index *adjncy, Index *perm,
                   Index *invp, Index *colcnt, Index nsuper, Index *xsuper,
                   Index *snode, Index nofsub, Index *xlindx, Index *lindx,
                   Index *xlnz, Index *mrglnk, Index *rchlnk, Index *marker,
                   Index &flag) {
            /* Local variables */
            Index i, knz, head, node, tail, pcol, newi, jptr, kptr, jsup, ksup, psup,
                    nzbeg, nzend, width, nexti, point, jnzbeg, knzbeg, length, jnzend, jwidth,
                    fstcol, knzend, lstcol;

            /* Parameter adjustments */
            --marker;
            --xlnz;
            --snode;
            --colcnt;
            --invp;
            --perm;
            --xadj;
            --adjncy;
            --mrglnk;
            --xlindx;
            --xsuper;
            --lindx;

            /* Function Body */
            flag = 0;
            if (neqns <= 0) {
                return 0;
            }

            /*       --------------------------------------------------- */
            /*       INITIALIZATIONS ... */
            /*           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX. */
            /*           TAIL   : END OF LIST INDICATOR */
            /*                    (IN RCHLNK(*), NOT MRGLNK(*)). */
            /*           MRGLNK : CREATE EMPTY LISTS. */
            /*           MARKER : "UNMARK" THE INDICES. */
            /*       --------------------------------------------------- */
            nzend = 0;
            head = 0;
            tail = neqns + 1;
            point = 1;
            for (i = 1; i <= neqns; ++i) {
                marker[i] = 0;
                xlnz[i] = point;
                point += colcnt[i];
            }
            xlnz[neqns + 1] = point;
            point = 1;
            for (ksup = 1; ksup <= nsuper; ++ksup) {
                mrglnk[ksup] = 0;
                fstcol = xsuper[ksup];
                xlindx[ksup] = point;
                point += colcnt[fstcol];
            }
            xlindx[nsuper + 1] = point;

            /*       --------------------------- */
            /*       FOR EACH SUPERNODE KSUP ... */
            /*       --------------------------- */
            for (ksup = 1; ksup <= nsuper; ++ksup) {

                /*           --------------------------------------------------------- */
                /*           INITIALIZATIONS ... */
                /*               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP. */
                /*               LSTCOL : LAST COLUMN OF SUPERNODE KSUP. */
                /*               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL. */
                /*               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL. */
                /*           --------------------------------------------------------- */
                fstcol = xsuper[ksup];
                lstcol = xsuper[ksup + 1] - 1;
                width = lstcol - fstcol + 1;
                length = colcnt[fstcol];
                knz = 0;
                rchlnk[head] = tail;
                jsup = mrglnk[ksup];

                /*           ------------------------------------------------- */
                /*           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ... */
                /*           ------------------------------------------------- */
                if (jsup > 0) {
                    /*               --------------------------------------------- */
                    /*               COPY THE INDICES OF THE FIRST CHILD JSUP INTO */
                    /*               THE LINKED LIST, AND MARK EACH WITH THE VALUE */
                    /*               KSUP. */
                    /*               --------------------------------------------- */
                    jwidth = xsuper[jsup + 1] - xsuper[jsup];
                    jnzbeg = xlindx[jsup] + jwidth;
                    jnzend = xlindx[jsup + 1] - 1;
                    for (jptr = jnzend; jptr >= jnzbeg; --jptr) {
                        newi = lindx[jptr];
                        ++knz;
                        marker[newi] = ksup;
                        rchlnk[newi] = rchlnk[head];
                        rchlnk[head] = newi;
                    }
                    /*               ------------------------------------------ */
                    /*               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ... */
                    /*               ------------------------------------------ */
                    jsup = mrglnk[jsup];
                    L300:
                    if (jsup != 0 && knz < length) {
                        /*                   ---------------------------------------- */
                        /*                   MERGE THE INDICES OF JSUP INTO THE LIST, */
                        /*                   AND MARK NEW INDICES WITH VALUE KSUP. */
                        /*                   ---------------------------------------- */
                        jwidth = xsuper[jsup + 1] - xsuper[jsup];
                        jnzbeg = xlindx[jsup] + jwidth;
                        jnzend = xlindx[jsup + 1] - 1;
                        nexti = head;
                        for (jptr = jnzbeg; jptr <= jnzend; ++jptr) {
                            newi = lindx[jptr];
                            L400:
                            i = nexti;
                            nexti = rchlnk[i];
                            if (newi > nexti) {
                                goto L400;
                            }
                            if (newi < nexti) {
                                ++knz;
                                rchlnk[i] = newi;
                                rchlnk[newi] = nexti;
                                marker[newi] = ksup;
                                nexti = newi;
                            }
                        }
                        jsup = mrglnk[jsup];
                        goto L300;
                    }
                }
                /*           --------------------------------------------------- */
                /*           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET. */
                /*           "SORT" ITS STRUCTURE INTO THE LINKED LIST, */
                /*           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE */
                /*           LIST. */
                /*           --------------------------------------------------- */
                if (knz < length) {
                    node = perm[fstcol];
                    knzbeg = xadj[node];
                    knzend = xadj[node + 1] - 1;
                    for (kptr = knzbeg; kptr <= knzend; ++kptr) {
                        newi = adjncy[kptr];
                        newi = invp[newi];
                        if (newi > fstcol && marker[newi] != ksup) {
                            /*                       -------------------------------- */
                            /*                       POSITION AND INSERT NEWI IN LIST */
                            /*                       AND MARK IT WITH KCOL. */
                            /*                       -------------------------------- */
                            nexti = head;
                            L600:
                            i = nexti;
                            nexti = rchlnk[i];
                            if (newi > nexti) {
                                goto L600;
                            }
                            ++knz;
                            rchlnk[i] = newi;
                            rchlnk[newi] = nexti;
                            marker[newi] = ksup;
                        }
                    }
                }
                /*           ------------------------------------------------------------ */
                /*           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LIST. */
                /*           ------------------------------------------------------------ */
                if (rchlnk[head] != fstcol) {
                    rchlnk[fstcol] = rchlnk[head];
                    rchlnk[head] = fstcol;
                    ++knz;
                }

                /*           -------------------------------------------- */
                /*           COPY INDICES FROM LINKED LIST INTO LINDX(*). */
                /*           -------------------------------------------- */
                nzbeg = nzend + 1;
                nzend += knz;
                if (nzend + 1 != xlindx[ksup + 1]) {
                    goto L8000;
                }
                i = head;
                for (kptr = nzbeg; kptr <= nzend; ++kptr) {
                    i = rchlnk[i];
                    lindx[kptr] = i;
                    /* L800: */
                }

                /*           --------------------------------------------------- */
                /*           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S */
                /*           "MERGE" LIST. */
                /*           --------------------------------------------------- */
                if (length > width) {
                    pcol = lindx[xlindx[ksup] + width];
                    psup = snode[pcol];
                    mrglnk[ksup] = mrglnk[psup];
                    mrglnk[psup] = ksup;
                }

            }

            return 0;

/*       ----------------------------------------------- */
/*       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED. */
/*       ----------------------------------------------- */
            L8000:
            flag = -2;
            return 0;

        } /* symfc2 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFCT ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE CALLS SYMFC2 WHICH PERFORMS SUPERNODAL SYMBOLIC */
/*       FACTORIZATION ON A REORDERED LINEAR SYSTEM. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */
/*       (I) IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                              -1 - INSUFFICIENT INTEGER WORKING SPACE. */
/*                              -2 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) IWORK       -   WORKING ARRAY OF LENGTH NSUPER+2*NEQNS. */

/* *********************************************************************** */

        template<typename Index = int64_t>
        void symfct(Index neqns, Index adjlen, Index *xadj,
                    Index *adjncy, Index *perm, Index *invp, Index *colcnt,
                    Index nsuper, Index *xsuper, Index *snode, Index nofsub,
                    Index *xlindx, Index *lindx, Index *xlnz, Index iwsiz, Index *iwork,
                    Index &flag) {
            flag = 0;
            if (iwsiz < nsuper + (neqns << 1) + 1) {
                flag = -1;
                std::cerr << "\n";
                std::cerr << "*** INTEGER WORK STORAGE = " << iwsiz << "\n";
                std::cerr << "*** IS SMALLER THAN REQUIRED = "
                          << nsuper + (neqns << 1) + 1 << "\n";
                std::cerr << "\n";
                return;
            }
            symfc2(neqns, adjlen, xadj, adjncy, perm, invp, colcnt,
                   nsuper, xsuper, snode, nofsub, xlindx, lindx, xlnz,
                   iwork, &iwork[nsuper], &iwork[nsuper + neqns + 1], flag);
            if (flag == -2) {
                std::cerr << "\n";
                std::cerr << "*** INCONSISTENCY IN THE INPUT\n";
                std::cerr << "*** TO SYMBOLIC FACTORIZATION\n";
                std::cerr << "\n";
                return;
            }
        } /* symfct */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     DSCAL .... SCALE A VECTOR                     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE COMPUTES A <-- AX, WHERE A IS A */
/*               SCALAR AND X IS A VECTOR. */

/*     INPUT PARAMETERS - */
/*        N - LENGTH OF THE VECTOR X. */
/*        A - SCALAR MULIPLIER. */
/*        X - VECTOR TO BE SCALED. */

/*     OUTPUT PARAMETERS - */
/*        X - REPLACED BY THE SCALED VECTOR, AX. */

/* *********************************************************************** */

        template<typename Scalar>
        void scal(int64_t n, Scalar a, Scalar *x) {
            for (int64_t i = 0; i < n; ++i) {
                x[i] = a * x[i];
            }
        } /* scal */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     SMXPY1 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '1' SIGNIFIES NO LOOP UNROLLING, I.E., */
/*               LOOP-UNROLLING TO LEVEL 1. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void smxpy1(Index m, Index n, Scalar *y, Index *apnt, Scalar *a) {
            /* System generated locals */
            Index i__3, i__4, i__5;
            Scalar z__2;

            /* Local variables */
            Index i, j, ii, jj;
            Scalar amult;

            /* Parameter adjustments */
            --y;
            --apnt;
            --a;

            /* Function Body */
            for (j = 1; j <= n; ++j) {
                jj = apnt[j];
                ii = apnt[j + 1] - m;
                z__2 = -a[jj];
                i__3 = ii;
                amult = z__2 * a[i__3];
                for (i = 1; i <= m; ++i) {
                    i__3 = i;
                    i__4 = i;
                    i__5 = ii;
                    z__2 = amult * a[i__5];
                    y[i__3] = y[i__4] + z__2;
                    ++ii;
                }
            }
        } /* smxpy1 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     PCHOL .... DENSE PARTIAL CHOLESKY             ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY */
/*               FACTORIZATION ON THE COLUMNS OF A SUPERNODE */
/*               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS */
/*               EXTERNAL TO THE SUPERNODE. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN). */
/*        N      - NUMBER OF COLUMNS IN THE SUPERNODE. */
/*        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END */
/*                 OF THE J-TH COLUMN OF THE SUPERNODE. */
/*        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO */
/*                 BE FACTORED. */
/*        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*     OUTPUT PARAMETERS - */
/*        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF */
/*                 THE SUPERNODE. */
/*        IFLAG  - UNCHANGED IF THERE IS NO ERROR. */
/*                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void pchol(Index m, Index n, Index *xpnt, Scalar *x, Index &iflag) {
            /* Local variables */
            Index mm;
            Scalar diag;
            Index jcol, jpnt;

            /*       ------------------------------------------ */
            /*       FOR EVERY COLUMN JCOL IN THE SUPERNODE ... */
            /*       ------------------------------------------ */
            /* Parameter adjustments */
            --x;
            --xpnt;

            /* Function Body */
            mm = m;
            jpnt = xpnt[1];
            for (jcol = 1; jcol <= n; ++jcol) {

                /*           ---------------------------------- */
                /*           UPDATE JCOL WITH PREVIOUS COLUMNS. */
                /*           ---------------------------------- */
                if (jcol > 1) {
                    smxpy1(mm, jcol - 1, &x[jpnt], &xpnt[1], &x[1]);
                }

                /*           --------------------------- */
                /*           COMPUTE THE DIAGONAL ENTRY. */
                /*           --------------------------- */
                diag = x[jpnt];
                if (diag == static_cast<Scalar>(0.0)) {
                    iflag = 1;
                    return;
                }
                diag = static_cast<Scalar>(1.0) / diag;

                /*           ---------------------------------------------------- */
                /*           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY. */
                /*           ---------------------------------------------------- */
                --mm;
                ++jpnt;
                scal(mm, diag, &x[jpnt]);
                jpnt += mm;
            }
        } /* pchol */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     MMPY1  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       LOOP UNROLLING: LEVEL 1 */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        int mmpy1(Index m, Index n, Index q, Index *xpnt, Scalar *x, Scalar *y,
                   Index ldy) {
            /* System generated locals */
            Index i__4, i__5, i__6;
            Scalar z__1, z__2;

            /* Local variables */
            Scalar a1;
            Index i1, j1, mm, iy, xcol, ycol, leny, iylast, iystop, iystrt;

            /* Parameter adjustments */
            --y;
            --x;
            --xpnt;

            /* Function Body */
            mm = m;
            iylast = 0;
            leny = ldy;
            /*       ------------------------------------ */
            /*       TO COMPUTE EACH COLUMN YCOL OF Y ... */
            /*       ------------------------------------ */
            for (ycol = 1; ycol <= q; ++ycol) {
                iystrt = iylast + 1;
                iystop = iystrt + mm - 1;
                iylast += leny;
                /*           -------------------------------------------------- */
                /*           ... PERFORM THE APPROPRIATE MATRIX VECTOR MULTIPLY: */
                /*               X * A(*,YCOL). */
                /*           -------------------------------------------------- */
                for (xcol = 1; xcol <= n; ++xcol) {
                    j1 = xpnt[xcol];
                    i1 = xpnt[xcol + 1] - mm;
                    z__2 = -x[j1];
                    z__1 = z__2 * x[i1];
                    a1 = z__1;
                    for (iy = iystrt; iy <= iystop; ++iy) {
                        i__4 = iy;
                        i__5 = iy;
                        i__6 = i1;
                        z__2 = a1 * x[i__6];
                        z__1 = y[i__5] + z__2;
                        y[i__4] = z__1;
                        ++i1;
                    }
                }
                --mm;
                --leny;
            }

            return 0;
        } /* mmpy1 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     CHLSUP .... DENSE CHOLESKY WITHIN SUPERNODE   ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY */
/*               FACTORIZATION ON THE COLUMNS OF A SUPERNODE */
/*               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS */
/*               EXTERNAL TO THE SUPERNODE. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN). */
/*        N      - NUMBER OF COLUMNS IN THE SUPERNODE. */
/*        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END */
/*                 OF THE J-TH COLUMN OF THE SUPERNODE. */
/*        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO */
/*                 BE FACTORED. */
/*        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*     OUTPUT PARAMETERS - */
/*        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF */
/*                 THE SUPERNODE. */
/*        IFLAG  - UNCHANGED IF THERE IS NO ERROR. */
/*                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void chlsup(Index m, Index n, Index *split, Index *xpnt, Scalar *x,
                   Index &iflag) {
            Index q, mm, nn, jblk, jpnt;
            Index fstcol, nxtcol;

            /* Parameter adjustments */
            --x;
            --xpnt;
            --split;

            /* Function Body */
            jblk = 0;
            fstcol = 1;
            mm = m;
            jpnt = xpnt[fstcol];

            /*       ---------------------------------------- */
            /*       FOR EACH BLOCK JBLK IN THE SUPERNODE ... */
            /*       ---------------------------------------- */
            while (fstcol <= n) {
                ++jblk;
                nn = split[jblk];
                /*           ------------------------------------------ */
                /*           ... PERFORM PARTIAL CHOLESKY FACTORIZATION */
                /*               ON THE BLOCK. */
                /*           ------------------------------------------ */
                pchol(mm, nn, &xpnt[fstcol], &x[1], iflag);
                if (iflag == 1) {
                    return;
                }
                /*           ---------------------------------------------- */
                /*           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS */
                /*               OF THE SUPERNODE REMAINING TO BE COMPUTED. */
                /*           ---------------------------------------------- */
                nxtcol = fstcol + nn;
                q = n - nxtcol + 1;
                mm -= nn;
                jpnt = xpnt[nxtcol];
                if (q > 0) {
                    mmpy1(mm, nn, q, &xpnt[fstcol], &x[1], &x[jpnt], mm);
                }
                fstcol = nxtcol;
            }
        } /* chlsup */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************     MMPY  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       SPLIT(*)        -   BLOCK PARTITIONING OF X. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void mmpy(Index m, Index n, Index q, Index *split, Index *xpnt, Scalar *x,
                  Scalar *y, Index ldy) {
            Index nn, blk, fstcol;
            blk = 0;
            fstcol = 0;
            while (fstcol < n) {
                nn = split[blk];
                mmpy1(m, nn, q, &xpnt[fstcol], x, y, ldy);
                fstcol += nn;
                ++blk;
            }
        } /* mmpy */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     MMPYI  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       MATRIX X HAS ONLY 1 COLUMN. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       D               -   DIAGONAL ENTRY IN COLUMN X BACK IN CHOLESKY. */
/*       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN */
/*       RELIND(*)       -   RELATIVE INDICES. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void mmpyi_(Index m, Index q, Index *xpnt, Scalar *x, Scalar d, Index *iy,
                    Scalar *y, Index *relind) {
            /* Local variables */
            Scalar a, z1, z2;
            Index i, k, col, isub, ylast;

            /* Parameter adjustments */
            --relind;
            --y;
            --iy;
            --x;
            --xpnt;

            /* Function Body */
            for (k = 1; k <= q; ++k) {
                col = xpnt[k];
                ylast = iy[col + 1] - 1;
                z2 = -(d);
                z1 = z2 * x[k];
                a = z1;
                for (i = k; i <= m; ++i) {
                    isub = xpnt[i];
                    isub = ylast - relind[isub];
                    y[isub] += a * x[i];
                }
            }
        } /* mmpyi_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKSLV ... BLOCK TRIANGULAR SOLUTIONS          ********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC */
/*       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE */
/*       TRIANGULAR SOLUTION.  IT USES OUTPUT FROM BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE. */
/*       (XLNZ,LNZ)      -   CHOLESKY FACTOR. */

/*   UPDATED PARAMETERS: */
/*       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON */
/*                           OUTPUT, CONTAINS THE SOLUTION. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void blkslv_(Index nsuper, Index *xsuper, Index *xlindx, Index *lindx,
                     Index *xlnz, Scalar *lnz, Scalar *rhs) {
            /* System generated locals */
            Index i__3, i__4, i__5, i__6;
            Scalar z__1, z__2;

            /* Local variables */
            Index i;
            Scalar t;
            Index ix, jcol, ipnt, jpnt, jsup, fjcol, ljcol, ixstop, ixstrt;

            /* Parameter adjustments */
            --rhs;
            --lnz;
            --xlnz;
            --lindx;
            --xlindx;
            --xsuper;

            /* Function Body */
            if (nsuper <= 0) {
                return;
            }

            /*       ------------------------ */
            /*       FORWARD SUBSTITUTION ... */
            /*       ------------------------ */
            fjcol = xsuper[1];
            for (jsup = 1; jsup <= nsuper; ++jsup) {
                ljcol = xsuper[jsup + 1] - 1;
                ixstrt = xlnz[fjcol];
                jpnt = xlindx[jsup];
                for (jcol = fjcol; jcol <= ljcol; ++jcol) {
                    ixstop = xlnz[jcol + 1] - 1;
                    t = rhs[jcol];
                    ipnt = jpnt + 1;
                    for (ix = ixstrt + 1; ix <= ixstop; ++ix) {
                        i = lindx[ipnt];
                        i__4 = i;
                        i__5 = i;
                        i__6 = ix;
                        z__2 = t * lnz[i__6];
                        z__1 = rhs[i__5] - z__2;
                        rhs[i__4] = z__1;
                        ++ipnt;
                    }
                    ixstrt = ixstop + 1;
                    ++jpnt;
                }
                fjcol = ljcol + 1;
            }

            /*       ------------------ */
            /*       DIAGONAL SOLVE ... */
            /*       ------------------ */
            for (jcol = 1; jcol < xsuper[nsuper + 1]; ++jcol) {
                rhs[jcol] = rhs[jcol] / lnz[xlnz[jcol]];
            }

            /*       ------------------------- */
            /*       BACKWARD SUBSTITUTION ... */
            /*       ------------------------- */
            ljcol = xsuper[nsuper + 1] - 1;
            for (jsup = nsuper; jsup >= 1; --jsup) {
                fjcol = xsuper[jsup];
                ixstop = xlnz[ljcol + 1] - 1;
                jpnt = xlindx[jsup] + (ljcol - fjcol);
                for (jcol = ljcol; jcol >= fjcol; --jcol) {
                    ixstrt = xlnz[jcol];
                    ipnt = jpnt + 1;
                    t = rhs[jcol];
                    for (ix = ixstrt + 1; ix <= ixstop; ++ix) {
                        i = lindx[ipnt];
                        i__3 = ix;
                        i__4 = i;
                        z__2 = lnz[i__3] * rhs[i__4];
                        z__1 = t - z__2;
                        t = z__1;
                        ++ipnt;
                    }
                    rhs[jcol] = t;
                    ixstop = ixstrt - 1;
                    --jpnt;
                }
                ljcol = fjcol - 1;
            }
        } /* blkslv_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFC2 .....  BLOCK GENERAL SPARSE CHOLESKY     ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE FACTORS A SPARSE POSITIVE DEFINITE MATRIX. */
/*       THE COMPUTATION IS ORGANIZED AROUND KERNELS THAT PERFORM */
/*       SUPERNODE-TO-SUPERNODE UPDATES, I.E., BLOCK-TO-BLOCK UPDATES. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING */
/*                           IT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING */
/*                           THE DIAGONAL ELEMENTS). */
/*       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED. */
/*       TMPSIZ          -   SIZE OF TEMPORARY WORKING STORAGE. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY. */
/*       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   ON OUTPUT, CONTAINS CHOLESKY FACTOR. */
/*       IFLAG           -   ERROR FLAG. */
/*                               0: SUCCESSFUL FACTORIZATION. */
/*                              -1: NONPOSITIVE DIAGONAL ENCOUNTERED, */
/*                                  MATRIX IS NOT POSITIVE DEFINITE. */
/*                              -2: INSUFFICIENT WORKING STORAGE */
/*                                  [TEMP(*)]. */

/*   WORKING PARAMETERS: */
/*       LINK            -   LINKS TOGETHER THE SUPERNODES IN A SUPERNODE */
/*                           ROW. */
/*       LENGTH          -   LENGTH OF THE ACTIVE PORTION OF EACH */
/*                           SUPERNODE. */
/*       INDMAP          -   VECTOR OF SIZE NEQNS INTO WHICH THE GLOBAL */
/*                           INDICES ARE SCATTERED. */
/*       RELIND          -   MAPS LOCATIONS IN THE UPDATING COLUMNS TO */
/*                           THE CORRESPONDING LOCATIONS IN THE UPDATED */
/*                           COLUMNS.  (RELIND IS GATHERED FROM INDMAP). */
/*       TEMP            -   REAL VECTOR FOR ACCUMULATING UPDATES.  MUST */
/*                           ACCOMODATE ALL COLUMNS OF A SUPERNODE. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void blkfc2(Index nsuper, Index *xsuper, Index *snode, Index *split,
                   Index *xlindx, Index *lindx, Index *xlnz, Scalar *lnz, Index *link,
                   Index *length, Index *indmap, Index *relind, Index tmpsiz,
                   Scalar *temp, Index &iflag) {
            /* Local variables */
            Index i, ilen, jlen, klen, jsup, ksup;
            Index fjcol, fkcol, ljcol;
            Index klast, kdpnt, ilpnt, jlpnt, klpnt, store;
            Index jxpnt, kxpnt, inddif;
            Index njcols, nkcols;
            Index ncolup, kfirst, nxtcol, nxksup, nxtsup;

            /* Parameter adjustments */
            --temp;
            --relind;
            --indmap;
            --length;
            --link;
            --lnz;
            --xlnz;
            --lindx;
            --xlindx;
            --split;
            --snode;
            --xsuper;

            /* Function Body */
            iflag = 0;

            /*       ----------------------------------------------------------- */
            /*       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*). */
            /*       ----------------------------------------------------------- */
            for (jsup = 1; jsup <= nsuper; ++jsup) {
                link[jsup] = 0;
            }
            for (i = 1; i <= tmpsiz; ++i) {
                temp[i] = 0;
            }

            /*       --------------------------- */
            /*       FOR EACH SUPERNODE JSUP ... */
            /*       --------------------------- */
            for (jsup = 1; jsup <= nsuper; ++jsup) {

                /*           ------------------------------------------------ */
                /*           FJCOL  ...  FIRST COLUMN OF SUPERNODE JSUP. */
                /*           LJCOL  ...  LAST COLUMN OF SUPERNODE JSUP. */
                /*           NJCOLS ...  NUMBER OF COLUMNS IN SUPERNODE JSUP. */
                /*           JLEN   ...  LENGTH OF COLUMN FJCOL. */
                /*           JXPNT  ...  POINTER TO INDEX OF FIRST */
                /*                       NONZERO IN COLUMN FJCOL. */
                /*           ------------------------------------------------ */
                fjcol = xsuper[jsup];
                njcols = xsuper[jsup + 1] - fjcol;
                ljcol = fjcol + njcols - 1;
                jlen = xlnz[fjcol + 1] - xlnz[fjcol];
                jxpnt = xlindx[jsup];

                /*           ----------------------------------------------------- */
                /*           SET UP INDMAP(*) TO MAP THE ENTRIES IN UPDATE COLUMNS */
                /*           TO THEIR CORRESPONDING POSITIONS IN UPDATED COLUMNS, */
                /*           RELATIVE THE THE BOTTOM OF EACH UPDATED COLUMN. */
                /*           ----------------------------------------------------- */
                ldindx(jlen, &lindx[jxpnt], &indmap[1]);

                /*           ----------------------------------------- */
                /*           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ... */
                /*           ----------------------------------------- */
                ksup = link[jsup];
                L300:
                if (ksup > 0) {
                    nxksup = link[ksup];

                    /*               -------------------------------------------------------
       */
                    /*               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE. */

                    /*               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP. */
                    /*               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP. */
                    /*               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FKCOL. */
                    /*               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN ACTIVE
       */
                    /*                           PORTION OF COLUMN FJCOL. */
                    /*               -------------------------------------------------------
       */
                    fkcol = xsuper[ksup];
                    nkcols = xsuper[ksup + 1] - fkcol;
                    klen = length[ksup];
                    kxpnt = xlindx[ksup + 1] - klen;

                    /*               ------------------------------------------- */
                    /*               PERFORM CMOD(JSUP,KSUP), WITH SPECIAL CASES */
                    /*               HANDLED DIFFERENTLY. */
                    /*               ------------------------------------------- */

                    if (klen != jlen) {

                        /*                   ------------------------------------------- */
                        /*                   SPARSE CMOD(JSUP,KSUP). */

                        /*                   NCOLUP ... NUMBER OF COLUMNS TO BE UPDATED. */
                        /*                   ------------------------------------------- */

                        for (i = 0; i < klen; ++i) {
                            nxtcol = lindx[kxpnt + i];
                            if (nxtcol > ljcol) {
                                goto L500;
                            }
                        }
                        i = klen;
                        L500:
                        ncolup = i;

                        if (nkcols == 1) {

                            /*                       ----------------------------------------------
           */
                            /*                       UPDATING TARGET SUPERNODE BY TRIVIAL */
                            /*                       SUPERNODE (WITH ONE COLUMN). */

                            /*                       KLPNT  ...  POINTER TO FIRST NONZERO IN
           * ACTIVE */
                            /*                                   PORTION OF COLUMN FKCOL. */
                            /*                       KDPNT  ...  POINTER TO DIAGONAL ENTRY OF */
                            /*                                   COLUMN FKCOL. */
                            /*                       ----------------------------------------------
           */
                            klpnt = xlnz[fkcol + 1] - klen;
                            kdpnt = xlnz[fkcol];
                            mmpyi_(klen, ncolup, &lindx[kxpnt], &lnz[klpnt], lnz[kdpnt],
                                   &xlnz[1], &lnz[1], &indmap[1]);

                        } else {

                            /*                       --------------------------------------------
           */
                            /*                       KFIRST ...  FIRST INDEX OF ACTIVE PORTION OF
           */
                            /*                                   SUPERNODE KSUP (FIRST COLUMN TO
           */
                            /*                                   BE UPDATED). */
                            /*                       KLAST  ...  LAST INDEX OF ACTIVE PORTION OF
           */
                            /*                                   SUPERNODE KSUP. */
                            /*                       --------------------------------------------
           */

                            kfirst = lindx[kxpnt];
                            klast = lindx[kxpnt + klen - 1];
                            inddif = indmap[kfirst] - indmap[klast];

                            if (inddif < klen) {

                                /*                           ---------------------------------------
             */
                                /*                           DENSE CMOD(JSUP,KSUP). */

                                /*                           ILPNT  ...  POINTER TO FIRST NONZERO IN
             */
                                /*                                       COLUMN KFIRST. */
                                /*                           ILEN   ...  LENGTH OF COLUMN KFIRST. */
                                /*                           ---------------------------------------
             */
                                ilpnt = xlnz[kfirst];
                                ilen = xlnz[kfirst + 1] - ilpnt;
                                mmpy(klen, nkcols, ncolup, &split[fkcol], &xlnz[fkcol], &lnz[1],
                                     &lnz[ilpnt], ilen);

                            } else {

                                /*                           ------------------------------- */
                                /*                           GENERAL SPARSE CMOD(JSUP,KSUP). */
                                /*                           COMPUTE CMOD(JSUP,KSUP) UPDATE */
                                /*                           IN WORK STORAGE. */
                                /*                           ------------------------------- */
                                store = klen * ncolup - ncolup * (ncolup - 1) / 2;
                                if (store > tmpsiz) {
                                    iflag = -2;
                                    return;
                                }
                                mmpy(klen, nkcols, ncolup, &split[fkcol], &xlnz[fkcol], &lnz[1],
                                     &temp[1], klen);
                                /*                           ----------------------------------------
             */
                                /*                           GATHER INDICES OF KSUP RELATIVE TO
             * JSUP. */
                                /*                           ----------------------------------------
             */
                                igathr(klen, &lindx[kxpnt], &indmap[1], &relind[1]);
                                /*                           --------------------------------------
             */
                                /*                           INCORPORATE THE CMOD(JSUP,KSUP) BLOCK
             */
                                /*                           UPDATE INTO THE TO APPROPRIATE COLUMNS
             */
                                /*                           OF L. */
                                /*                           --------------------------------------
             */
                                assmb(klen, ncolup, &temp[1], &relind[1], &xlnz[fjcol], &lnz[1],
                                      jlen);
                            }
                        }

                    } else {

                        /*                   ---------------------------------------------- */
                        /*                   DENSE CMOD(JSUP,KSUP). */
                        /*                   JSUP AND KSUP HAVE IDENTICAL STRUCTURE. */

                        /*                   JLPNT  ...  POINTER TO FIRST NONZERO IN COLUMN */
                        /*                               FJCOL. */
                        /*                   ---------------------------------------------- */
                        jlpnt = xlnz[fjcol];
                        mmpy(klen, nkcols, njcols, &split[fkcol], &xlnz[fkcol], &lnz[1],
                             &lnz[jlpnt], jlen);
                        ncolup = njcols;
                        if (klen > njcols) {
                            nxtcol = lindx[jxpnt + njcols];
                        }
                    }

                    /*               ------------------------------------------------ */
                    /*               LINK KSUP INTO LINKED LIST OF THE NEXT SUPERNODE */
                    /*               IT WILL UPDATE AND DECREMENT KSUP'S ACTIVE */
                    /*               LENGTH. */
                    /*               ------------------------------------------------ */
                    if (klen > ncolup) {
                        nxtsup = snode[nxtcol];
                        link[ksup] = link[nxtsup];
                        link[nxtsup] = ksup;
                        length[ksup] = klen - ncolup;
                    } else {
                        length[ksup] = 0;
                    }

                    /*               ------------------------------- */
                    /*               NEXT UPDATING SUPERNODE (KSUP). */
                    /*               ------------------------------- */
                    ksup = nxksup;
                    goto L300;
                }

                /*           ---------------------------------------------- */
                /*           APPLY PARTIAL CHOLESKY TO THE COLUMNS OF JSUP. */
                /*           ---------------------------------------------- */
                chlsup(jlen, njcols, &split[fjcol], &xlnz[fjcol], &lnz[1], iflag);
                if (iflag != 0) {
                    iflag = -1;
                    return;
                }

                /*           ----------------------------------------------- */
                /*           INSERT JSUP INTO LINKED LIST OF FIRST SUPERNODE */
                /*           IT WILL UPDATE. */
                /*           ----------------------------------------------- */
                if (jlen > njcols) {
                    nxtcol = lindx[jxpnt + njcols];
                    nxtsup = snode[nxtcol];
                    link[jsup] = link[nxtsup];
                    link[nxtsup] = jsup;
                    length[jsup] = jlen - njcols;
                } else {
                    length[jsup] = 0;
                }
            }
        } /* blkfc2 */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFCT .....  BLOCK GENERAL SPARSE LDL'         ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE CALLS THE BLOCK GENERAL SPARSE LDL' ROUTINE, */
/*       BLKFC2. */

/*   INPUT PARAMETERS: */
/*       OUTUNT          -   OUTPUT UNIT. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       NUNROL          -   LOOP UNROLLING LEVEL. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING */
/*                           IT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING */
/*                           THE DIAGONAL ELEMENTS). */
/*       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED. */
/*       IWSIZ           -   SIZE OF INTEGER WORKING STORAGE */
/*       TMPSIZ          -   SIZE OF FLOATING POINT WORKING STORAGE. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY. */
/*       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   ON OUTPUT, CONTAINS TRIANGULAR FACTOR */
/*                           AND DIAGONAL MATRIX OF LDL' DECOMPOSITION. */
/*       DIAG            -   DIAGONAL MATRIX OF LDL' DECOMPOSITION. */
/*       IFLAG           -   ERROR FLAG. */
/*                               0: SUCCESSFUL FACTORIZATION. */
/*                              -1: ZERO DIAGONAL ENCOUNTERED, */
/*                                  MATRIX IS SINGULAR. */
/*                              -2: INSUFFICIENT WORKING STORAGE */
/*                                  [TEMP(*)]. */
/*                              -3: INSUFFICIENT WORKING STORAGE */
/*                                  [IWORK(*)]. */

/*   WORKING PARAMETERS: */
/*       IWORK           -   INTEGER WORKING STORAGE OF LENGTH */
/*                           2*NEQNS + 2*NSUPER. */
/*       TMPVEC          -   COMPLEX*16 WORKING STORAGE OF LENGTH */
/*                           NEQNS. */

/* *********************************************************************** */

        template<typename Scalar, typename Index>
        void blkfct(Index neqns, Index nsuper, Index *xsuper, Index *snode, Index *split,
                    Index *xlindx, Index *lindx, Index *xlnz, Scalar *lnz, Scalar *diag,
                    Index iwsiz, Index *iwork, Index tmpsiz, Scalar *tmpvec,
                    Index &iflag) {
            /* Parameter adjustments */
            --diag;
            --lnz;
            --xlnz;

            /* Function Body */
            iflag = 0;
            if (iwsiz < (neqns << 1) + (nsuper << 1)) {
                iflag = -3;
                std::cerr << "\n";
                std::cerr << "*** INTEGER WORK SPACE = " << iwsiz << "\n";
                std::cerr << "*** IS SMALLER THAN REQUIRED = "
                          << (nsuper << 1) + (neqns << 1) << "\n";
                std::cerr << "\n";
                return;
            }
            iflag = 0;
            if (iwsiz < (neqns << 1) + (nsuper << 1)) {
                iflag = -3;
                return;
            }

            blkfc2(nsuper, xsuper, snode, split, xlindx, lindx,
                   &xlnz[1], &lnz[1], iwork, &iwork[nsuper],
                   &iwork[(nsuper << 1)], &iwork[(nsuper << 1) + neqns], tmpsiz,
                   tmpvec, iflag);

            if (iflag == -1) {
                std::cerr << "\n";
                std::cerr << "*** MATRIX IS SINGULAR ***\n";
                std::cerr << "*** ZERO DIAGONAL ENTRY ENCOUNTERED ***\n";
                std::cerr << "\n";
                return;
            } else if (iflag == -2) {
                std::cerr << "\n";
                std::cerr << "*** INSUFFICIENT WORK STORAGE [TMPVEC(*)] ***\n";
                std::cerr << "\n";
                return;
            }
            if (iflag == 0) {
                for (Index jcol = 1; jcol <= neqns; ++jcol) {
                    diag[jcol] = lnz[xlnz[jcol]];
                }
            }
        } /* blkfct */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP */
/*       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION. */

/*   NOTE: */
/*       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E., */
/*       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES). */

/*   INPUT PARAMETERS: */
/*       NEQNS       -   NUMBER OF EQUATIONS. */
/*       NNZA        -   LENGTH OF ADJACENCY STRUCTURE. */
/*       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                       TO THE ADJACENCY STRUCTURE. */
/*       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                       THE ADJACENCY STRUCTURE. */
/*       PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       POSTORDERING. */
/*       INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       INVERSE OF THE POSTORDERING. */
/*       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                       OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                       INCLUDING THE DIAGONAL ENTRY. */
/*       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                       THE DIAGONAL ENTRIES. */
/*       NSUB        -   NUMBER OF SUBSCRIPTS. */
/*       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                       SUPERNODE MEMBERSHIP. */
/*       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                       SUPERNODE PARTITIONING. */
/*       IFLAG(*)    -   ERROR FLAG. */
/*                          0: SUCCESSFUL SF INITIALIZATION. */
/*                         -1: INSUFFICENT WORKING STORAGE */
/*                             [IWORK(*)]. */

/*   WORK PARAMETERS: */
/*       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3. */

/*   FIRST CREATED ON    NOVEMEBER 14, 1994. */
/*   LAST UPDATED ON     January 12, 1995. */

/* *********************************************************************** */

        template<typename Index>
        void sfinit(Index neqns, Index nnza, Index *xadj, Index *adjncy,
                    Index *perm, Index *invp, Index *colcnt, Index &nnzl, Index &nsub,
                    Index &nsuper, Index *snode, Index *xsuper, Index iwsiz,
                    Index *iwork, Index &iflag) {
            iflag = 0;
            if (iwsiz < neqns * 7 + 3) {
                iflag = -1;
                std::cerr << "\n";
                std::cerr << "*** INTEGER WORK STORAGE = " << iwsiz << "\n";
                std::cerr << "*** IS SMALLER THAN REQUIRED = " << 3 + 7 * (neqns) << "\n";
                std::cerr << "\n";
                return;
            }

            /*       ------------------------------------------ */
            /*       COMPUTE ELIMINATION TREE AND POSTORDERING. */
            /*       ------------------------------------------ */
            etordr(neqns, xadj, adjncy, perm, invp, iwork,
                   &iwork[neqns], &iwork[(neqns << 1)], &iwork[neqns * 3]);

            /*       --------------------------------------------- */
            /*       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS. */
            /*       --------------------------------------------- */
            fcnthn(neqns, xadj, adjncy, perm, invp, iwork,
                   snode, colcnt, nnzl, &iwork[neqns],
                   &iwork[(neqns << 1)], xsuper, &iwork[neqns * 3],
                   &iwork[(neqns << 2) + 1], &iwork[neqns * 5 + 2],
                   &iwork[neqns * 6 + 3]);

            /*       --------------------------------------------------------- */
            /*       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM */
            /*       NUMBER OF NONZEROS IN ITS COLUMN OF L. */
            /*       --------------------------------------------------------- */
            chordr(neqns, xadj, adjncy, perm, invp, colcnt, iwork,
                   &iwork[neqns], &iwork[(neqns << 1)], &iwork[neqns * 3]);

            /*       ---------------- */
            /*       FIND SUPERNODES. */
            /*       ---------------- */
            fsup1(neqns, iwork, colcnt, nsub, nsuper, snode);
            fsup2(neqns, nsuper, snode, xsuper);
        } /* sfinit */

    } // namespace details

//
// Definition of templated functions for class SymmetricSparse
//

    template<typename Scalar>
    SymmetricSparse<Scalar>::SymmetricSparse(int n_,
                                             const int *colptr_, const int *rowind_,
                                             int &Lnnz_,
                                             const int order_, const int *perm_) {

        int nnz = 0, nnza = 0, ibegin = 0, iend = 0, i = 0, j = 0, iwsiz = 0,
                iflag = 0;
        bool sfiflg = true, notfound = true;

        int pjend = 0, pibeg = 0, irow = 0, jcol = 0;

#ifdef TIMING
        double t0, t1;
#endif

        n = n_;
        nnz = colptr_[n] - 1;

        /* Check to see whether a full representation the symmetric
     matrix is used */
        j = 1;
        while (notfound && j < n) {
            /* pointer to the last entry of column j */
            pjend = colptr_[j] - 1;
            if (pjend >= colptr_[j - 1]) {
                /* there is at least one off-diagonal entry,
         find its row index */
                irow = rowind_[pjend - 1];

                /* get column pointer to the first entry of irow-th column */
                pibeg = colptr_[irow - 1];

                /* check to see if upper triangular part is present */
                if (rowind_[pibeg - 1] != irow) {
                    fullrep = true;
                }
                notfound = false;
            }
            j++;
        }

        /* nnza = number of nonzeros off diagonal */
        nnza = (fullrep) ? nnz - n : 2 * (nnz - n);
        iwsiz = 7 * n + 3;

        /* --------------------------------------------------
        Allocate matrix storage required for reordering
        and symbolic factorization.
     -------------------------------------------------- */

        /* adj, xadj contain the structure of the
     full representation of the original matrix */

        xadj.resize(n + 1);
        adj.resize(nnza + n);

        /* ----------------------------------------------------------
     convert the indices and pointers of the lower triangular
     representation of a symmetric HB matrix to its full
     representation.
     ---------------------------------------------------------- */

        iwork.resize(iwsiz);

        if (fullrep) {
            /* make a copy of the non-zero structure, take out the
       diagonals */
            xadj[0] = 1;
            for (i = 0; i < n; i++) {
                xadj[i + 1] = xadj[i] + (colptr_[i + 1] - colptr_[i] - 1);
            }
            if (xadj[n] - 1 != nnza) {
                std::cerr << " Something is wrong with the matrix!\n";
            }

            i = 0;
            for (jcol = 1; jcol <= n; jcol++) {
                ibegin = colptr_[jcol - 1];
                iend = colptr_[jcol] - 1;
                for (irow = ibegin; irow <= iend; irow++) {
                    if (rowind_[irow - 1] != jcol) {
                        adj[i] = rowind_[irow - 1];
                        i++;
                    }
                }
            }
        } else {
            /* convert */
            details::ilo2ho<int>(n, nnza, colptr_, rowind_, &xadj[0], &adj[0], &iwork[0]);
        }

        /* -------------------------------------------------
     Make a copy of the converted matrix for reordering
     ------------------------------------------------- */

        /* Since reordering will destroy the original matrix
     we need to make an extra copy of the indices and pointers
     of the full representation of the original matrix.  */

        std::vector<int> adj2(adj), xadj2(xadj);

        /* ----------------------------------------
     Allocate storage for supernode partition
     ----------------------------------------*/

        factor = std::unique_ptr<LDLt_factor < Scalar> > (new LDLt_factor<Scalar>());
        factor->snodes.resize(n);
        factor->xsuper.resize(n + 1);
        factor->colcnt.resize(n);
        factor->perm.resize(n);
        factor->invp.resize(n);

#ifdef TIMING
        t0 = getime();
#endif

#ifdef METIS
        if (order_ < 0 || order_ > 3) {
#else
        if (order_ < 0 || order_ > 1) {
#endif

            /* ------------------------
        Multiple Minimum Degree
       ------------------------ */

            details::ordmmd<int>(n, &xadj2[0], &adj2[0], &factor->invp[0],
                                 &factor->perm[0], iwsiz, &iwork[0], factor->nnzl,
                                 factor->nsub, &factor->colcnt[0], factor->nsuper,
                                 &factor->xsuper[0], &factor->snodes[0], sfiflg, iflag);

        } else if (order_ == 1) {
            printf(" AMD not implemented yet, use MMD !\n");

            details::ordmmd<int>(n, &xadj2[0], &adj2[0], &factor->invp[0],
                                 &factor->perm[0], iwsiz, &iwork[0], factor->nnzl,
                                 factor->nsub, &factor->colcnt[0], factor->nsuper,
                                 &factor->xsuper[0], &factor->snodes[0], sfiflg, iflag);
        }
#ifdef METIS
                                                                                                                                    else if (order_ == 2) {

    /* ----------------------
       Node Nested Dissection
       ---------------------- */

    int numflag = 1, options[8];
    options[0] = 0;

    /* use default options now
      options[1] = 3;
      options[2] = 1;
      options[3] = 2;
      options[4] = 0;
      options[5] = 3;
      options[6] = 0;
      options[7] = 1;
    */

    METIS_NodeND(&n, &xadj2[0], &adj2[0], &numflag, options, &factor->perm[0],
                 &factor->invp[0]);
  } else if (order_ == 3) {
    /* ----------------------
       Edge Nested Dissection
       ---------------------- */
    int numflag = 1, options[8];
    options[0] = 0;
    METIS_EdgeND(&n, &xadj2[0], &adj2[0], &numflag, options, &factor->perm[0],
                 &factor->invp[0]);
  }
#endif
        else if (order_ == 0) {
            /* use the input permutation vector */
            for (i = 0; i < n; i++) {
                factor->perm[i] = perm_[i];
                factor->invp[factor->perm[i] - 1] = i + 1;
            }
        }
        xadj2.clear();
        adj2.clear();

#ifdef TIMING
                                                                                                                                t1 = getime();
  t1 = t1 - t0;
  printf("\n");

  if (order_ < 0 || order_ > 4) {
    printf(" MMD   :");
  } else if (*order == 1)
    printf(" AMD         :");
  else if (*order == 2) {
    printf(" NodeND      :");
  } else if (*order == 3) {
    printf(" EdgeND      :");
  } else if (*order == 0) {
    printf(" User defined:");
  }
  printf(" Time reordering    = %9.3e\n", t1);
#endif

#ifdef TIMING
        t0 = getime();
#endif
        /* ----------------------
     Symbolic factorization
     ---------------------- */
        if (order_ >= 0 && order_ <= 3) {
            /* not needed when MMD has been called */
            for (i = 0; i < iwsiz; i++)
                iwork[i] = 0;
            details::sfinit(n, nnza, &xadj[0], &adj[0],
                            &factor->perm[0], &factor->invp[0], &factor->colcnt[0],
                            factor->nnzl, factor->nsub, factor->nsuper,
                            &factor->snodes[0], &factor->xsuper[0], iwsiz, &iwork[0],
                            iflag);
        }

        Lnnz_ = factor->nnzl;

        factor->xlindx.resize(n + 1);
        factor->lindx.resize(factor->nsub);
        factor->xlnz.resize(n + 1);

        for (i = 0; i < iwsiz; i++)
            iwork[i] = 0;

        details::symfct(n, nnza, &xadj[0], &adj[0],
                        &factor->perm[0], &factor->invp[0], &factor->colcnt[0],
                        factor->nsuper, &factor->xsuper[0], &factor->snodes[0],
                        factor->nsub, &factor->xlindx[0], &factor->lindx[0],
                        &factor->xlnz[0], iwsiz, &iwork[0], iflag);

#ifdef TIMING
                                                                                                                                t1 = getime();
  t1 = t1 - t0;
  printf(" Time Symbolic factorization      = %9.3e\n\n", t1);
#endif

        /* ---------------------------
     Prepare for Numerical factorization
     and triangular solve.
     --------------------------- */

        factor->split.resize(n);

        int cachsz_ = 700;
        details::bfinit(n, factor->nsuper, &factor->xsuper[0], &factor->snodes[0],
                        &factor->xlindx[0], &factor->lindx[0], cachsz_, factor->tmpsiz,
                        &factor->split[0]);

    }


/// \brief CSC matrix
/// colptr_ Array of size 'neqns + 1' for pointing entries of column
/// rowind_ Array of size 'nnz + 1' for indicating row indices
/// nzvals_
    template<typename Scalar>
    template<typename Index>
    void SymmetricSparse<Scalar>::ldlTFactorize(const Index *colptr_,
                                                const Index *rowind_,
                                                const Scalar *nzvals_) {
        // Assume that rowind[*] is 1-based
        // Same for colptr
        Index i, nnz, nnzl, neqns, nsuper, tmpsiz, iflag, iwsiz;
#ifdef TIMING
        double t0, t1;
#endif
        Index maxsup, jsup, nnzlplus, supsize;
        auto &myMat = *(factor.get());

        neqns = n;
        nnz = colptr_[neqns] - 1;
        nnzl = myMat.nnzl;
        nsuper = myMat.nsuper;
        tmpsiz = myMat.tmpsiz;
        iwsiz = 7 * neqns + 3;

#ifdef TIMING
        t0 = getime();
#endif
        anz.resize(2 * nnz - neqns);

        /* Aug 30, 2008, allocate extra space for diagonal extraction
     the extra space is not used in the LDLT factorization */
        maxsup = 0;
        nnzlplus = nnzl;
        for (jsup = 0; jsup < nsuper; jsup++) {
            supsize = myMat.xsuper[jsup + 1] - myMat.xsuper[jsup];
            maxsup = std::max(maxsup, supsize);
            nnzlplus = nnzlplus + supsize * (supsize - 1) / 2;
        }
        printf(" nnzlplus = %d\n", nnzlplus);

        myMat.lnz.resize(nnzlplus);
        myMat.tmat.resize(tmpsiz);
        myMat.diag.resize(neqns);
        myMat.newrhs.resize(neqns);

        if (fullrep) {
            for (i = 0; i < neqns + 1; i++)
                xadj[i] = colptr_[i];
            for (i = 0; i < nnz; i++)
                adj[i] = rowind_[i];
            for (i = 0; i < nnz; i++)
                anz[i] = nzvals_[i];
        } else {
            details::flo2ho(neqns, colptr_, rowind_, nzvals_, &xadj[0], &adj[0],
                            &anz[0], &iwork[0]);
        }

        details::inpnv(neqns, &xadj[0], &adj[0], &anz[0],
                        &myMat.perm[0], &myMat.invp[0], myMat.nsuper,
                        &myMat.xsuper[0], &myMat.xlindx[0], &myMat.lindx[0],
                        &myMat.xlnz[0], &myMat.lnz[0], iwsiz, &iwork[0], iflag);

        myMat.tmat.assign(tmpsiz, 0);

        details::blkfct(neqns, nsuper, &myMat.xsuper[0], &myMat.snodes[0],
                        &myMat.split[0], &myMat.xlindx[0], &myMat.lindx[0],
                        &myMat.xlnz[0], &myMat.lnz[0], &myMat.diag[0], iwsiz,
                        &iwork[0], tmpsiz, &myMat.tmat[0], iflag);
#ifdef TIMING
                                                                                                                                t1 = getime();
  t1 = t1 - t0;
  printf(" Time numerical factorization      = %9.3e\n\n", t1);
#endif

    }

    template<typename Scalar>
    void SymmetricSparse<Scalar>::solve(Scalar *x, const Scalar *rhs) {

        if (factor == nullptr) {
            std::cerr << "\n !! Error -- Needs to call factorization first !! \n\n";
            exit(-1);
        }

        int nsuper = factor->nsuper;
        const auto perm = factor->perm;
        const auto invp = factor->invp;
        auto &newrhs = factor->newrhs;

        for (int i = 0; i < n; i++)
            newrhs[i] = rhs[perm[i] - 1];

        details::blkslv_(nsuper, &factor->xsuper[0], &factor->xlindx[0], &factor->lindx[0],
                         &factor->xlnz[0], &factor->lnz[0], &factor->newrhs[0]);

        for (int i = 0; i < n; i++)
            x[i] = newrhs[invp[i] - 1];
    }

} // namespace NgPeytonCpp

#endif // CPPSELINV_SYMMETRICSPARSE_IMPL_H
