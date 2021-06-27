/*************************************************************************
	> File Name: local_alignment.c
	> Author: 
	> Mail: 
	> Created Time: Sun May  3 15:20:52 2020
 ************************************************************************/

#include<stdio.h>
#include "ksw2.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>


unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline void update_max_zdrop(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int *max_zdrop)
{
	if (score < *max) {
		int li = i - *max_i;
		int lj = j - *max_j;
		int diff = li > lj? li - lj : lj - li;
		int z = *max - score - diff * e;
		if (z > *max_zdrop) {
			*max_zdrop = z;
		}
	} else *max = score, *max_i = i, *max_j = j;
}


static inline int test_zdrop(const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat, const int8_t q, const int8_t e)
{
    uint32_t k;
	int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
	int q_len, t_len;

	// find the score and the region where score drops most along diagonal
	for (k = 0, score = 0; k < n_cigar; ++k) {
		uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 0) {
			for (l = 0; l < len; ++l) {
				score += mat[tseq[i + l] * 5 + qseq[j + l]];
				update_max_zdrop(score, i+l, j+l, &max, &max_i, &max_j, e, &max_zdrop);
			}
			i += len, j += len;
		} else if (op == 1 || op == 2 || op == 3) {
			score -= q + e * len;
			if (op == 1) j += len; // insertion
			else i += len;         // deletion
			update_max_zdrop(score, i, j, &max, &max_i, &max_j, e, &max_zdrop);
		}
	}

	return max_zdrop;
}

static PyObject* ksw2_aligner(PyObject *self, PyObject *args)
{
    char *query;
    char *target;
    int x_score;
    if(! PyArg_ParseTuple(args, "ssi", &query, &target, &x_score))
        return NULL;

    int i;
    //int a = 2, b = -4;
    //int gapo = 4, gapo2 = 44, gape = 2, gape2 = 1;
    int a = 1, b = -x_score;
    int gapo = 4, gapo2 = 44, gape = 2, gape2 = 1;  //gapo=8
    
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };

    int lenq = strlen(query);
    int lent = strlen(target);
    uint8_t *qs = (uint8_t *)malloc(lenq + 10);
    uint8_t *ts = (uint8_t *)malloc(lent + 10);
    for(i = 0; i < lenq; ++i) qs[i] = seq_nt4_table[(int)query[i]];
    for(i = 0; i < lent; ++i) ts[i] = seq_nt4_table[(int)target[i]];

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    int FLAG = 0 | KSW_EZ_APPROX_MAX;
    //int FLAG = 0 | KSW_EZ_EXTZ_ONLY;
    //int FLAG = 0 |KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR;
    //int FLAG = 0 |KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT;

    //ksw_extd2_sse(0, lenq, qs, lent, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, -1, FLAG, &ez);
    ksw_extz2_sse(0, lenq, qs, lent, ts, 5, mat, gapo, gape, -1, -1, -1, FLAG, &ez);
    
    int droped = test_zdrop(qs, ts, ez.n_cigar, ez.cigar, mat, gapo, gape);

    char *cr;
    int offset = 0;
    cr = (char*)malloc(ez.n_cigar<<1 + 2);
    for(i = 0;i < ez.n_cigar; ++i)
    {
        offset += sprintf(cr+offset, "%d%c", ez.cigar[i]>>4, "MIDNS"[ez.cigar[i]&0xf]);
    }

    if(qs != NULL) free(qs);
    if(ts != NULL) free(ts);
    //PyObject *rList = Py_BuildValue("[i,s]", ez.score, cr);
    PyObject *rList = Py_BuildValue("[i,s]", droped, cr);

    if(ez.cigar != NULL) free(ez.cigar);
    free(cr);
    
    return rList;
}

/*define functions in module */
static PyMethodDef myMethods[] = {
    {"ksw2_aligner", ksw2_aligner, METH_VARARGS, "Execute other function."},
    {NULL, NULL, 0, NULL}
};

//#if PY_MAJOR_VERSION >= 3
/*module initialization */
/*Python version 3 */
static struct PyModuleDef ksw2ModuleDem = 
{
    PyModuleDef_HEAD_INIT,
    "ksw2_module", "Some documentation",
    -1,
    myMethods
};
PyMODINIT_FUNC PyInit_ksw2_module(void)
{
    return PyModule_Create(&ksw2ModuleDem);
}

//void gap_align(const uint8_t *ts, const int tl, const uint8_t *qs, const int ql, int sc_mch, int sc_mis, int gapo, int gapo2, int gape, int gape2, int flag, ksw_extz_t *ez)
//{
//	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
//	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
//
//
//    int FLAG = 0;
//    if (flag == -1) //left ext
//        FLAG = 0|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT;
//    else if(flag == 0) // end-to-end
//        FLAG = 0|KSW_EZ_APPROX_MAX;
//    else
//        FLAG = 0|KSW_EZ_EXTZ_ONLY;
//    
//    //ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, -1, FLAG, ez);
//    
//    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, -1, FLAG, ez);
//}
//
