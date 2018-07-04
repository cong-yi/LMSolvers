#if defined _WIN32 || defined __CYGWIN__
#ifdef DLL_EXPORT
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllexport))
#else
#define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#endif
#else
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllimport))
#else
#define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
#endif
#endif
#define DLL_LOCAL
#else
#if __GNUC__ >= 4
#define DLL_PUBLIC __attribute__ ((visibility ("default")))
#define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif
#endif

#define SPLM_OPTS_SZ        6 /* max(5, 6) */
#define SPLM_INFO_SZ        10

namespace sparselm
{
	/* sparse LM for functions with CRS and CCS Jacobians */

	DLL_PUBLIC int sparselm_dercrs(
		void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
		void(*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
		double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
		int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

	DLL_PUBLIC int sparselm_derccs(
		void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
		void(*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
		double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
		int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

	DLL_PUBLIC int sparselm_difcrs(
		void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
		void(*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
		double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
		int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

	DLL_PUBLIC int sparselm_difccs(
		void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
		void(*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
		double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
		int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

	/* error checking for CRS and CCS Jacobians */
	DLL_PUBLIC void sparselm_chkjaccrs(
		void(*func)(double *p, double *hx, int m, int n, void *adata),
		void(*jacf)(double *p, struct splm_crsm *jac, int m, int n, void *adata),
		double *p, int m, int n, int jnnz, void *adata, double *err);

	DLL_PUBLIC void sparselm_chkjacccs(
		void(*func)(double *p, double *hx, int m, int n, void *adata),
		void(*jacf)(double *p, struct splm_ccsm *jac, int m, int n, void *adata),
		double *p, int m, int n, int jnnz, void *adata, double *err);

	/* CRS sparse matrices manipulation routines */
	DLL_PUBLIC void splm_crsm_alloc(struct splm_crsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_crsm_alloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_crsm_alloc_values(struct splm_crsm *sm);
	DLL_PUBLIC void splm_crsm_realloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_crsm_free(struct splm_crsm *sm);
	DLL_PUBLIC int splm_crsm_elmidx(struct splm_crsm *sm, int i, int j);
	DLL_PUBLIC int splm_crsm_elmrow(struct splm_crsm *sm, int idx);
	DLL_PUBLIC int splm_crsm_row_elmidxs(struct splm_crsm *sm, int i, int *vidxs, int *jidxs);
	DLL_PUBLIC int splm_crsm_row_maxnelms(struct splm_crsm *sm);
	DLL_PUBLIC int splm_crsm_col_elmidxs(struct splm_crsm *sm, int j, int *vidxs, int *iidxs);
	DLL_PUBLIC void splm_crsm2ccsm(struct splm_crsm *crs, struct splm_ccsm *ccs);
	DLL_PUBLIC void splm_crsm_row_sort(struct splm_crsm *sm);

	/* CCS sparse matrices manipulation routines */
	DLL_PUBLIC void splm_ccsm_alloc(struct splm_ccsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_ccsm_alloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_ccsm_alloc_values(struct splm_ccsm *sm);
	DLL_PUBLIC void splm_ccsm_realloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz);
	DLL_PUBLIC void splm_ccsm_free(struct splm_ccsm *sm);
	DLL_PUBLIC int splm_ccsm_elmidx(struct splm_ccsm *sm, int i, int j);
	DLL_PUBLIC int splm_crsm_elmcol(struct splm_ccsm *sm, int idx);
	DLL_PUBLIC int splm_ccsm_row_elmidxs(struct splm_ccsm *sm, int i, int *vidxs, int *jidxs);
	DLL_PUBLIC int splm_ccsm_col_elmidxs(struct splm_ccsm *sm, int j, int *vidxs, int *iidxs);
	DLL_PUBLIC int splm_ccsm_col_maxnelms(struct splm_ccsm *sm);
	DLL_PUBLIC void splm_ccsm2crsm(struct splm_ccsm *ccs, struct splm_crsm *crs);
	DLL_PUBLIC int splm_ccsm_drop_cols(struct splm_ccsm *A, int ncols);
	DLL_PUBLIC void splm_ccsm_restore_cols(struct splm_ccsm *A, int ncols, int ncnnz);
	DLL_PUBLIC void splm_ccsm_col_sort(struct splm_ccsm *sm);

	DLL_PUBLIC double splm_gettime(void);


	/* Sparse matrix representation using Sparse Triplet (ST) format.
	* Note that the matrix might have an allocated size (maxnnz) larger
	* than the number of elements it actually holds (nnz).
	* Primarily intended to be used for setting up the structure
	* of a corresponding CCS matrix.
	* See http://people.sc.fsu.edu/~jburkardt/data/st/st.html
	*/

	DLL_PUBLIC void splm_stm_alloc(struct splm_stm *sm, int nr, int nc, int maxnnz);
	DLL_PUBLIC void splm_stm_allocval(struct splm_stm *sm, int nr, int nc, int maxnnz);
	DLL_PUBLIC void splm_stm_free(struct splm_stm *sm);
	DLL_PUBLIC int splm_stm_nonzero(struct splm_stm *sm, int i, int j);
	DLL_PUBLIC int splm_stm_nonzeroval(struct splm_stm *sm, int i, int j, double val);
	DLL_PUBLIC void splm_stm2ccsm(struct splm_stm *st, struct splm_ccsm *ccs);
	DLL_PUBLIC void splm_tri2ccsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_ccsm *ccs);
	DLL_PUBLIC void splm_stm2crsm(struct splm_stm *st, struct splm_crsm *crs);
	DLL_PUBLIC void splm_tri2crsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_crsm *crs);

}