/*
C->D glpk header interface file is a direct conversion from the C header file
Authors: Chibisi Chima-Okereke
*/

module glpk;
import core.stdc.stdarg;

extern (C) {

/* library version numbers: */
enum GLP_MAJOR_VERSION = 4;
enum GLP_MINOR_VERSION = 57;

struct glp_prob{
}
/* LP/MIP problem object */

/* optimization direction flag: */
enum GLP_MIN = 1;  /* minimization */
enum GLP_MAX = 2;  /* maximization */

/* kind of structural variable: */
enum GLP_CV = 1;  /* continuous variable */
enum GLP_IV = 2;  /* integer variable */
enum GLP_BV = 3;  /* binary variable */

/* type of auxiliary/structural variable: */
enum GLP_FR = 1;  /* free (unbounded) variable */
enum GLP_LO = 2;  /* variable with lower bound */
enum GLP_UP = 3;  /* variable with upper bound */
enum GLP_DB = 4;  /* double-bounded variable */
enum GLP_FX = 5;  /* fixed variable */

/* status of auxiliary/structural variable: */
enum GLP_BS = 1;  /* basic variable */
enum GLP_NL = 2;  /* non-basic variable on lower bound */
enum GLP_NU = 3;  /* non-basic variable on upper bound */
enum GLP_NF = 4;  /* non-basic free (unbounded) variable */
enum GLP_NS = 5;  /* non-basic fixed variable */

/* scaling options: */
enum GLP_SF_GM = 0x01;  /* perform geometric mean scaling */
enum GLP_SF_EQ = 0x10;  /* perform equilibration scaling */
enum GLP_SF_2N  = 0x20;  /* round scale factors to power of two */
enum GLP_SF_SKIP = 0x40;  /* skip if problem is well scaled */
enum GLP_SF_AUTO = 0x80;  /* choose scaling options automatically */

/* solution indicator: */
enum GLP_SOL =  1;  /* basic solution */
enum GLP_IPT =  2;  /* interior-point solution */
enum GLP_MIP =  3;  /* mixed integer solution */

/* solution status: */
enum GLP_UNDEF = 1;  /* solution is undefined */
enum GLP_FEAS  = 2;  /* solution is feasible */
enum GLP_INFEAS = 3;  /* solution is infeasible */
enum GLP_NOFEAS = 4;  /* no feasible solution exists */
enum GLP_OPT  = 5;  /* solution is optimal */
enum GLP_UNBND = 6;  /* solution is unbounded */

struct glp_bfcp
{     /* basis factorization control parameters */
      int msg_lev;            /* (not used) */
      int type;               /* factorization type: */
      enum GLP_BF_LUF = 0x00;  /* plain LU-factorization */
      enum GLP_BF_BTF = 0x10;  /* block triangular LU-factorization */
      enum GLP_BF_FT = 0x01;  /* Forrest-Tomlin (LUF only) */
      enum GLP_BF_BG = 0x02;  /* Schur compl. + Bartels-Golub */
      enum GLP_BF_GR = 0x03;  /* Schur compl. + Givens rotation */
      int lu_size;            /* (not used) */
      double piv_tol;         /* sgf_piv_tol */
      int piv_lim;            /* sgf_piv_lim */
      int suhl;               /* sgf_suhl */
      double eps_tol;         /* sgf_eps_tol */
      double max_gro;         /* (not used) */
      int nfs_max;            /* fhvint.nfs_max */
      double upd_tol;         /* (not used) */
      int nrs_max;            /* scfint.nn_max */
      int rs_size;            /* (not used) */
      double foo_bar[38];     /* (reserved) */
}

struct glp_smcp
{     /* simplex method control parameters */
      int msg_lev;            /* message level: */
      enum GLP_MSG_OFF = 0;  /* no output */
      enum GLP_MSG_ERR = 1;  /* warning and error messages only */
      enum GLP_MSG_ON = 2;  /* normal output */
      enum GLP_MSG_ALL = 3;  /* full output */
      enum GLP_MSG_DBG = 4;  /* debug output */
      int meth;               /* simplex method option: */
      enum GLP_PRIMAL = 1;  /* use primal simplex */
      enum GLP_DUALP = 2;  /* use dual; if it fails, use primal */
      enum GLP_DUAL = 3;  /* use dual simplex */
      int pricing;            /* pricing technique: */
      enum GLP_PT_STD = 0x11;  /* standard (Dantzig's rule) */
      enum GLP_PT_PSE = 0x22;  /* projected steepest edge */
      int r_test;             /* ratio test technique: */
      enum GLP_RT_STD = 0x11;  /* standard (textbook) */
      enum GLP_RT_HAR = 0x22;  /* Harris' two-pass ratio test */
      double tol_bnd;         /* spx.tol_bnd */
      double tol_dj;          /* spx.tol_dj */
      double tol_piv;         /* spx.tol_piv */
      double obj_ll;          /* spx.obj_ll */
      double obj_ul;          /* spx.obj_ul */
      int it_lim;             /* spx.it_lim */
      int tm_lim;             /* spx.tm_lim (milliseconds) */
      int out_frq;            /* spx.out_frq */
      int out_dly;            /* spx.out_dly (milliseconds) */
      int presolve;           /* enable/disable using LP presolver */
      double foo_bar[36];     /* (reserved) */
}

struct glp_iptcp
{     /* interior-point solver control parameters */
      int msg_lev;            /* message level (see glp_smcp) */
      int ord_alg;            /* ordering algorithm: */
      enum GLP_ORD_NONE = 0;  /* natural (original) ordering */
      enum GLP_ORD_QMD = 1;  /* quotient minimum degree (QMD) */
      enum GLP_ORD_AMD = 2;  /* approx. minimum degree (AMD) */
      enum GLP_ORD_SYMAMD = 3;  /* approx. minimum degree (SYMAMD) */
      double foo_bar[48];     /* (reserved) */
}

struct glp_tree;
/* branch-and-bound tree */

struct glp_iocp
{     /* integer optimizer control parameters */
      int msg_lev;            /* message level (see glp_smcp) */
      int br_tech;            /* branching technique: */
      enum GLP_BR_FFV = 1;  /* first fractional variable */
      enum GLP_BR_LFV = 2;  /* last fractional variable */
      enum GLP_BR_MFV = 3;  /* most fractional variable */
      enum GLP_BR_DTH = 4;  /* heuristic by Driebeck and Tomlin */
      enum GLP_BR_PCH = 5;  /* hybrid pseudocost heuristic */
      int bt_tech;            /* backtracking technique: */
      enum GLP_BT_DFS = 1;  /* depth first search */
      enum GLP_BT_BFS = 2;  /* breadth first search */
      enum GLP_BT_BLB = 3;  /* best local bound */
      enum GLP_BT_BPH = 4;  /* best projection heuristic */
      double tol_int;         /* mip.tol_int */
      double tol_obj;         /* mip.tol_obj */
      int tm_lim;             /* mip.tm_lim (milliseconds) */
      int out_frq;            /* mip.out_frq (milliseconds) */
      int out_dly;            /* mip.out_dly (milliseconds) */
      alias cb_func = void function(glp_tree *T, void *info); /* mip.cb_func */
      void *cb_info;          /* mip.cb_info */
      int cb_size;            /* mip.cb_size */
      int pp_tech;            /* preprocessing technique: */
      enum GLP_PP_NONE = 0;  /* disable preprocessing */
      enum GLP_PP_ROOT = 1;  /* preprocessing only on root level */
      enum GLP_PP_ALL = 2;  /* preprocessing on all levels */
      double mip_gap;         /* relative MIP gap tolerance */
      int mir_cuts;           /* MIR cuts       (GLP_ON/GLP_OFF) */
      int gmi_cuts;           /* Gomory's cuts  (GLP_ON/GLP_OFF) */
      int cov_cuts;           /* cover cuts     (GLP_ON/GLP_OFF) */
      int clq_cuts;           /* clique cuts    (GLP_ON/GLP_OFF) */
      int presolve;           /* enable/disable using MIP presolver */
      int binarize;           /* try to binarize integer variables */
      int fp_heur;            /* feasibility pump heuristic */
      int ps_heur;            /* proximity search heuristic */
      int ps_tm_lim;          /* proxy time limit, milliseconds */
      int sr_heur;            /* simple rounding heuristic */
      int use_sol;            /* use existing solution */
      const(char)* save_sol;   /* filename to save every new solution */
      int alien;              /* use alien solver */
      double foo_bar[24];     /* (reserved) */
}

struct glp_attr
{     /* additional row attributes */
      int level;
      /* subproblem level at which the row was added */
      int origin;
      /* row origin flag: */
      enum GLP_RF_REG = 0;  /* regular constraint */
      enum GLP_RF_LAZY = 1;  /* "lazy" constraint */
      enum GLP_RF_CUT = 2;  /* cutting plane constraint */
      int klass;  /* row class descriptor: */
      enum GLP_RF_GMI = 1;  /* Gomory's mixed integer cut */
      enum GLP_RF_MIR = 2;  /* mixed integer rounding cut */
      enum GLP_RF_COV = 3;  /* mixed cover cut */
      enum GLP_RF_CLQ = 4;  /* clique cut */
      double foo_bar[7];
      /* (reserved) */
}

/* enable/disable flag: */
enum GLP_ON = 1;  /* enable something */
enum GLP_OFF = 0;  /* disable something */

/* reason codes: */
enum GLP_IROWGEN = 0x01;  /* request for row generation */
enum GLP_IBINGO = 0x02;  /* better integer solution found */
enum GLP_IHEUR = 0x03;  /* request for heuristic solution */
enum GLP_ICUTGEN = 0x04;  /* request for cut generation */
enum GLP_IBRANCH = 0x05;  /* request for branching */
enum GLP_ISELECT = 0x06;  /* request for subproblem selection */
enum GLP_IPREPRO = 0x07;  /* request for preprocessing */

/* branch selection indicator: */
enum GLP_NO_BRNCH = 0;  /* select no branch */
enum GLP_DN_BRNCH = 1;  /* select down-branch */
enum GLP_UP_BRNCH = 2;  /* select up-branch */

/* return codes: */
enum GLP_EBADB = 0x01;  /* invalid basis */
enum GLP_ESING = 0x02;  /* singular matrix */
enum GLP_ECOND = 0x03;  /* ill-conditioned matrix */
enum GLP_EBOUND = 0x04;  /* invalid bounds */
enum GLP_EFAIL = 0x05;  /* solver failed */
enum GLP_EOBJLL = 0x06;  /* objective lower limit reached */
enum GLP_EOBJUL = 0x07;  /* objective upper limit reached */
enum GLP_EITLIM = 0x08;  /* iteration limit exceeded */
enum GLP_ETMLIM = 0x09;  /* time limit exceeded */
enum GLP_ENOPFS = 0x0A;  /* no primal feasible solution */
enum GLP_ENODFS = 0x0B;  /* no dual feasible solution */
enum GLP_EROOT = 0x0C;  /* root LP optimum not provided */
enum GLP_ESTOP = 0x0D;  /* search terminated by application */
enum GLP_EMIPGAP = 0x0E;  /* relative mip gap tolerance reached */
enum GLP_ENOFEAS = 0x0F;  /* no primal/dual feasible solution */
enum GLP_ENOCVG = 0x10;  /* no convergence */
enum GLP_EINSTAB = 0x11;  /* numerical instability */
enum GLP_EDATA = 0x12;  /* invalid data */
enum GLP_ERANGE = 0x13;  /* result out of range */

/* condition indicator: */
enum GLP_KKT_PE = 1;  /* primal equalities */
enum GLP_KKT_PB = 2;  /* primal bounds */
enum GLP_KKT_DE = 3;  /* dual equalities */
enum GLP_KKT_DB = 4;  /* dual bounds */
enum GLP_KKT_CS = 5;  /* complementary slackness */

/* MPS file format: */
enum GLP_MPS_DECK = 1;  /* fixed (ancient) */
enum GLP_MPS_FILE = 2;  /* free (modern) */

struct glp_mpscp
{     /* MPS format control parameters */
      int blank;
      /* character code to replace blanks in symbolic names */
      char *obj_name;
      /* objective row name */
      double tol_mps;
      /* zero tolerance for MPS data */
      double foo_bar[17];
      /* (reserved for use in the future) */
}

struct glp_cpxcp
{     /* CPLEX LP format control parameters */
      double foo_bar[20];
      /* (reserved for use in the future) */
}

struct glp_tran;
/* MathProg translator workspace */

glp_prob* glp_create_prob();
/* create problem object */

void glp_set_prob_name(glp_prob *P, const(char)* name);
/* assign (change) problem name */

void glp_set_obj_name(glp_prob *P, const(char)* name);
/* assign (change) objective function name */

void glp_set_obj_dir(glp_prob *P, int dir);
/* set (change) optimization direction flag */

int glp_add_rows(glp_prob *P, int nrs);
/* add new rows to problem object */

int glp_add_cols(glp_prob *P, int ncs);
/* add new columns to problem object */

void glp_set_row_name(glp_prob *P, int i, const(char)* name);
/* assign (change) row name */

void glp_set_col_name(glp_prob *P, int j, const(char)* name);
/* assign (change) column name */

void glp_set_row_bnds(glp_prob *P, int i, int type, double lb,
      double ub);
/* set (change) row bounds */

void glp_set_col_bnds(glp_prob *P, int j, int type, double lb,
      double ub);
/* set (change) column bounds */

void glp_set_obj_coef(glp_prob *P, int j, double coef);
/* set (change) obj. coefficient or constant term */

void glp_set_mat_row(glp_prob *P, int i, int len, const(int)* ind,
      const(double)* val);
/* set (replace) row of the constraint matrix */

void glp_set_mat_col(glp_prob *P, int j, int len, const(int)* ind,
      const(double)* val);
/* set (replace) column of the constraint matrix */

void glp_load_matrix(glp_prob *P, int ne, const(int)* ia, const(int) * ja, const(double)* ar);
/* load (replace) the whole constraint matrix */

int glp_check_dup(int m, int n, int ne, const(int)* ia, const(int)* ja);
/* check for duplicate elements in sparse matrix */

void glp_sort_matrix(glp_prob *P);
/* sort elements of the constraint matrix */

void glp_del_rows(glp_prob *P, int nrs, const(int)* num);
/* delete specified rows from problem object */

void glp_del_cols(glp_prob *P, int ncs, const(int)* num);
/* delete specified columns from problem object */

void glp_copy_prob(glp_prob *dest, glp_prob *prob, int names);
/* copy problem object content */

void glp_erase_prob(glp_prob *P);
/* erase problem object content */

void glp_delete_prob(glp_prob *P);
/* delete problem object */

const(char)* glp_get_prob_name(glp_prob *P);
/* retrieve problem name */

const(char)* glp_get_obj_name(glp_prob *P);
/* retrieve objective function name */

int glp_get_obj_dir(glp_prob *P);
/* retrieve optimization direction flag */

int glp_get_num_rows(glp_prob *P);
/* retrieve number of rows */

int glp_get_num_cols(glp_prob *P);
/* retrieve number of columns */

const(char)* glp_get_row_name(glp_prob *P, int i);
/* retrieve row name */

const(char)* glp_get_col_name(glp_prob *P, int j);
/* retrieve column name */

int glp_get_row_type(glp_prob *P, int i);
/* retrieve row type */

double glp_get_row_lb(glp_prob *P, int i);
/* retrieve row lower bound */

double glp_get_row_ub(glp_prob *P, int i);
/* retrieve row upper bound */

int glp_get_col_type(glp_prob *P, int j);
/* retrieve column type */

double glp_get_col_lb(glp_prob *P, int j);
/* retrieve column lower bound */

double glp_get_col_ub(glp_prob *P, int j);
/* retrieve column upper bound */

double glp_get_obj_coef(glp_prob *P, int j);
/* retrieve obj. coefficient or constant term */

int glp_get_num_nz(glp_prob *P);
/* retrieve number of constraint coefficients */

int glp_get_mat_row(glp_prob *P, int i, int* ind, double* val);
/* retrieve row of the constraint matrix */

int glp_get_mat_col(glp_prob *P, int j, int* ind, double* val);
/* retrieve column of the constraint matrix */

void glp_create_index(glp_prob *P);
/* create the name index */

int glp_find_row(glp_prob *P, const(char)* name);
/* find row by its name */

int glp_find_col(glp_prob *P, const(char)* name);
/* find column by its name */

void glp_delete_index(glp_prob *P);
/* delete the name index */

void glp_set_rii(glp_prob *P, int i, double rii);
/* set (change) row scale factor */

void glp_set_sjj(glp_prob *P, int j, double sjj);
/* set (change) column scale factor */

double glp_get_rii(glp_prob *P, int i);
/* retrieve row scale factor */

double glp_get_sjj(glp_prob *P, int j);
/* retrieve column scale factor */

void glp_scale_prob(glp_prob *P, int flags);
/* scale problem data */

void glp_unscale_prob(glp_prob *P);
/* unscale problem data */

void glp_set_row_stat(glp_prob *P, int i, int stat);
/* set (change) row status */

void glp_set_col_stat(glp_prob *P, int j, int stat);
/* set (change) column status */

void glp_std_basis(glp_prob *P);
/* construct standard initial LP basis */

void glp_adv_basis(glp_prob *P, int flags);
/* construct advanced initial LP basis */

void glp_cpx_basis(glp_prob *P);
/* construct Bixby's initial LP basis */

int glp_simplex(glp_prob *P, const(glp_smcp)* parm);
/* solve LP problem with the simplex method */

int glp_exact(glp_prob *P, const(glp_smcp)* parm);
/* solve LP problem in exact arithmetic */

void glp_init_smcp(glp_smcp *parm);
/* initialize simplex method control parameters */

int glp_get_status(glp_prob *P);
/* retrieve generic status of basic solution */

int glp_get_prim_stat(glp_prob *P);
/* retrieve status of primal basic solution */

int glp_get_dual_stat(glp_prob *P);
/* retrieve status of dual basic solution */

double glp_get_obj_val(glp_prob *P);
/* retrieve objective value (basic solution) */

int glp_get_row_stat(glp_prob *P, int i);
/* retrieve row status */

double glp_get_row_prim(glp_prob *P, int i);
/* retrieve row primal value (basic solution) */

double glp_get_row_dual(glp_prob *P, int i);
/* retrieve row dual value (basic solution) */

int glp_get_col_stat(glp_prob *P, int j);
/* retrieve column status */

double glp_get_col_prim(glp_prob *P, int j);
/* retrieve column primal value (basic solution) */

double glp_get_col_dual(glp_prob *P, int j);
/* retrieve column dual value (basic solution) */

int glp_get_unbnd_ray(glp_prob *P);
/* determine variable causing unboundedness */

int glp_get_it_cnt(glp_prob *P);
/* get simplex solver iteration count */

void glp_set_it_cnt(glp_prob *P, int it_cnt);
/* set simplex solver iteration count */

int glp_interior(glp_prob *P, const(glp_iptcp)* parm);
/* solve LP problem with the interior-point method */

void glp_init_iptcp(glp_iptcp *parm);
/* initialize interior-point solver control parameters */

int glp_ipt_status(glp_prob *P);
/* retrieve status of interior-point solution */

double glp_ipt_obj_val(glp_prob *P);
/* retrieve objective value (interior point) */

double glp_ipt_row_prim(glp_prob *P, int i);
/* retrieve row primal value (interior point) */

double glp_ipt_row_dual(glp_prob *P, int i);
/* retrieve row dual value (interior point) */

double glp_ipt_col_prim(glp_prob *P, int j);
/* retrieve column primal value (interior point) */

double glp_ipt_col_dual(glp_prob *P, int j);
/* retrieve column dual value (interior point) */

void glp_set_col_kind(glp_prob *P, int j, int kind);
/* set (change) column kind */

int glp_get_col_kind(glp_prob *P, int j);
/* retrieve column kind */

int glp_get_num_int(glp_prob *P);
/* retrieve number of integer columns */

int glp_get_num_bin(glp_prob *P);
/* retrieve number of binary columns */

int glp_intopt(glp_prob *P, const(glp_iocp)* parm);
/* solve MIP problem with the branch-and-bound method */

void glp_init_iocp(glp_iocp *parm);
/* initialize integer optimizer control parameters */

int glp_mip_status(glp_prob *P);
/* retrieve status of MIP solution */

double glp_mip_obj_val(glp_prob *P);
/* retrieve objective value (MIP solution) */

double glp_mip_row_val(glp_prob *P, int i);
/* retrieve row value (MIP solution) */

double glp_mip_col_val(glp_prob *P, int j);
/* retrieve column value (MIP solution) */

void glp_check_kkt(glp_prob *P, int sol, int cond, double *ae_max,
      int *ae_ind, double *re_max, int *re_ind);
/* check feasibility/optimality conditions */

int glp_print_sol(glp_prob *P, const(char)* fname);
/* write basic solution in printable format */

int glp_read_sol(glp_prob *P, const(char)* fname);
/* read basic solution from text file */

int glp_write_sol(glp_prob *P, const(char)* fname);
/* write basic solution to text file */

int glp_print_ranges(glp_prob *P, int len, const(int)* list,
      int flags, const(char)* fname);
/* print sensitivity analysis report */

int glp_print_ipt(glp_prob *P, const(char)* fname);
/* write interior-point solution in printable format */

int glp_read_ipt(glp_prob *P, const(char)* fname);
/* read interior-point solution from text file */

int glp_write_ipt(glp_prob *P, const(char)* fname);
/* write interior-point solution to text file */

int glp_print_mip(glp_prob *P, const(char)* fname);
/* write MIP solution in printable format */

int glp_read_mip(glp_prob *P, const(char)* fname);
/* read MIP solution from text file */

int glp_write_mip(glp_prob *P, const(char)* fname);
/* write MIP solution to text file */

int glp_bf_exists(glp_prob *P);
/* check if LP basis factorization exists */

int glp_factorize(glp_prob *P);
/* compute LP basis factorization */

int glp_bf_updated(glp_prob *P);
/* check if LP basis factorization has been updated */

void glp_get_bfcp(glp_prob *P, glp_bfcp* parm);
/* retrieve LP basis factorization control parameters */

void glp_set_bfcp(glp_prob *P, const(glp_bfcp)* parm);
/* change LP basis factorization control parameters */

int glp_get_bhead(glp_prob *P, int k);
/* retrieve LP basis header information */

int glp_get_row_bind(glp_prob *P, int i);
/* retrieve row index in the basis header */

int glp_get_col_bind(glp_prob *P, int j);
/* retrieve column index in the basis header */

void glp_ftran(glp_prob *P, double* x);
/* perform forward transformation (solve system B*x = b) */

void glp_btran(glp_prob *P, double* x);
/* perform backward transformation (solve system B'*x = b) */

int glp_warm_up(glp_prob *P);
/* "warm up" LP basis */

int glp_eval_tab_row(glp_prob *P, int k, int* ind, double* val);
/* compute row of the simplex tableau */

int glp_eval_tab_col(glp_prob *P, int k, int* ind, double* val);
/* compute column of the simplex tableau */

int glp_transform_row(glp_prob *P, int len, int* ind, double* val);
/* transform explicitly specified row */

int glp_transform_col(glp_prob *P, int len, int* ind, double* val);
/* transform explicitly specified column */

int glp_prim_rtest(glp_prob *P, int len, const(int)* ind,
      const(double)* val, int dir, double eps);
/* perform primal ratio test */

int glp_dual_rtest(glp_prob *P, int len, const(int)* ind,
      const(double)* val, int dir, double eps);
/* perform dual ratio test */

void glp_analyze_bound(glp_prob *P, int k, double *value1, int *var1,
      double *value2, int *var2);
/* analyze active bound of non-basic variable */

void glp_analyze_coef(glp_prob *P, int k, double *coef1, int *var1,
      double *value1, double *coef2, int *var2, double *value2);
/* analyze objective coefficient at basic variable */

int glp_ios_reason(glp_tree *T);
/* determine reason for calling the callback routine */

glp_prob *glp_ios_get_prob(glp_tree *T);
/* access the problem object */

void glp_ios_tree_size(glp_tree *T, int *a_cnt, int *n_cnt,
      int *t_cnt);
/* determine size of the branch-and-bound tree */

int glp_ios_curr_node(glp_tree *T);
/* determine current active subproblem */

int glp_ios_next_node(glp_tree *T, int p);
/* determine next active subproblem */

int glp_ios_prev_node(glp_tree *T, int p);
/* determine previous active subproblem */

int glp_ios_up_node(glp_tree *T, int p);
/* determine parent subproblem */

int glp_ios_node_level(glp_tree *T, int p);
/* determine subproblem level */

double glp_ios_node_bound(glp_tree *T, int p);
/* determine subproblem local bound */

int glp_ios_best_node(glp_tree *T);
/* find active subproblem with best local bound */

double glp_ios_mip_gap(glp_tree *T);
/* compute relative MIP gap */

void *glp_ios_node_data(glp_tree *T, int p);
/* access subproblem application-specific data */

void glp_ios_row_attr(glp_tree *T, int i, glp_attr *attr);
/* retrieve additional row attributes */

int glp_ios_pool_size(glp_tree *T);
/* determine current size of the cut pool */

int glp_ios_add_row(glp_tree *T,
      const(char)* name, int klass, int flags, int len, const(int)* ind,
      const(double) val, int type, double rhs);
/* add row (constraint) to the cut pool */

void glp_ios_del_row(glp_tree *T, int i);
/* remove row (constraint) from the cut pool */

void glp_ios_clear_pool(glp_tree *T);
/* remove all rows (constraints) from the cut pool */

int glp_ios_can_branch(glp_tree *T, int j);
/* check if can branch upon specified variable */

void glp_ios_branch_upon(glp_tree *T, int j, int sel);
/* choose variable to branch upon */

void glp_ios_select_node(glp_tree *T, int p);
/* select subproblem to continue the search */

int glp_ios_heur_sol(glp_tree *T, const(double) x);
/* provide solution found by heuristic */

void glp_ios_terminate(glp_tree *T);
/* terminate the solution process */

void glp_init_mpscp(glp_mpscp *parm);
/* initialize MPS format control parameters */

int glp_read_mps(glp_prob *P, int fmt, const(glp_mpscp)* parm,
      const(char)* fname);
/* read problem data in MPS format */

int glp_write_mps(glp_prob *P, int fmt, const(glp_mpscp)* parm,
      const(char)* fname);
/* write problem data in MPS format */

void glp_init_cpxcp(glp_cpxcp *parm);
/* initialize CPLEX LP format control parameters */

int glp_read_lp(glp_prob *P, const(glp_cpxcp)* parm, const(char)* fname);
/* read problem data in CPLEX LP format */

int glp_write_lp(glp_prob *P, const(glp_cpxcp)* parm, const(char)* fname);
/* write problem data in CPLEX LP format */

int glp_read_prob(glp_prob *P, int flags, const(char)* fname);
/* read problem data in GLPK format */

int glp_write_prob(glp_prob *P, int flags, const(char)* fname);
/* write problem data in GLPK format */

glp_tran *glp_mpl_alloc_wksp();
/* allocate the MathProg translator workspace */

int glp_mpl_read_model(glp_tran *tran, const(char)* fname, int skip);
/* read and translate model section */

int glp_mpl_read_data(glp_tran *tran, const(char)* fname);
/* read and translate data section */

int glp_mpl_generate(glp_tran *tran, const(char)* fname);
/* generate the model */

void glp_mpl_build_prob(glp_tran *tran, glp_prob *prob);
/* build LP/MIP problem instance from the model */

int glp_mpl_postsolve(glp_tran *tran, glp_prob *prob, int sol);
/* postsolve the model */

void glp_mpl_free_wksp(glp_tran *tran);
/* free the MathProg translator workspace */

int glp_main(int argc, const(char)** argv);
/* stand-alone LP/MIP solver */

int glp_read_cnfsat(glp_prob *P, const(char)* fname);
/* read CNF-SAT problem data in DIMACS format */

int glp_check_cnfsat(glp_prob *P);
/* check for CNF-SAT problem instance */

int glp_write_cnfsat(glp_prob *P, const(char)* fname);
/* write CNF-SAT problem data in DIMACS format */

int glp_minisat1(glp_prob *P);
/* solve CNF-SAT problem with MiniSat solver */

int glp_intfeas1(glp_prob *P, int use_bound, int obj_bound);
/* solve integer feasibility problem */

int glp_init_env();
/* initialize GLPK environment */

const(char)* glp_version();
/* determine library version */

int glp_free_env();
/* free GLPK environment */

void glp_puts(const(char)* s);
/* write string on terminal */

void glp_printf(const(char)* fmt, ...);
/* write formatted output on terminal */

void glp_vprintf(const(char)* fmt, va_list arg);
/* write formatted output on terminal */

int glp_term_out(int flag);
/* enable/disable terminal output */

void glp_term_hook(int function(void *info, const(char)* s) func, void *info);
/* install hook to intercept terminal output */

int glp_open_tee(const(char)* name);
/* start copying terminal output to text file */

int glp_close_tee();
/* stop copying terminal output to text file */

alias glp_errfunc = void function(const(char)* fmt, ...);

glp_errfunc glp_error_(const(char)* file, int line);
/* display fatal error message and terminate execution */

//alias glp_error = glp_error_(__FILE__, __LINE__);

void glp_error(const(char)* file, int line){
  glp_error_(__FILE__, __LINE__);
}

int glp_at_error();

void glp_assert_(const(char)* expr, const(char)* file, int line);
/* check for logical condition */

void glp_error_hook(void function(void *info) func, void *info);
/* install hook to intercept abnormal termination */

void *glp_alloc(int n, int size);
/* allocate memory block */

//alias glp_malloc(size) = glp_alloc(1, size);
alias glp_calloc = glp_alloc;

void glp_malloc(int size){
  glp_error_(1, size);
}

void *glp_realloc(void *ptr, int n, int size);
/* reallocate memory block */

void glp_free(void *ptr);
/* free (deallocate) memory block */

void glp_mem_limit(int limit);
/* set memory usage limit */

void glp_mem_usage(int *count, int *cpeak, size_t *total, size_t *tpeak);
/* get memory usage information */

//alias glp_graph = void*;
//alias glp_vertex = void*;
//alias glp_arc = void*;

struct glp_graph
{     /* graph descriptor */
      void *pool; /* DMP *pool; */
      /* memory pool to store graph components */
      char *name;
      /* graph name (1 to 255 chars); NULL means no name is assigned
         to the graph */
      int nv_max;
      /* length of the vertex list (enlarged automatically) */
      int nv;
      /* number of vertices in the graph, 0 <= nv <= nv_max */
      int na;
      /* number of arcs in the graph, na >= 0 */
      glp_vertex** v; /* glp_vertex *v[1+nv_max]; */
      /* v[i], 1 <= i <= nv, is a pointer to i-th vertex */
      void *index; /* AVL *index; */
      /* vertex index to find vertices by their names; NULL means the
         index does not exist */
      int v_size;
      /* size of data associated with each vertex (0 to 256 bytes) */
      int a_size;
      /* size of data associated with each arc (0 to 256 bytes) */
}

struct glp_vertex
{     /* vertex descriptor */
      int i;
      /* vertex ordinal number, 1 <= i <= nv */
      char *name;
      /* vertex name (1 to 255 chars); NULL means no name is assigned
         to the vertex */
      void *entry; /* AVLNODE *entry; */
      /* pointer to corresponding entry in the vertex index; NULL means
         that either the index does not exist or the vertex has no name
         assigned */
      void *data;
      /* pointer to data associated with the vertex */
      void *temp;
      /* working pointer */
      glp_arc* in_;
      /* pointer to the (unordered) list of incoming arcs */
      glp_arc* out_;
      /* pointer to the (unordered) list of outgoing arcs */
}

struct glp_arc
{     /* arc descriptor */
      glp_vertex *tail;
      /* pointer to the tail endpoint */
      glp_vertex *head;
      /* pointer to the head endpoint */
      void *data;
      /* pointer to data associated with the arc */
      void *temp;
      /* working pointer */
      glp_arc* t_prev;
      /* pointer to previous arc having the same tail endpoint */
      glp_arc* t_next;
      /* pointer to next arc having the same tail endpoint */
      glp_arc* h_prev;
      /* pointer to previous arc having the same head endpoint */
      glp_arc* h_next;
      /* pointer to next arc having the same head endpoint */
}

glp_graph* glp_create_graph(int v_size, int a_size);
/* create graph */

void glp_set_graph_name(glp_graph *G, const(char)* name);
/* assign (change) graph name */

int glp_add_vertices(glp_graph *G, int nadd);
/* add new vertices to graph */

void glp_set_vertex_name(glp_graph *G, int i, const(char)* name);
/* assign (change) vertex name */

glp_arc* glp_add_arc(glp_graph *G, int i, int j);
/* add new arc to graph */

void glp_del_vertices(glp_graph *G, int ndel, const(int)* num);
/* delete vertices from graph */

void glp_del_arc(glp_graph *G, glp_arc *a);
/* delete arc from graph */

void glp_erase_graph(glp_graph *G, int v_size, int a_size);
/* erase graph content */

void glp_delete_graph(glp_graph *G);
/* delete graph */

void glp_create_v_index(glp_graph *G);
/* create vertex name index */

int glp_find_vertex(glp_graph *G, const(char)* name);
/* find vertex by its name */

void glp_delete_v_index(glp_graph *G);
/* delete vertex name index */

int glp_read_graph(glp_graph *G, const(char)* fname);
/* read graph from plain text file */

int glp_write_graph(glp_graph *G, const(char)* fname);
/* write graph to plain text file */

void glp_mincost_lp(glp_prob *P, glp_graph *G, int names, int v_rhs,
      int a_low, int a_cap, int a_cost);
/* convert minimum cost flow problem to LP */

int glp_mincost_okalg(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, double *sol, int a_x, int v_pi);
/* find minimum-cost flow with out-of-kilter algorithm */

int glp_mincost_relax4(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, int crash, double *sol, int a_x, int a_rc);
/* find minimum-cost flow with Bertsekas-Tseng relaxation method */

void glp_maxflow_lp(glp_prob *P, glp_graph *G, int names, int s,
      int t, int a_cap);
/* convert maximum flow problem to LP */

int glp_maxflow_ffalg(glp_graph *G, int s, int t, int a_cap,
      double *sol, int a_x, int v_cut);
/* find maximal flow with Ford-Fulkerson algorithm */

int glp_check_asnprob(glp_graph *G, int v_set);
/* check correctness of assignment problem data */

/* assignment problem formulation: */
enum GLP_ASN_MIN = 1;  /* perfect matching (minimization) */
enum GLP_ASN_MAX = 2;  /* perfect matching (maximization) */
enum GLP_ASN_MMP = 3;  /* maximum matching */

int glp_asnprob_lp(glp_prob *P, int form, glp_graph *G, int names, int v_set, int a_cost);
/* convert assignment problem to LP */

int glp_asnprob_okalg(int form, glp_graph *G, int v_set, int a_cost, double *sol, int a_x);
/* solve assignment problem with out-of-kilter algorithm */

int glp_asnprob_hall(glp_graph *G, int v_set, int a_x);
/* find bipartite matching of maximum cardinality */

double glp_cpp(glp_graph *G, int v_t, int v_es, int v_ls);
/* solve critical path problem */

int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const(char)* fname);
/* read min-cost flow problem data in DIMACS format */

int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const(char)* fname);
/* write min-cost flow problem data in DIMACS format */

int glp_read_maxflow(glp_graph *G, int *s, int *t, int a_cap,
      const(char)* fname);
/* read maximum flow problem data in DIMACS format */

int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
      const(char)* fname);
/* write maximum flow problem data in DIMACS format */

int glp_read_asnprob(glp_graph *G, int v_set, int a_cost, const(char)* fname);
/* read assignment problem data in DIMACS format */

int glp_write_asnprob(glp_graph *G, int v_set, int a_cost, const(char)* fname);
/* write assignment problem data in DIMACS format */

int glp_read_ccdata(glp_graph *G, int v_wgt, const(char)* fname);
/* read graph in DIMACS clique/coloring format */

int glp_write_ccdata(glp_graph *G, int v_wgt, const(char)* fname);
/* write graph in DIMACS clique/coloring format */

int glp_netgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
      const(int) parm[1+15]);
/* Klingman's network problem generator */

void glp_netgen_prob(int nprob, int parm[1+15]);
/* Klingman's standard network problem instance */

int glp_gridgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
      const(int) parm[1+14]);
/* grid-like network problem generator */

int glp_rmfgen(glp_graph *G, int *s, int *t, int a_cap,
      const(int) parm[1+5]);
/* Goldfarb's maximum flow problem generator */

int glp_weak_comp(glp_graph *G, int v_num);
/* find all weakly connected components of graph */

int glp_strong_comp(glp_graph *G, int v_num);
/* find all strongly connected components of graph */

int glp_top_sort(glp_graph *G, int v_num);
/* topological sorting of acyclic digraph */

int glp_wclique_exact(glp_graph *G, int v_wgt, double *sol, int v_set);
/* find maximum weight clique with exact algorithm */

}

