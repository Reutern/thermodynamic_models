#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H

#include "SeqAnnotator.h"
#include <signal.h>

enum ModelType {
    LOGISTIC,   // logistic regression
    DIRECT,     // direct interaction between TF and BTM, repressor works through BTM
    QUENCHING,      // repressor stops activator from interacting with BTM
    CHRMOD_UNLIMITED,     // repressor works by chromatin modification (making it unaccessible), unlimited activation
    CHRMOD_LIMITED        // repressor works by chromatin modification (making it unaccessible), limited activation
};

ModelType getModelOption( const string& modelOptionStr );
string getModelOptionStr( ModelType modelOption ); 

enum FactorIntType {
    BINARY,     // Binary model of interaction
	LINEAR,		// Linear model of interaction
    GAUSSIAN    // Gaussian model of interaction
}; 

FactorIntType getIntOption( const string& intOptionStr );
string getIntOptionStr( FactorIntType intOption );

enum ObjType {
    SSE,    // sum of squared error
    CORR,   // Pearson correlation
    PGP         // PGP score
};

enum PenaltyType {
    NONE,   // no parameter penalty
    L1,    // L1 norm 
    L2,    // L2 norm
};

ObjType getObjOption( const string& objOptionStr );
string getObjOptionStr( ObjType objOption );
PenaltyType getPenaltyOption( const string& PenaltyOptionStr );
string getPenaltyOptionStr( PenaltyType PenaltyOption );

enum SearchType {
    UNCONSTRAINED,  // unconstrained search
    CONSTRAINED     // constrained search
};

string getSearchOptionStr( SearchType searchOption );

/*****************************************************
* Expression Model and Parameters
******************************************************/

/* ExprPar class: the parameters of the expression model */
class ExprPar {
public:
    // constructors 
    ExprPar() : factorIntMat() {}
    ExprPar( int _nFactors, int _nSeqs );		// default values of parameters
    ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >&  _basalTxps, int _nSeqs, double _acc_scale, double _par_penalty );
    ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs );	// construct from a "flat" vector of free parameters (assuming they are in the correct/uniform scale)
    void copy( const ExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; basalTxps = other.basalTxps; nSeqs = basalTxps.size(); acc_scale = other.acc_scale; par_penalty = other.par_penalty; }
    ExprPar( const ExprPar& other ) { copy( other ); }

    // assignment
    ExprPar& operator=( const ExprPar& other ) { copy( other ); return *this; }	
	
    // access methods
    int nFactors() const { return maxBindingWts.size(); }

	// parameter norm
	double parameter_L2_norm() const;
	double parameter_L1_norm() const;
	
    // get the free parameters (in the correct/uniform scale)
    void getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const; 
	
    // print the parameters
    void print( ostream& os, const vector< string >& motifNames, const vector< string >& seqNames, const IntMatrix& coopMat ) const;
 		
    // load the parameter values from a file, assuming the parameter has the correct dimensions (and initialized)
    int load( const string& file ); 
    int load( const string& file, const vector <string>& seqNames, const vector <string>& motifNames);	// extended version of load with automated seqName allocation for q_btm 

    // adjust the values of parameters: if the value is close to min or max allowed value, slightly change it s.t. it is away from the boundary
    void adjust(); 

    // constrain the parameters between the minimum and maximum value
    void constrain_parameters();

    // parameters
    vector< double > maxBindingWts;			// binding weight of the strongest site for each TF: K(S_max) [TF_max]
    Matrix factorIntMat; 		// (maximum) interactions between pairs of factors: omega(f,f')
    vector< double > txpEffects;    // transcriptional effects: alpha for Direct and Quenching model, exp(alpha) for Logistic model (so that the same default values can be used). Equal to 1 if a TF is not an activator under the Quenching model
    vector< double > repEffects;    // repression effects: beta under ChrMod models (the equlibrium constant of nucleosome association with chromatin). Equal to 0 if a TF is not a repressor. 
    vector < double > basalTxps;        // basal transcription: q_p for Direct and Quenching model, exp(alpha_0) for Logistic model (so that the same default value can be used)
	//double expRatio; 		// constant factor of measurement to prediction 
    double acc_scale;
	double par_penalty;

	int nSeqs;

    static ModelType modelOption;     // model option
    static SearchType searchOption;    // search option: 0 - unconstrained search; 1 - constrained search
    static int estBindingOption;    // whether to estimate binding parameters
    static bool one_qbtm_per_crm;
    
    static double default_acc_scale;	// default accessibility scaling parameter
    static double default_par_penalty;	// default parameter penalty
    static double default_weight;	// default binding weight
    static double default_interaction;		// default factor interaction
    static double default_effect_Logistic;   // default transcriptional effect under Logistic model
    static double default_effect_Thermo;     // default transcriptional effect under thermo. models
    static double default_repression;   // default repression
    static double default_basal_Logistic;       // default basal transcription under Logistic model
    static double default_basal_Thermo;         // default basal transcriptional under thermo. models
    static double min_acc_scale;	// min accessibility scaling parameter
    static double max_acc_scale;	// max accessibility scaling parameter
    static double min_weight;		// min. binding weight
    static double max_weight;		// max. binding weight
    static double min_interaction;	    // min. interaction
    static double max_interaction;	    // max. interaction
    static double min_effect_Logistic;   // min. transcriptional effect under Logistic model
    static double max_effect_Logistic;   // max. transcriptional effect under Logistic model
//     static double min_effect_Direct;   // min. transcriptional effect under Direct model
    static double min_effect_Thermo;    // min. transcriptional effect under thermo. models
    static double max_effect_Thermo;   // max. transcriptional effect under thermo. models   
    static double min_repression;   // min. repression
    static double max_repression;   // max. repression
    static double min_basal_Logistic;    // min. basal transcription under Logistic model
    static double max_basal_Logistic;    // max. basal transcription under Logistic model
    static double min_basal_Thermo;   // min. basal transcription under thermo. models
    static double max_basal_Thermo;   // max. basal transcription under thermo. models
    static double delta;        // a small number for testing the range of parameter values
// 	static double wt_step;		// step of maxExprWt (log10)
// 	static double int_step;		// step of interaction (log10)
// 	static double ratio_step;	// step of expRatio
};

/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc {
public:
    // constructors
    ExprFunc( const vector< Motif >& _motifs,  const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr,  int _coopDistThr, const ExprPar& _par, const vector< Sequence >& _seqs  );

    // access methods
    const vector< Motif >& getMotifs() const {
        return motifs;
    }
    
    // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
    double predictExpr( int length, const vector< double >& factorConcs, int seq_num );
    double predictExpr( int length, const vector< double >& factorConcs, int seq_num, std::ofstream& fout );
    double predictExpr_scanning_mode( int length, const vector< double >& factorConcs, int seq_num );
    // Returns the efficiency Z_ON/Z_OFF
    double predictExpr_scalefree( int length, const vector< double >& factorConcs, int seq_num );    
    static ModelType modelOption;     // model option   
	static FactorIntType FactorIntOption;
    static bool one_qbtm_per_crm;

    // Access functions to privat variables
    void set_sites(SiteVec _sites); 
    void set_boundaries( vector< int > _boundaries) {boundaries = _boundaries;}
    void set_bindingWts( vector< double > _bindingWts) {bindingWts = _bindingWts;}

    // TF binding motifs
    const vector< Motif >& motifs; 	

    // control parameters
    const vector< bool >& actIndicators;    // 1 if the TF is in the activator set
    int maxContact;     // the maximum contact     
    const vector< bool >& repIndicators;    // 1 if the TF is in the repressor set
    const IntMatrix& repressionMat;    // repression matrix: R(f,f') = 1 if f can repress f'
    int repressionDistThr;   // distance threshold for repression: d_R
    int coopDistThr;    // distance threshold for interaction	
		
    // model parameters
    const ExprPar& par;

    // the sequences
    const vector< Sequence >& seqs;
		    
    // the sequence whose expression is to be predicted
    SiteVec sites;
    vector< int > boundaries;   // left boundary of each site beyond which there is no interaction    

    // intermediate computational results
    vector< double > bindingWts; 	// binding weights without the concentration exp(-E)
		
    // compute the partition function when the basal transcriptional machinery (BTM) is not bound
    double compPartFuncOff( const vector< double >& factorConcs) const;

    void compProb_scanning_mode( const vector< double >& factorConcs, vector< double >& p_bound);

    // compute the partition function with BTM bound and unbound on a sequence basis 
    int compPartFunc_seq(double &result_Z_on, double &result_Z_off, int seq_num, const vector< double >& factorConcs) const;

    // compute the partition function with BTM bound and unbound on a sequence basis (with interfactor interaction)
    int compPartFunc_seq_interfactor(double &result_Z_on, double &result_Z_off, int seq_num, const vector< double >& factorConcs) const;

    // compute the partition function when the BTM is not bound: ChrMod model 
    double compPartFuncOffChrMod( const vector< double >& factorConcs ) const; 
    
    // compute the partition function when the BTM is bound 
    double compPartFuncOn( const vector< double >& factorConcs) const;

    // compute the paritition function when the BTM is bound: Direct model
    double compPartFuncOnDirect( const vector< double >& factorConcs) const;
    
    // compute the paritition function when the BTM is bound: Quenching model
    double compPartFuncOnQuenching( const vector< double >& factorConcs) const;

     // compute the paritition function when the BTM is bound: ChrMod_Unlimited model
    double compPartFuncOnChrMod_Unlimited( const vector< double >& factorConcs) const;

     // compute the paritition function when the BTM is bound: ChrMod_Limited model
    double compPartFuncOnChrMod_Limited( const vector< double >& factorConcs) const;    


   private: 
    // compute the TF-TF interaction between two occupied sites
    double compFactorInt( const Site& a, const Site& b ) const;
    double compFactorInt( int t_1, int t_2, int _dist  ) const;

    // test if one site represses another site
    bool testRepression( const Site& a, const Site& b ) const;
};

/*****************************************************
* Model Training and Testing
******************************************************/

/* ExprPredictor class: the thermodynamic sequence-to-expression predictor */
class ExprPredictor {
public:
    // constructors
    ExprPredictor( const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr,const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector <string>& _seqNames, const vector< Sequence >& _seqs  );

    // access methods
    int nSeqs() const {
        return seqSites.size();
    }
    int nFactors() const { 
        return motifs.size(); 
    }
    int nConds() const {
         return exprData.nCols();
    }
    static const IntMatrix& getCoopMat() {
        return coopMat;
    }
    const vector< bool >& getActIndicators() const {
        return actIndicators;
    }
    const vector< bool >& getRepIndicators() const {
        return repIndicators;
    }
    const IntMatrix& getRepressionMat() const {
        return repressionMat;
    }
    const ExprPar& getPar() { return par_model; }
    const void setPar(ExprPar& par) { par_model = par; }
    double getObj() const { return obj_model; }
    
    // the objective function to be minimized
    double objFunc( const ExprPar& par ) ;
    double objFunc( const ExprPar& par, int crm ) ;
    // training the model
    int train( const ExprPar& par_init ); 	// training with the initial values given
    int train( const ExprPar& par_init, const gsl_rng* rng );   // training with the initial values and allowing random starts
    int train();	// automatic training: first estimate the initial values, then train
            
    // predict expression values of a sequence (across the same conditions)
    int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const; 

    double comp_impact( const ExprPar& par, int tf );		// The impact of the parameter
    double comp_impact_coop( const ExprPar& par, int tf );		// The impact of all cooperativity parameters with tf involved
    double comp_impact_acc( const ExprPar& par );		// The impact of the accessibility
    double comp_impact( const ExprPar& par, int tf, int crm );		// The impact of the parameter on one crm
    double comp_impact_coop( const ExprPar& par, int tf, int crm );		// The impact of all cooperativity parameters with tf involved on one crm
    double comp_impact_coop_pair( const ExprPar& par, int tf1, int tf2 );		// The impact of the cooperativity parameter

    // test the model, perfOption = 0: RMSE
// 	double test( const vector< Sequence  >& testSeqs, const Matrix& testExprData, Matrix& predictions ) const;    

	std::ofstream gene_crm_fout;

    static ModelType modelOption;     // model option
    static int estBindingOption;    // whether estimate binding parameters
    static ObjType objOption;       // option of the objective function
    static PenaltyType PenaltyOption;       // option of the objective function

    // the similarity between two expression patterns, using cross-correlation
    static double exprSimCrossCorr( const vector< double >& x, const vector< double >& y ); 
    static int maxShift;    // maximum shift when computing cross correlation
    static double shiftPenalty;     // the penalty for shift (when weighting different positions)

    // the parameters for the optimizer
    static int nAlternations;   // number of alternations (between two optimization methods)
    static int nRandStarts;     // number of random starts
    static double min_delta_f_SSE;      // the minimum change of the objective function under SSE
    static double min_delta_f_Corr;     // the minimum change of the objective function under correlation
    static double min_delta_f_PGP;            // the minimum change of the objective function under PGP
    static int nSimplexIters;       // maximum number of iterations for Simplex optimizer
    static int nGradientIters;      // maximum number of iterations for Gradient optimizer
    static bool one_qbtm_per_crm;
    vector < bool > indicator_bool;	// States if Parameters are free or fixed for training
    static vector <string> motifNames;
    static vector <string> seqNames;
    vector < double > fix_pars;
    vector < double > free_pars;

    // print the parameter values (the ones that are estimated) in a single line
    void printPar( const ExprPar& par ) const;
    static ExprPar par_curr;

private:
    // training data
    const vector< SiteVec >& seqSites;		// the extracted sites for all sequences
    const vector< int >& seqLengths;           // lengths of all sequences
    const Matrix& exprData;		// expressions of the corresponding sequences across multiple conditions         
    const vector< Motif >& motifs;		// TF binding motifs
    const Matrix& factorExprData;		// [TF] of all factors over multiple conditions	    

    // control parameters 
    static IntMatrix coopMat;       // cooperativity matrix: C(f,f') = 1 if f and f' bind cooperatively    
    const vector< bool >& actIndicators;   // 1 if the TF is in the activator set
    int maxContact;     // the maximum contact     
    const vector< bool >& repIndicators;    // 1 if the TF is in the repressor set
    const IntMatrix& repressionMat;    // repression matrix: R(f,f') = 1 if f can repress f'
    int repressionDistThr;   // distance threshold for repression: d_R
    int coopDistThr;   // distance threshold for cooperativity
    
    // model parameters and the value of the objective function
    ExprPar par_model;
    double obj_model;	
    double obj_pgp;
    double obj_corr;	
    double obj_sse;	

    // the sequence
    const vector< Sequence >& seqs;

    // randomly sample parameter values (only those free parameters), the parameters should be initialized
    int randSamplePar( const gsl_rng* rng, ExprPar& par ) const; 

    // check if some parameter combination is valid
    bool testPar( const ExprPar& par ) const; 
    
    // create the expression function
    ExprFunc* createExprFunc( const ExprPar& par ) const;
    
    // objective functions
    double compAvgCorr( const ExprPar& par );     	// the average Pearson correlation		Outdated !!
    double compAvgCrossCorr( const ExprPar& par );    	// the average cross correlation -based similarity		Outdated !!

    // minimize the objective function, using the current model parameters as initial values
    int simplex_minimize( ExprPar& par_result, double& obj_result );	// simplex	
    int gradient_minimize( ExprPar& par_result, double& obj_result );	// gradient: BFGS or conjugate gradient
    int simulated_annealing( ExprPar& par_result, double& obj_result ); 
    // function to save parameters to file
    static int save_param();
    // Signal handler
    static void catch_signal(int sig_num);
 //   static void catch_param(int sig_num);		
};

double gsl_obj_f( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad ); 
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad ); 

#endif
