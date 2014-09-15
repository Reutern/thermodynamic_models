#ifndef OCC_PREDICTOR_H
#define OCC_PREDICTOR_H

#include "SeqAnnotator.h"
#include "ExprPredictor.h"
#include <signal.h>


/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class OccPredictor {
public:
    // constructors
    OccPredictor(const SiteVec & _sites, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, int _coopDistThr, const ExprPar& _par );

    // access methods
    const vector< Motif >& getMotifs() const {
        return motifs;
    }
    
    // predict the occupancy value of a given sequence under given TF concentrations
    double predictOcc( );

    // Access functions to privat variables
    void set_sites( SiteVec _sites) {sites = _sites;}
    void set_boundaries( vector< int > _boundaries) {boundaries = _boundaries;}
    void set_bindingWts( vector< double > _bindingWts) {bindingWts = _bindingWts;}

private:
    // TF binding motifs
    const vector< Motif >& motifs; 	

    // TF concentrations of all factors over multiple conditions	
    const Matrix& factorExprData;		    

    // control parameters
    static IntMatrix coopMat;       // cooperativity matrix: C(f,f') = 1 if f and f' bind cooperatively    
    int coopDistThr;    // distance threshold for interaction	
		
    // model parameters
    const ExprPar& par;

    // the sequence sites whose expression is to be predicted
    SiteVec sites;
    vector< int > boundaries;   // left boundary of each site beyond which there is no interaction    

    // intermediate computational results
    vector< double > bindingWts; 	// binding weights without the concentration exp(-E)
		
    // compute the partition functions
    double compPartFunc( const vector< double >& factorConcs, int centre) const;
    double compPartFunc_total( const vector< double >& factorConcs) const;
   
    // compute the TF-TF interaction between two occupied sites
    double compFactorInt( const Site& a, const Site& b ) const;
    double compFactorInt( int t_1, int t_2, int _dist  ) const;


};

#endif
