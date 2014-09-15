#include "OccPredictor.h"
#include "ExprPredictor.h"
#include "param.h"
#include <omp.h>
#include <unistd.h>
#include <set>


OccPredictor::OccPredictor(const SiteVec & _sites, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, int _coopDistThr, const ExprPar& _par ) : sites( _sites ), motifs( _motifs ), factorExprData( _factorExprData ), coopDistThr(_coopDistThr), par( _par )
{
    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors ); 
    int n = sites.size();
    vector< double > _bindingWts (n,0.0);
    vector< int > _boundaries (n,0);
	
    // Determin the boundaries for occpred
    _boundaries[0] = 0;
    int range = coopDistThr;
    for ( int k = 1; k < n; k++ ) {
    	int l; 
	for ( l=0; l < k; l++ ) {
	    if ( ( sites[k].start - sites[l].start ) <= range ) {break;} 
	}
    _boundaries[k] = l ;
    }	
    set_boundaries(_boundaries);
    
    // compute the Boltzman weights of binding for all sites for occpred
    _bindingWts[0] = 1.0;
    for ( int k = 1; k < n; k++ ) {
        _bindingWts[k] = par.maxBindingWts[ sites[k].factorIdx ] * sites[k].wtRatio ;	
    }
    set_bindingWts(_bindingWts); 

}

double OccPredictor::predictOcc( )
{
    vector< double > factorConcs = factorExprData.getCol( 50 );
    double Z_total = 0;
    double Z = 0;

    Z_total = compPartFunc_total(factorConcs);
    Z = compPartFunc(factorConcs, 10);


    // compute the expression (promoter occupancy)
    double occupancy = Z / Z_total;
    cout << occupancy <<endl;
    return occupancy;
}

// Calculate the binding weight for all configuration
double OccPredictor::compPartFunc_total(const vector< double >& factorConcs) const
{
        
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence 
    for ( int i = 1; i <= n; i++ ) {
	double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];	
        }
        Z[i] = bindingWts[ i ] * factorConcs[sites[ i ].factorIdx] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }
    return Zt[n];
}


// Calculate the binding weight for all configuration with binding site centre occupied
double OccPredictor::compPartFunc(const vector< double >& factorConcs, int centre) const
{
        
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence 
    for ( int i = 1; i <= n; i++ ) {
	double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];	
        }
        Z[i] = bindingWts[ i ] * factorConcs[sites[ i ].factorIdx] * sum;
	if(i != centre)
      		  Zt[i] = Z[i] + Zt[i - 1];		// normal partition function calculation  
	else
      		  Zt[i] = Z[i];				// binding site centre is always occupied 
    }

    return Zt[n];
}


double OccPredictor::compFactorInt( const Site& a, const Site& b ) const
{

// 	assert( !siteOverlap( a, b, motifs ) );
    if(a.factorIdx != b.factorIdx)	return 1.0;	// Only TF of the same type interact 	
    double maxInt = par.factorIntMat( a.factorIdx, b.factorIdx );
    unsigned dist = abs( a.start - b.start );
    assert( dist >= 0 );

    #if FactorIntFunc
    double spacingTerm = ( dist < coopDistThr ? maxInt *  (1 - dist/coopDistThr ) + 1.0: 1.0 );
    #else 
    double spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 );
    #endif // FactorIntFunc

    #if ORIENTATION
    double orientationTerm = ( a.strand == b.strand ) ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
    #else
    return spacingTerm;
    #endif //ORIENTATION
}

double OccPredictor::compFactorInt( int t_1, int t_2, int _dist  ) const
{

    double maxInt = par.factorIntMat( t_1, t_2 );
    int dist = _dist;
    assert( dist >= 0 );

    #if FactorIntFunc
    double spacingTerm = ( dist < coopDistThr ? maxInt *  (1 - dist/coopDistThr ) + 1.0: 1.0 );
    #else
    double spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 );
    #endif // FactorIntFunc

    #if ORIENTATION
    double orientationTerm = ( a.strand == b.strand ) ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
    #else
    return spacingTerm;
    #endif //ORIENTATION
}

