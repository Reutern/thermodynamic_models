#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include "ExprPredictor.h"
#include "param.h"
#include <sys/time.h>


ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;

    cerr << "modelOptionStr is not a valid model option" << endl; 
    exit(1);
}

string getModelOptionStr( ModelType modelOption )
{
    if ( modelOption == LOGISTIC ) return "Logisitic";
    if ( modelOption == DIRECT ) return "Direct";
    if ( modelOption == QUENCHING ) return "Quenching";
    if ( modelOption == CHRMOD_UNLIMITED ) return "ChrMod_Unlimited";
    if ( modelOption == CHRMOD_LIMITED ) return "ChrMod_Limited";

    return "Invalid";
}

string getIntOptionStr( FactorIntType intOption )
{
    if ( intOption == BINARY ) return "Binary";
    if ( intOption == GAUSSIAN ) return "Gaussian";

    return "Invalid";
}

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "NORM_CORR" ) return NORM_CORR;

    cerr << "objOptionStr is not a valid option of objective function" << endl; 
    exit(1);
}

string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == NORM_CORR ) return "Norm_Corr";

    return "Invalid";
}

string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}

double FactorIntFuncBinary::compFactorInt( double normalInt, int dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;	
}

double FactorIntFuncGaussian::compFactorInt( double normalInt, int dist, bool orientation ) const
{
    assert( dist >= 0 );

    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );    
}

double FactorIntFuncGeometric::compFactorInt( double normalInt, int dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) ); 
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;
}

ExprPar::ExprPar( int _nFactors, int _nSeqs ) : factorIntMat()
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );       

    for ( int i = 0; i < _nFactors; i++ ) {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
        repEffects.push_back( ExprPar::default_repression );
    }

	nSeqs = _nSeqs;
	if( one_qbtm_per_crm  ){
		for( int i = 0; i < nSeqs; i++ ){
    			double basalTxp_val = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo; 
    			basalTxps.push_back( basalTxp_val ); 
    		}
	}
	else{
		double basalTxp_val = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo; 
    		basalTxps.push_back( basalTxp_val ); 
	}
}
	
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, int _nSeqs ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), nSeqs( _nSeqs  )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() ); 	
    assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
    	if ( one_qbtm_per_crm ){
    		assert( basalTxps.size() == nSeqs );
	}
	else{
		assert( basalTxps.size() == 1 );
	}
}

ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs ) : factorIntMat()
{	

    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
//     assert( pars.size() == ( _nFactors * ( _nFactors + 1 ) / 2 + 2 * _nFactors + 2 ); 
    int counter = 0;
	
    // set maxBindingWts 
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    
    // set the interaction matrix
    factorIntMat.setDimensions( _nFactors, _nFactors );
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
                double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );
                factorIntMat( i, j ) = interaction;  
            }
            else factorIntMat( i, j ) = ExprPar::default_interaction;
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {
            factorIntMat( i, j ) = factorIntMat( j, i );
        }
    }       

    // set the transcriptional effects
    for ( int i = 0; i < _nFactors; i++ ) {
//         double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo; 
        if ( modelOption == LOGISTIC ) {
            double effect = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_effect_Logistic, max_effect_Logistic ) : pars[counter++];
            txpEffects.push_back( effect );
        } else if ( modelOption == DIRECT ) {
            double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
            txpEffects.push_back( effect ); 
        } else {
            if ( actIndicators[i] ) {
                double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
                txpEffects.push_back( effect );
            } else {
                txpEffects.push_back( ExprPar::default_effect_Thermo );
            }
        }
    }

    // set the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            if ( repIndicators[i] ) {
                double repression = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_repression ), log( max_repression ) ) ) : exp( pars[counter++] );
                repEffects.push_back( repression );
            } else {
                repEffects.push_back( ExprPar::default_repression ); 
            }
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) repEffects.push_back( ExprPar::default_repression );
    }
    
    // set the basal transcription
if( one_qbtm_per_crm ){	
	nSeqs = _nSeqs;
	for( int i = 0; i < nSeqs; i++ ){
		if ( modelOption == LOGISTIC ) {
        		double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
        		basalTxps.push_back( basal );
    		} else {
        		double basal = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_basal_Thermo ), log( max_basal_Thermo ) ) ) : exp( pars[counter++] );
        		basalTxps.push_back( basal );
    		}
	}
	}
	else{
	
		if ( modelOption == LOGISTIC ) {
        		double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
        		basalTxps.push_back( basal );
    		} else {
        		double basal = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_basal_Thermo ), log( max_basal_Thermo ) ) ) : exp( pars[counter++] );
        		basalTxps.push_back( basal );
    		}
	}
}

void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
		
    // write maxBindingWts
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }

    // write the interaction matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( coopMat( i, j ) ) {
                    double interaction = searchOption == CONSTRAINED ? infty_transform( log( factorIntMat( i, j ) ), log( min_interaction ), log( max_interaction ) ) : log( factorIntMat( i, j ) ); 
                   pars.push_back( interaction );
                }
            }
        }	       
    }
	
    // write the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC ) {
            double effect = searchOption == CONSTRAINED ? infty_transform( txpEffects[i], min_effect_Logistic, max_effect_Logistic ) : txpEffects[i];
            pars.push_back( effect );
        } else if ( modelOption == DIRECT ) {
            double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
            pars.push_back( effect );
        } else {
            if ( actIndicators[i] ) {
                double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
                pars.push_back( effect );
            }
        }
    }

    // write the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) {
                double repression = searchOption == CONSTRAINED ? infty_transform( log( repEffects[i] ), log( min_repression ), log( max_repression ) ) : log( repEffects[i] );
                pars.push_back( repression );
            }
        }
    }

    for( int i = 0; i < basalTxps.size(); i++ ){
    // write the basal transcription
    	if ( modelOption == LOGISTIC ) {
        	double basal = searchOption == CONSTRAINED ? infty_transform( basalTxps[ i ], min_basal_Logistic, max_basal_Logistic ) : basalTxps[ i ];
        	pars.push_back( basal );
    	} else {
        	double basal = searchOption == CONSTRAINED ? infty_transform( log( basalTxps[ i ] ), log( min_basal_Thermo ), log( max_basal_Thermo ) ) : log( basalTxps[ i ] );
        	pars.push_back( basal );
    	}
    }
}

void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const
{
//     os.setf( ios::fixed );
//     os.precision( 3 );
    
    // print the factor information
    os << "Motif" << "\t" << "max Binding Wights K[A]" << "\t" << "transcriptional factor";	
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t" << "repression Effects";
    os << endl;    

    for ( int i = 0; i < nFactors(); i++ ) {
        os << motifNames[i] << "\t \t " << maxBindingWts[i] << "\t \t " << txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t \t " << repEffects[i];
        os << endl;
    }

    // print the basal transcription
    os << "basal_transcription = " << basalTxps[ 0 ] << endl;
    for( int _i = 1; _i < basalTxps.size(); _i++ ){
    	os << basalTxps[ _i ] << endl;
    }
    
    // print the cooperative interactions
    os << "Cooperativity Factor:"  << endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << motifNames[i] << "\t" << motifNames[j] << "\t" << factorIntMat( i, j ) << endl;
        }
    }
}

int ExprPar::load( const string& file )
{
    // open the file
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl;	exit( 1 ); } 
    
    // read the factor information
    vector< string > motifNames( nFactors() );
    for ( int i = 0; i < nFactors(); i++ ) {
        fin >> motifNames[i] >> maxBindingWts[i] >> txpEffects[i]; 
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) fin >> repEffects[i];
    }

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    
    // read the basal transcription
    string symbol, eqSign, value;
    fin >> symbol >> eqSign >> value;
    if ( symbol != "basal_transcription" || eqSign != "=" ) return RET_ERROR;
    double basalTxp_val = atof( value.c_str() );
    basalTxps[ 0 ] =  basalTxp_val ;
    if( one_qbtm_per_crm ){
    	for( int _i = 1; _i < nSeqs; _i++ ){
    		fin >> value;
    		double basalTxp_val = atof( value.c_str() );
    		basalTxps[ _i ] = basalTxp_val;
    	}
    }
    
    // read the cooperative interactions
    string factor1, factor2;
    double coopVal;
    while ( fin >> factor1 >> factor2 >> coopVal ) {
        if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorIntMat( idx1, idx2 ) = coopVal;
        factorIntMat( idx2, idx1 ) = coopVal;
    }

    return 0;
}

void ExprPar::adjust()
{
    // adjust binding paramters
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) maxBindingWts[i] *= 2.0;
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) maxBindingWts[i] /= 2.0;
    }

    // adjust the interaction matrix 
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) ) { 
                factorIntMat( i, j ) *= 2.0; 
                factorIntMat( j, i ) = factorIntMat( i, j ); 
            }
            if ( factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
                factorIntMat( i, j ) /= 2.0;
                factorIntMat( j, i ) = factorIntMat( i, j ); 
            }
        }
    }
    
    // adjust transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC ) {
            if ( txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta ) txpEffects[i] /= 2.0;
            if ( txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta ) txpEffects[i] /= 2.0;
        } else {
            if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) txpEffects[i] *= 2.0;
            if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) txpEffects[i] /= 2.0;
        }
        
    }

    // adjust the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) ) repEffects[i] *= 10.0;
            if ( repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) ) repEffects[i] /= 2.0;
        }
    }

    // adjust the basl transcription
    for( int _i = 0; _i < basalTxps.size(); _i ++ )
    {
    	if ( modelOption == LOGISTIC ) {
        	if ( basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta ) basalTxps[ _i ] /= 2.0;
        	if ( basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta ) basalTxps[ _i ] *= 2.0;
    	} else {
        	if ( basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) ) basalTxps[ _i ] *= 2.0;
        	if ( basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) ) basalTxps[ _i ] /= 2.0;
    	}
    }
}

ModelType ExprPar::modelOption = DIRECT;
SearchType ExprPar::searchOption = UNCONSTRAINED;
int ExprPar::estBindingOption = 1;  // 1. estimate binding parameters; 0. not estimate binding parameters
 
double ExprPar::default_weight = 1.0;
double ExprPar::default_interaction = 1.0;
double ExprPar::default_effect_Logistic = 0.0;
double ExprPar::default_effect_Thermo = 1.0;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -5.0;
double ExprPar::default_basal_Thermo = 0.01;
double ExprPar::min_weight = 0.01;		
double ExprPar::max_weight = 100;//500;		
double ExprPar::min_interaction = 0.01;	
double ExprPar::max_interaction = 100;//500;
double ExprPar::min_effect_Logistic = -5;	
double ExprPar::max_effect_Logistic = 5;
// double ExprPar::min_effect_Direct = 0.01;
double ExprPar::min_effect_Thermo = 0.01;	
double ExprPar::max_effect_Thermo = 500;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Logistic = -9.0;	
double ExprPar::max_basal_Logistic = -1.0;
double ExprPar::min_basal_Thermo = 1.0E-4;	
double ExprPar::max_basal_Thermo = 0.1;
double ExprPar::delta = 0.0001;


bool ExprPar::one_qbtm_per_crm = true;
bool ExprFunc::one_qbtm_per_crm = true;

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, const ExprPar& _par, const vector< Sequence >& _seqs ) : motifs( _motifs ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ), par( _par ), seqs( _seqs )
{

    

    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
    bindingWts.clear();
    boundaries.clear();
	
    int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    int range = max(coopDistThr, repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] / ( 1.0 + bindingWts[i] );
        }
        double totalEffect = 0;
//         cout << "factor\toccupancy\ttxp_effect" << endl;
        for ( int i = 0; i < motifs.size(); i++ ) {
            double effect = par.txpEffects[i] * factorOcc[i];
            totalEffect += effect;
//             cout << i << "\t" << factorOcc[i] << "\t" << effect << endl;

            // length correction
//             totalEffect = totalEffect / (double)length;
        }
//         return par.expRatio * logistic( log( par.basalTxp ) + totalEffect );
        return logistic( par.basalTxps[ promoter_number ] + totalEffect );
    }

    #if TOSCA  
   //timeval start, end;
  //gettimeofday(&start, 0);
    double Z_off = 0;
    double Z_on = 0;
    compPartFunc_seq(Z_on, Z_off, seq_num, factorConcs);
  //gettimeofday(&end, 0);
  //cout <<end.tv_usec-start.tv_usec << endl;
    #else
    double Z_off = compPartFuncOff();
    double Z_on = compPartFuncOn();
    #endif //TOSCA


    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    double promoterOcc = efficiency * par.basalTxps[ promoter_number ] / ( 1.0 + efficiency * par.basalTxps[ promoter_number ] );
    return promoterOcc;
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num, std::ofstream& fout )
{
    bindingWts.clear();
    boundaries.clear();
	
    int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    int range = max( coopDistThr, repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] / ( 1.0 + bindingWts[i] );
        }
        double totalEffect = 0;
//         cout << "factor\toccupancy\ttxp_effect" << endl;
        for ( int i = 0; i < motifs.size(); i++ ) {
            double effect = par.txpEffects[i] * factorOcc[i];
            totalEffect += effect;
//             cout << i << "\t" << factorOcc[i] << "\t" << effect << endl;

            // length correction
//             totalEffect = totalEffect / (double)length;
        }
//         return par.expRatio * logistic( log( par.basalTxp ) + totalEffect );
        return logistic( par.basalTxps[ promoter_number ] + totalEffect );
    }

    // Thermodynamic models: Direct, Quenching, ChrMod_Unlimited and ChrMod_Limited
    // compute the partition functions
    #if TOSCA  
 // timeval start, end;
  //gettimeofday(&start, 0);
    double Z_off = 0;
    double Z_on = 0;
    compPartFunc_seq(Z_on, Z_off, seq_num, factorConcs);
  //gettimeofday(&end, 0);
  //cout <<end.tv_usec-start.tv_usec << endl;
    #else
    double Z_off = compPartFuncOff();
    double Z_on = compPartFuncOn();
    #endif //TOSCA

    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    double promoterOcc = efficiency * par.basalTxps[ promoter_number ] / ( 1.0 + efficiency * par.basalTxps[ promoter_number ] );
    fout << Z_on << "\t" << Z_off << "\t" << par.basalTxps [ promoter_number ] << endl;
    return promoterOcc;
}

ModelType ExprFunc::modelOption = DIRECT;

double ExprFunc::compPartFuncOff() const
{
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) return compPartFuncOffChrMod(); 
        
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
        Z[i] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }
    return Zt[n];
}

int ExprFunc::compPartFunc_seq(double &result_Z_on, double &result_Z_off, int _seq_num, const vector< double >& factorConcs) const
{

    // The distance index
    int dmax = max( coopDistThr, repressionDistThr );

    // The binder index
    int tmax = motifs.size();	
    vector< int > motif_length( tmax ); 
    for(  int t = 0; t < tmax; t++  ) { motif_length[t] = motifs[ t ].length(); } 
   
    // The Sequenz index
    int imax = *max_element(motif_length.begin(),motif_length.end());
    int n = seqs[_seq_num].size() - imax;

    // initialization of the partition sums with value 1
    vector<vector<vector<double> > > Z_off (tmax,vector<vector<double> >(dmax + 1,vector <double>(imax,1.0)));
    vector<vector<vector<double> > > Z_on  (tmax,vector<vector<double> >(dmax + 1,vector <double>(imax,1.0))); 

    // recurrence 
    int idx = 0;
    int idx_1 = 0;
    for ( int i = 1; i < n; i++ ) {

	idx = i % imax;
	idx_1 = (i-1) % imax;

	for ( int t = 0; t < tmax; t++ ) {
		   int idx_t = (idx - motif_length[t]) % imax;
		   
		   double sum_on = 0;
		   double sum_off = 0;
		   for (int d = 1; d < motif_length[t]; d++ ){
			Z_on[t][d][idx] = Z_on[t][d-1][idx_1];
			Z_off[t][d][idx] = Z_off[t][d-1][idx_1];
		   } 
		   for ( int d = motif_length[t]; d < dmax; d++ ) {
			Z_on[t][d][idx] = Z_on[t][d-1][idx_1];
			Z_off[t][d][idx] = Z_off[t][d-1][idx_1];
			for ( int t_alt = 0; t_alt < tmax; t_alt++ ){ 
				sum_on  += compFactorInt( t_alt, t, d ) * Z_on[t_alt][d - motif_length[t]][idx_t] ;
				sum_off += compFactorInt( t_alt, t, d ) * Z_off[t_alt][d - motif_length[t]][idx_t] ;
			}
		   }

      		   Sequence elem( seqs[_seq_num], i, motif_length[t] , true );
	           Sequence elem_comp( seqs[_seq_num], i, motif_length[t] , false );
        	   double binding_weight = exp( -motifs[ t ].energy( elem ) ) + exp( -motifs[ t ].energy( elem_comp ) );

	           Z_on[t][0][idx] = par.txpEffects[ t ] *  par.maxBindingWts[ t ] * factorConcs[t] * binding_weight * sum_on;
	           Z_off[t][0][idx] =  par.maxBindingWts[ t ] * factorConcs[t] * binding_weight * sum_off;
 
		   Z_on[t][dmax][idx] = Z_on[t][dmax][idx_1] + Z_on[t][dmax-1][idx_1];
		   Z_off[t][dmax][idx] = Z_off[t][dmax][idx_1] + Z_off[t][dmax-1][idx_1];
	}
    }
    

    result_Z_on = 0;
    result_Z_off = 0;
    for( int t = 0; t < tmax; t++ ) { 
	for(int d = 0; d < dmax+1; d++ ) {
		result_Z_on += Z_on[t][d][idx];
		result_Z_off += Z_off[t][d][idx];
	}        
    }
    return 1;
}


double ExprFunc::compPartFuncOffChrMod() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 ); 
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;
    
    // recurrence
    for ( int i = 1; i <= n; i++ ) {
        double sum = Zt[boundaries[i]]; 
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            int dist = sites[i].start - sites[j].start;
            
            // sum for Z0 
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j]; 

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] ) {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1; 
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n]; 
}

double ExprFunc::compPartFuncOn() const
{
    if ( modelOption == DIRECT ) return compPartFuncOnDirect();
    if ( modelOption == QUENCHING ) return compPartFuncOnQuenching();
    if ( modelOption == CHRMOD_UNLIMITED) return compPartFuncOnChrMod_Unlimited(); 
    if ( modelOption == CHRMOD_LIMITED ) return compPartFuncOnChrMod_Limited();
}

double ExprFunc::compPartFuncOnDirect() const
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
        Z[i] = bindingWts[ i ] * par.txpEffects[ sites[i].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }
     return Zt[n];
}

double ExprFunc::compPartFuncOnQuenching() const
{
    int n = sites.size() - 1;
    int N0 = maxContact;
    //vector< vector< double > > Z1( n + 1, N0 + 1 );
    ////vector< vector< double > > Z1( n + 1, vector < double > ( N0 + 1 ) );
    vector < vector < double > >Z1;
    for( int i = 0; i < n + 1; i++ ){
    	Z1.push_back( vector < double > ( N0 + 1, 0 ));
    }
    //vector< vector< double > > Z0( n + 1, N0 + 1 );
    ////vector< vector< double > > Z0( n + 1, vector < double > ( N0 + 1 ) );

	vector < vector <double> > Z0;
	for( int i = 0; i < n + 1; i++ ){
		Z0.push_back( vector < double > ( N0 + 1, 0 ) );
	}

    // k = 0
    for ( int i = 0; i <= n; i++ ) {
        double sum1 = 1, sum0 = 0;
        for ( int j = 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            bool R = testRepression( sites[j], sites[i] );
            double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1[j][0] + Z0[j][0] );
            sum1 += ( 1 - R )* term;
            sum0 += R * term;
        }
        Z1[i][0] = bindingWts[i] * sum1;
        Z0[i][0] = bindingWts[i] * sum0;
    }
    
    // k >= 1
    for ( int k = 1; k <= N0; k++ ) {
        for ( int i = 0; i <= n; i++ ) {
            if ( i < k ) {
                Z1[i][k] = 0; 
                Z1[i][k] = 0;
                continue;
            }
            double sum1 = 0, sum0 = 0;
            for ( int j = 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                bool R = testRepression( sites[j], sites[i] );
                double effect = actIndicators[sites[j].factorIdx] * ( 1 - testRepression( sites[i], sites[j] ) ) * Z1[j][k - 1] * par.txpEffects[sites[j].factorIdx];
                double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1[j][k] + Z0[j][k] + effect );
                sum1 += ( 1 - R )* term;
                sum0 += R * term;
            }
            Z1[i][k] = bindingWts[i] * sum1;
            Z0[i][k] = bindingWts[i] * sum0;
        }       
    }

//     for ( int i = 1; i <= n; i++ ) {
//         for ( int k = 0; k <= N0; k++ ) {
//             cout << "Z1(" << i << ", " << k << ") = " << Z1[i][k] << "\t";
//             cout << "Z0(" << i << ", " << k << ") = " << Z0[i][k] << endl;
//         }
//         cout << endl;
//     }
    
    // the partition function 
    double Z_on = 1;
    for ( int i = 1; i <= n; i++ ) {
        for ( int k = 0; k <= N0; k++ ) {
            double term = Z1[i][k] + Z0[i][k];
            Z_on += term;
        }	
        for ( int k = 0; k <= N0 - 1; k++ ) {
            Z_on += actIndicators[sites[i].factorIdx] * Z1[i][k] * par.txpEffects[sites[i].factorIdx];
        }
    }
    return Z_on;
}

double ExprFunc::compPartFuncOnChrMod_Unlimited() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 ); 
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;
    
    // recurrence
    for ( int i = 1; i <= n; i++ ) {
        double sum = Zt[boundaries[i]]; 
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            int dist = sites[i].start - sites[j].start;
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

            // sum for Z0
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j]; 

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] ) {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * par.txpEffects[ sites[i].factorIdx ] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1; 
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];     
}

double ExprFunc::compPartFuncOnChrMod_Limited() const
{
    int n = sites.size()- 1;

    // initialization
    int N0 = maxContact;
    ////vector< vector< double > > Z0( n + 1, N0 + 1 );
    ////vector< vector< double > > Z1( n + 1, N0 + 1 );
    ////vector< vector< double > > Zt( n + 1, N0 + 1 );
    
    //vector< vector< double > > Z0( n + 1, vector< double > ( N0 + 1 ) );
    vector < vector < double > > Z0;
    for( int i = 0; i < n + 1; i++ ){
    	Z0.push_back( vector < double >( N0 + 1, 0 ) );
    }
    
    //vector< vector< double > > Z1( n + 1, vector< double > ( N0 + 1 ) );
    vector < vector < double > > Z1;
    for( int i = 0; i < n + 1; i++ ){
    	Z1.push_back( vector < double >( N0 + 1, 0 ) );
    }
    //vector< vector< double > > Zt( n + 1, vector< double > ( N0 + 1 ) );
    vector < vector < double > > Zt;
    for( int i = 0; i < n + 1; i++ ){
    	Zt.push_back( vector < double >( N0 + 1, 0 ) );
    }
    Z0[0][0] = 0;
    Z1[0][0] = 0;
    Zt[0][0] = 1.0;
    for ( int k = 1; k <= N0; k++ ) {
        Z0[0][k] = 0;
        Z1[0][k] = 0;
        Zt[0][k] = 0;
    }

    // recurrence
    for ( int k = 0; k <= N0; k++ ) {
        for ( int i = 1; i <= n; i++ ) {
//             cout << "k = " << k << " i = " << i << endl;
            double sum0 = Zt[boundaries[i]][k], sum0A = k > 0 ? Zt[boundaries[i]][k-1] : 0, sum1 = sum0;
            for ( int j = boundaries[i] + 1; j < i; j++ ) {
                int dist = sites[i].start - sites[j].start;
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

                // sum for Z0
                sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j][k];
                sum0A += k > 0 ? compFactorInt( sites[i], sites[j] ) * Z0[j][k-1] : 0; 
                if ( dist > repressionDistThr ) {
                    sum0 += Z1[j][k]; 
                    sum0A += k > 0 ? Z1[j][k-1] : 0;
                }

                // sum for Z1
                if ( repIndicators[ sites[i].factorIdx ] ) {
                    sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j][k];
                    if ( dist > repressionDistThr ) sum1 += Z0[j][k];
                }
            }
            Z0[i][k] = bindingWts[i] * sum0;
            if ( actIndicators[sites[i].factorIdx] ) Z0[i][k] += k > 0 ? bindingWts[i] * par.txpEffects[sites[i].factorIdx] * sum0A : 0;
            if ( repIndicators[ sites[i].factorIdx ] ) Z1[i][k] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1; 
            else Z1[i][k] = 0;
            Zt[i][k] = Z0[i][k] + Z1[i][k] + Zt[i - 1][k];  
//             cout << "i = " << i << " k = " << k << " Z0 = " << Z0[i][k] << " Z1 = " << Z1[i][k] << " Zt = " << Zt[i][k] << endl;
        }
    }

    // the partition function
//     cout << "Zt[n] = " << Zt[n] << endl;
    return sum( Zt[n] );         
}

double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
// 	assert( !siteOverlap( a, b, motifs ) );
	
    double maxInt = par.factorIntMat( a.factorIdx, b.factorIdx );
    int dist = abs( a.start - b.start );
    assert( dist >= 0 );

    assert(  modelOption == DIRECT  );	// For now only Direct model implemented
    #if FactorIntFunc == 1 
    double spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 );
    #endif // FactorIntFunc

    #ifdef ORIENTATION
    double orientationTerm = ( a.strand == b.strand ) ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
    #else
    return spacingTerm;
    #endif //ORIENTATION
}

double ExprFunc::compFactorInt( int t_1, int t_2, int _dist  ) const
{

	
    double maxInt = par.factorIntMat( t_1, t_2 );
    int dist = _dist;
    assert( dist >= 0 );

    assert(  modelOption == DIRECT  );	// For now only Direct model implemented
    #if FactorIntFunc == 1 
    double spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 );
    #endif // FactorIntFunc

    #ifdef ORIENTATION
    double orientationTerm = ( a.strand == b.strand ) ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
    #else
    return spacingTerm;
    #endif //ORIENTATION
}


bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
// 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}

ExprPredictor::ExprPredictor( const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts, const vector< Sequence >& _seqs  ) : seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ) ,indicator_bool ( _indicator_bool ),  axis_start ( _axis_start ), axis_end( _axis_end ), axis_wts( _axis_wts ), seqs(_seqs)
{

    motifNames = _motifNames;
    coopMat = _coopMat;

    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
    assert( coopMat.isSquare() && coopMat.isSymmetric() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() );
    assert( maxContact > 0 );
    assert( repIndicators.size() == nFactors() );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors() );
    assert( repressionDistThr >= 0 );
	
	gene_crm_fout.open( "./output/gene_crm_fout.txt" );


    // set the model option for ExprPar and ExprFunc
    ExprPar::modelOption = modelOption;
    ExprFunc::modelOption = modelOption;

    // set the values of the parameter range according to the model option
    if ( modelOption != LOGISTIC && modelOption != DIRECT ) {
        ExprPar::min_effect_Thermo = 0.99;
        ExprPar::min_interaction = 0.99;
    }

    // set the option of parameter estimation
    ExprPar::estBindingOption = estBindingOption;
}

double ExprPredictor::objFunc( const ExprPar& par ) 
{
    if ( objOption == SSE ) return compRMSE( par );	
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par ); 
    if ( objOption == NORM_CORR ) return -compNormCorr( par ); 
}

int ExprPredictor::train( const ExprPar& par_init )
{
    par_model = par_init;
    par_curr = par_init;
    obj_model = objFunc( par_model );   

    signal(SIGINT, catch_signal);


    cout << "*** Diagnostic printing BEFORE adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_curr );
    cout << endl;
    cout << "Objective function value: " << obj_model << endl;
    cout << "*******************************************" << endl << endl;

   if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
    
    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
//    cout << "Objective function value: " << obj_model << endl;
    cout << "*******************************************" << endl << endl;



    if ( nAlternations == 0 ) return 0;
    
    // alternate between two different methods
    ExprPar par_result;
    double obj_result;
    for ( int i = 0; i < nAlternations; i++ ) {
	cout << "Minimisation step " << i+1 << " of " << nAlternations << endl; 
	cout << "Simplex minimisation: " << endl; 
	simplex_minimize( par_result, obj_result );
	cout << endl;

	// save result
        par_model = par_result; 
	save_param();	

//         par_model.adjust();
	cout << "Gradient minimisation: " << endl; 
        gradient_minimize( par_result, obj_result );
	cout << endl;

	// save result
        par_model = par_result;
	save_param();

//         par_model.adjust();
    }
	
    // commit the parameters and the value of the objective function
    par_model = par_result; 
    obj_model = obj_result;
		
    return 0;	
}

int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
/*
	//for random starts: 
	ExprPar par_rand_start = par_init;
	randSamplePar( rng, par_rand_start );
	train( par_rand_start );*/
    // training using the initial values
    train( par_init );
    
    cout << "Initial training:\tParameters = "; printPar( par_model ); 
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    // training with random starts
	ExprPar par_best = par_model;
	double obj_best = obj_model;
	for ( int i = 0; i < nRandStarts; i++ ) {
        	ExprPar par_curr = par_init; 
		randSamplePar( rng, par_curr ); 
		train( par_curr );
        	cout << "Random start " << i + 1 << ":\tParameters = "; printPar( par_model );
        	cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
		if ( obj_model < obj_best ) {
			par_best = par_model;
			obj_best = obj_model;	
		}
	}    

    // training using the best parameters so far
    if ( nRandStarts ) train( par_best ); 
    cout << "Final training:\tParameters = "; printPar( par_model ); 
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl; 
    

	gene_crm_fout.close();

    return 0;
}

int ExprPredictor::train()
{	
    // random number generator
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
    
    // training using the default initial values with random starts
    ExprPar par_default( nFactors(), nSeqs() );
    train( par_default, rng ); 
    
    return 0;	
}

int ExprPredictor::predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
    targetExprs.clear();
    
    // create site representation of the target sequence
//     SiteVec targetSites;
//     SeqAnnotator ann( motifs, energyThrs );	
//     ann.annot( targetSeq, targetSites );
            
    // predict the expression
    ExprFunc* func = createExprFunc( par_model );	
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr( targetSites, targetSeqLength, concs, seq_num );
        targetExprs.push_back( predicted );
    }
    
    return 0; 
}

// double ExprPredictor::test( const vector< Sequence >& testSeqs, const Matrix& testExprData, int perfOption ) const
// {
// 	assert( perfOption == 0 );	
// 	assert( testExprData.nRows() == testSeqs.size() && testExprData.nCols() == nConds() ); 
// 	
// 	// make predictions
// 	Matrix predicted( testSeqs.size(), nConds() );
// 	for ( int i = 0; i < testSeqs.size(); i++ ) {
//         vector< double > targetExprs;
//         predict( targetSeqs[i], targetExprs );
//         predicted.setRow( i, targetExprs );
// 	}
// 	
// 	// RMSE between predictions and observations
// 	if ( perfOption == 0 ) {
// 		vector< double > corrs; 
// 		for ( int i = 0; i < nExps; i++ ) {
// 			corrs.push_back( correlation( predicted[ i ], testExprData[ i ] ) );
// 		}
// 		
// 		return mean( corrs );	
// 	}
// }

ModelType ExprPredictor::modelOption = CHRMOD_LIMITED;
int ExprPredictor::estBindingOption = 1;    // 1. estimate binding parameters; 0. not estimate binding parameters
ObjType ExprPredictor::objOption = OBJECTIVE_FUNCTION;

double ExprPredictor::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
{
    vector< int > shifts; 
    for ( int s = -maxShift; s <= maxShift; s++ ) {
        shifts.push_back( s ); 
    }

    vector< double > cov; 
    vector< double > corr; 
    cross_corr( x, y, shifts, cov, corr ); 
    double result = 0, weightSum = 0; 
//     result = corr[maxShift]; 
    result = *max_element( corr.begin(), corr.end() );
//     for ( int i = 0; i < shifts.size(); i++ ) {
//         double weight = pow( shiftPenalty, abs( shifts[i] ) ); 
//         weightSum += weight; 
//         result += weight * corr[i]; 
//     }
//     result /= weightSum; 

    return result; 
}

int ExprPredictor::maxShift = 5; 
double ExprPredictor::shiftPenalty = 0.8; 

int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 5;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
double ExprPredictor::min_delta_f_NormCorr = 1.0E-8;
int ExprPredictor::nSimplexIters = 200;
int ExprPredictor::nGradientIters = 50;
bool ExprPredictor::one_qbtm_per_crm = true;

// Initialise static members as empty
ExprPar ExprPredictor::par_curr;
IntMatrix ExprPredictor::coopMat = IntMatrix();
vector <string> ExprPredictor::motifNames = vector <string>();

int ExprPredictor::randSamplePar( const gsl_rng* rng, ExprPar& par ) const
{
    // sample binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double rand_weight = exp( gsl_ran_flat( rng, log( ExprPar::min_weight ), log( ExprPar::max_weight ) ) ); 
            par.maxBindingWts[i] = rand_weight;
        }        
    }

    // sample the interaction matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interaction ), log( ExprPar::max_interaction ) ) );
                if ( coopMat( i, j ) ) par.factorIntMat( i, j ) = rand_interaction;
            }
        }
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = i + 1; j < nFactors(); j++ ) {
                par.factorIntMat( i, j ) = par.factorIntMat( j, i );
            }
        }       
    }

    // sample the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC ) {
            double rand_effect = gsl_ran_flat( rng, ExprPar::min_effect_Logistic, ExprPar::max_effect_Logistic );
            par.txpEffects[i] = rand_effect;
        } else if ( modelOption == DIRECT ) {
            double rand_effect = exp( gsl_ran_flat( rng, log( ExprPar::min_effect_Thermo ), log( ExprPar::max_effect_Thermo ) ) );
            par.txpEffects[i] = rand_effect;
        } else {
            if ( actIndicators[i] ) {
                double rand_effect = exp( gsl_ran_flat( rng, log( ExprPar::min_effect_Thermo ), log( ExprPar::max_effect_Thermo ) ) );
                par.txpEffects[i] = rand_effect;
            }
        }
    }

    // sample the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) {
                double rand_repression = exp( gsl_ran_flat( rng, log( ExprPar::min_repression ), log( ExprPar::max_repression ) ) );
                par.repEffects[i] = rand_repression;
            }
        }
    }
    
    // sample the basal transcription
    double rand_basal;
    for( int _i = 0; _i < par.basalTxps.size(); _i ++ ){
    	if ( modelOption == LOGISTIC ) 
        	rand_basal = gsl_ran_flat( rng, ExprPar::min_basal_Logistic, ExprPar::max_basal_Logistic );
    	else
        	rand_basal = exp( gsl_ran_flat( rng, log( ExprPar::min_basal_Thermo ), log( ExprPar::max_basal_Thermo ) ) );
    	par.basalTxps[ _i ] = rand_basal;
    }
    return 0;
}

bool ExprPredictor::testPar( const ExprPar& par ) const
{
    // test binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
                return false; 
        }        
    }

    // test the interaction matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( par.factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    // test the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC ) {
            if ( par.txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta || par.txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta ) 
                return false;
        } else if ( modelOption == DIRECT ) {
            if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        } else {
            if ( actIndicators[i] ) {
                if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    // test the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) {
                if ( par.repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) || par.repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }
    
    // test the basal transcription
    for( int _i = 0; _i < par.basalTxps.size(); _i++ ){
    	if ( modelOption == LOGISTIC ) {
        	if ( par.basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta || par.basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta ) 
            	return false;
    	} else {
        	if ( par.basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) || par.basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) )
            	return false;
    	}
    }
    return true;    
}

void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
//    cout.precision( 3 ); 
//     cout.width( 8 ); 
    
    // print binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            cout << par.maxBindingWts[i] << "\t"; 
        }        
    }

    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) cout << par.factorIntMat( i, j ) << "\t";
        }
    }

    // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC || modelOption == DIRECT ) cout << par.txpEffects[i] << "\t";
        else {
            if ( actIndicators[i] ) cout << par.txpEffects[i] << "\t";
        }
    }

    // print the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) cout << par.repEffects[i] << "\t"; 
        }
    }
    
    // print the basal transcriptions
    for( int _i = 0; _i < par.basalTxps.size(); _i ++ ){
    	cout << par.basalTxps[ _i ] << "\t"; 
    }
    cout << endl;
}

ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{	
    return new ExprFunc( motifs, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, par, seqs );
}

int indices_of_crm_in_gene[] = {
	5, 11, 17, 23, 29
};

double ExprPredictor::compRMSE( const ExprPar& par ) 
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );

    // error of each sequence
    double squaredErr = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
		for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
			if( i == indices_of_crm_in_gene[ _i ] ){
				gene_crm_fout << i << "\t" << j << "\t";
            			predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i, gene_crm_fout );
				break;
			}
		}	
		if ( predicted < 0 ){
            		predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i );
		}
		
            
            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        double beta;
        squaredErr += least_square( predictedExprs, observedExprs, beta );
	

    }	

    double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
    return rmse;
}

double ExprPredictor::compAvgCorr( const ExprPar& par ) 
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
            
    // Pearson correlation of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
		for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
			if( i == indices_of_crm_in_gene[ _i ] ){
				gene_crm_fout << i << "\t" << j << "\t";
            			predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i, gene_crm_fout );
				break;
			}
		}	
		if ( predicted < 0 ){
            		predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i );
		}
		
            
            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += corr( predictedExprs, observedExprs ); 
//         cout << "Sequence " << i << "\t" << corr( predictedExprs, observedExprs ) << endl;
    }	

    return totalSim / nSeqs();
}

double ExprPredictor::compAvgCrossCorr( const ExprPar& par ) 
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
            
    // cross correlation similarity of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
		for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
			if( i == indices_of_crm_in_gene[ _i ] ){
				gene_crm_fout << i << "\t" << j << "\t";
            			predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i, gene_crm_fout );
				break;
			}
		}	
		
		if ( predicted < 0 ){
            		predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i );
		}
            
            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs ); 
    }	

    return totalSim / nSeqs();
}

double ExprPredictor::compNormCorr( const ExprPar& par ) 
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
            
    // Normalised correlation of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
		for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
			if( i == indices_of_crm_in_gene[ _i ] ){
				gene_crm_fout << i << "\t" << j << "\t";
            			predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i, gene_crm_fout );
				break;
			}
		}	
		if ( predicted < 0 ){
            		predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i );
		}
		
            
            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += norm_corr( predictedExprs, observedExprs ); 
      //  cout << "Sequence " << i << "\t" << norm_corr( predictedExprs, observedExprs ) << endl;
    }	

    return totalSim / nSeqs();
}


int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result ) 
{
// 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector < double > pars;
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
        
	//Hassan start:
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
	
	pars.clear();
	pars = free_pars;

	//Hassan end
    // set the objective function
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    // set the initial values to be searched
    gsl_vector* x = vector2gsl( pars );
//     for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( x, i, v[ i ] );

	// CHECK POINT: evaluate gsl_obj_f() function
// 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
    // choose the method of optimization and set its parameters
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );
            
    // create the minimizer
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );
    gsl_multimin_fminimizer_set( s, &my_func, x, ss );
    
    // iteration
    size_t iter = 0;
    int status;
    double size;	
    do {
        double f_prev = iter ? s->fval : 1.0E6;     // the function starts with some very large number
        
        iter++;
        status = gsl_multimin_fminimizer_iterate( s );
     
        // check for error
        if ( status ) break;
	//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}


        par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators, nSeqs() );

	//Hassan end
	// check if the current values of parameters are valid
        //the following line should be uncommented if you remove all the changes by Hassan
	//ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
        
	
	if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;
        
        // check for stopping condition
//         double f_curr = s->fval;
//         double delta_f = abs( f_curr - f_prev ); 
//         if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
//         if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
//         if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
        
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size, 1e-1 ); 
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

        // print the current parameter and function values
	#ifdef SHORT_OUTPUT
	if(iter % SHORT_OUTPUT == 0){
      		printf( "\r %i \t f() = %8.5f", iter, s->fval );	
		fflush(stdout);
	}	
	#else
	cout << "======================================" << endl;
	cout << "======================================" << endl;
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
	cout << "======================================" << endl;
	par_curr.print( cout, motifNames, coopMat);
	cout << "======================================" << endl;
	cout << "======================================" << endl << endl;
	#endif // SHORT_OUTPUT
    } while ( status == GSL_CONTINUE && iter < nSimplexIters );

    // get the results
//     vector< double > expv;
//     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );	
	//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}

        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators, nSeqs() );
	//Hassan end
	//uncomment the following line if you remove all the changes by Hassan
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->fval;
    
    // free the minimizer
    gsl_vector_free( x );
    gsl_vector_free( ss );
    gsl_multimin_fminimizer_free( s );	
    
    return 0;
}

int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result ) 
{
// 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
            
	//Hassan start:
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			//cout << "testing 1: " << pars[ index ] << endl;
			free_pars.push_back( pars[ index ]);
		}
		else{
			//cout << "testing 2: " << pars[ index ] << endl;
			fix_pars.push_back( pars[ index ] );
		}
	}
	
	pars.clear();
	pars = free_pars;
	//Hassan end
    // set the objective function and its gradient
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    // set the initial values to be searched 
    gsl_vector* x = vector2gsl( pars ); 

	// CHECK POINT: evaluate gsl_obj_f() function
// 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
	// choose the method of optimization and set its parameters
// 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2; // Chose Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm 
		
    // create the minimizer
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = 0.02, tol = 0.1;
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
            
    // iteration
    size_t iter = 0;
    int status;
    do {
        double f_prev = iter ? s->f : 1.0E6;     // the function starts with some very large number	
        
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
// 	    if ( prev_f - curr_f < 0.001 ) break;
         
        // check for error
        if ( status ) break;

        // check if the current values of parameters are valid
	//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}

        par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators, nSeqs() );
	
	//Hassan end
        //ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;
     
        // check for stopping condition
        double f_curr = s->f;
        double delta_f = abs( f_curr - f_prev ); 
        if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
        if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
        if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
        if ( objOption == NORM_CORR && delta_f < min_delta_f_NormCorr ) break;
        
        status = gsl_multimin_test_gradient( s->gradient, 5e-4 );
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

         // print the current parameter and function values
	#ifdef SHORT_OUTPUT
	if(iter % SHORT_OUTPUT == 0){
        	printf( "\r %i \t f() = %8.5f", iter, s->f );
		fflush(stdout);
	}
	#else
	cout << "========================================" << endl;
	cout << "========================================" << endl;
        cout << iter << "\t";
        //printPar( par_curr );           //Redundant information 
        printf( "\tf() = %8.5f\n", s->f );
	cout << "========================================" << endl;
	/*vector <string> motifNames;
	motifNames.push_back( "bcd" );
	motifNames.push_back( "cad" );
	motifNames.push_back( "gt" );
	motifNames.push_back( "hb" );
	motifNames.push_back( "kni" );
	motifNames.push_back( "Kr" );
	motifNames.push_back( "nub" );*/
	par_curr.print( cout, motifNames, coopMat);
	cout << "========================================" << endl;
	cout << "========================================" << endl << endl;
	#endif // SHORT_OUTPUT
    } while ( status == GSL_CONTINUE && iter < nGradientIters );

    // get the results
//     vector< double > expv;
//     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );	
	//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}

        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators, nSeqs() );
	//Hassan end
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->f;
    
    // free the minimizer
    gsl_vector_free( x );    
    gsl_multimin_fdfminimizer_free( s );
    
    return 0;
}

// function to save parameters to file
int ExprPredictor::save_param()
{
	ofstream fparam_sm( "param.save" );
	par_curr.print( fparam_sm, motifNames, getCoopMat() );
        fparam_sm.close();
	return 0;
}


// Signal handler functions 
void ExprPredictor::catch_signal(int sig_num)
{
	std::cout << std::endl;
	std::cout << "Exit signal thrown!" << std::endl;
	save_param();
	std::cout << "Parameters are saved!" << std::endl;
	signal(SIGINT, SIG_DFL);
}



double gsl_obj_f( const gsl_vector* v, void* params )
{ 
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}



    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(), predictor -> nSeqs() );
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
            
    // call the ExprPredictor object to evaluate the objective function 
    double obj = predictor->objFunc( par );	
    return obj;
}

/*
int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result ) const
{
// 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
            
    // set the objective function
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    // set the initial values to be searched
    gsl_vector* x = vector2gsl( pars );
//     for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( x, i, v[ i ] );

	// CHECK POINT: evaluate gsl_obj_f() function
// 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
    // choose the method of optimization and set its parameters
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );
            
    // create the minimizer
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );
    gsl_multimin_fminimizer_set( s, &my_func, x, ss );
    
    // iteration
    size_t iter = 0;
    int status;
    double size;	
    do {
        double f_prev = iter ? s->fval : 1.0E6;     // the function starts with some very large number
        
        iter++;
        status = gsl_multimin_fminimizer_iterate( s );
     
        // check for error
        if ( status ) break;

        // check if the current values of parameters are valid
        ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;
        
        // check for stopping condition
//         double f_curr = s->fval;
//         double delta_f = abs( f_curr - f_prev ); 
//         if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
//         if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
//         if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
        
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size, 1e-1 ); 
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

        // print the current parameter and function values
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
    } while ( status == GSL_CONTINUE && iter < nSimplexIters );

    // get the results
//     vector< double > expv;
//     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );	
    par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->fval;
    
    // free the minimizer
    gsl_vector_free( x );
    gsl_vector_free( ss );
    gsl_multimin_fminimizer_free( s );	
    
    return 0;
}

int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result ) const
{
// 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
            
    // set the objective function and its gradient
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    // set the initial values to be searched 
    gsl_vector* x = vector2gsl( pars ); 

	// CHECK POINT: evaluate gsl_obj_f() function
// 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
	// choose the method of optimization and set its parameters
// 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs;
		
    // create the minimizer
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = 0.02, tol = 0.1;
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
            
    // iteration
    size_t iter = 0;
    int status;
    do {
        double f_prev = iter ? s->f : 1.0E6;     // the function starts with some very large number	
        
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
// 	    if ( prev_f - curr_f < 0.001 ) break;
         
        // check for error
        if ( status ) break;

        // check if the current values of parameters are valid
        ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;
     
        // check for stopping condition
        double f_curr = s->f;
        double delta_f = abs( f_curr - f_prev ); 
        if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
        if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
        if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
        
        status = gsl_multimin_test_gradient( s->gradient, 5e-4 );
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

         // print the current parameter and function values
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f\n", s->f );
    } while ( status == GSL_CONTINUE && iter < nGradientIters );

    // get the results
//     vector< double > expv;
//     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );	
    par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->f;
    
    // free the minimizer
    gsl_vector_free( x );    
    gsl_multimin_fdfminimizer_free( s );
    
    return 0;
}

double gsl_obj_f( const gsl_vector* v, void* params )
{ 
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp = gsl2vector(v);
	cout << "*********samee start**************" << endl;
	for(int temp_i = 0; temp_i < temp.size(); temp_i ++)
		cout << temp[temp_i] << endl;
	cout << "*********samee end****************" << endl;


    ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
            
    // call the ExprPredictor object to evaluate the objective function 
    double obj = predictor->objFunc( par );	
    return obj;
}
*/
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-3;
    numeric_deriv( grad, gsl_obj_f, v, params, step );	
}

void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params ); 
    gsl_obj_df( v, params, grad );		
}
