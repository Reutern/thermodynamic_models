#include </usr/local/include/gsl/gsl_math.h>
#include </usr/local/include/gsl/gsl_multimin.h>

#include "siman.h"
#include "ExprPredictor.h"
#include "param.h"
#include <sys/time.h>
#include <omp.h>
#include <unistd.h>
#include <set>

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
    if ( toupperStr( objOptionStr ) == "SSE_scale" ) return SSE_SCALE;
    if ( toupperStr( objOptionStr ) == "SSE_V" ) return SSE_V;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "NORM_CORR" ) return NORM_CORR;
    if ( toupperStr( objOptionStr ) == "NORM_CORR_V" ) return NORM_CORR_V;
    if ( toupperStr( objOptionStr ) == "PGP" ) return PGP;

    cerr << "objOptionStr is not a valid option of objective function" << endl; 
    exit(1);
}

string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == SSE_SCALE ) return "SSE_scale";
    if ( objOption == SSE_V ) return "SSE_V";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == NORM_CORR ) return "Norm_Corr";
    if ( objOption == NORM_CORR_V ) return "Norm_Corr_V";
    if ( objOption == PGP ) return "PGP";

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

    acc_scale = ExprPar::default_acc_scale;
    acc_base = ExprPar::default_acc_base;

}
	
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, int _nSeqs, double _acc_scale, double _acc_base ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), nSeqs( _nSeqs  ), acc_scale(_acc_scale), acc_base(_acc_base)
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
            //double effect = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++],  min_effect_Thermo ,  max_effect_Thermo  ) :  pars[counter++] ;
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

    #if ACCESSIBILITY
    // Write the accessibility parameter
    acc_scale = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_acc_scale ), log( max_acc_scale ) ) ) : exp( pars[counter++] );
    acc_base = searchOption == CONSTRAINED ?  exp( inverse_infty_transform( pars[counter++], log( min_acc_base  ), log( max_acc_base  ) ) ) : exp( pars[counter++] );
    #endif // ACCESSIBILITY

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
	    //double effect = searchOption == CONSTRAINED ? infty_transform( txpEffects[i] , min_effect_Thermo , max_effect_Thermo  ) :  txpEffects[i] ;
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
    
    #if ACCESSIBILITY
    // Write the accessibility parameter
    double scale = searchOption == CONSTRAINED ?  infty_transform( log( acc_scale ), log( min_acc_scale ), log( max_acc_scale ) ) : log( acc_scale );
    pars.push_back( scale );

    double base = searchOption == CONSTRAINED ? infty_transform( log( acc_base ), log( min_acc_base ), log( max_acc_base ) ) : log( acc_base );
    pars.push_back( base );
  
    #endif // ACCESSIBILITY

}

void ExprPar::print( ostream& os, const vector< string >& motifNames, const vector< string >& seqNames, const IntMatrix& coopMat ) const
{
//     os.setf( ios::fixed );
//     os.precision( 3 );
    
    // print the factor information
    //os << "Motif" << "\t" << "max Binding Wights K[A]" << "\t" << "transcriptional factor";	
    //if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t" << "repression Effects";
    //os << endl;    

    for ( int i = 0; i < nFactors(); i++ ) {
        os << motifNames[i] << "\t \t " << maxBindingWts[i] << "\t \t " << txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t \t " << repEffects[i];
        os << endl;
    }

    #if ACCESSIBILITY
    os << "Accessibility = " << acc_scale << "\t " << acc_base << endl;
    #endif // ACCESSIBILITY

    // print the basal transcription
    if (one_qbtm_per_crm == false){
        os << "basal_transcription = " << basalTxps[ 0 ] << endl;
    }
    else{
	for( int _i = 0; _i < basalTxps.size(); _i++ ){
	    os << seqNames[_i] << "\t" << basalTxps[ _i ] << endl;
	}	
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
    
    string symbol, eqSign, value;

    #if ACCESSIBILITY
    string value1, value2;
    fin >> symbol >> eqSign >> value1 >> value2;
    if (symbol != "Accessibility" || eqSign != "=") return RET_ERROR;
    acc_scale = atof( value1.c_str() );
    acc_base = atof( value2.c_str() );
    #endif // ACCESSIBILITY

    // read the basal transcription
    if( one_qbtm_per_crm ){
    	for( int _i = 0; _i < nSeqs; _i++ ){
    		fin >> symbol >> value;
    		double basalTxp_val = atof( value.c_str() );
    		basalTxps[ _i ] = basalTxp_val;
    	}
    }
    else{
	fin >> symbol >> eqSign >> value;
	if ( symbol != "basal_transcription" || eqSign != "=" ) return RET_ERROR;
	double basalTxp_val = atof( value.c_str() );
	basalTxps[ 0 ] =  basalTxp_val ;
    }

    fin >> symbol >> eqSign;
    if (symbol != "Cooperativity" || eqSign != "Factor:") return RET_ERROR;

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

    fin.close();
    return fin ? RET_ERROR : 0;
}

int ExprPar::load( const string& file, const vector <string>& seqNames )
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

    string symbol, eqSign, value;

    #if ACCESSIBILITY
    string value1, value2;
    fin >> symbol >> eqSign >> value1 >> value2;
    if (symbol != "Accessibility" || eqSign != "=") return RET_ERROR;
    acc_scale = atof( value1.c_str() );
    acc_base = atof( value2.c_str() );
    #endif // ACCESSIBILITY

    // read the basal transcription
    if( one_qbtm_per_crm ){
    	if( seqNames.size() != nSeqs ) return RET_ERROR;
    	for( int _i = 0; _i < nSeqs; _i++ ){
    		fin >> symbol >> value;
    		double basalTxp_val = atof( value.c_str() );
		// allocate basalTxps to seqence names
		for( int _l = 0; _l < nSeqs; _l++ ){
    			if (symbol == seqNames[_l] ){
				basalTxps[ _l ] = basalTxp_val;
				break;
			}
		}
    	}
    }
    else{
	fin >> symbol >> eqSign >> value;
	if ( symbol != "basal_transcription" || eqSign != "=" ) return RET_ERROR;
	double basalTxp_val = atof( value.c_str() );
	basalTxps[ 0 ] =  basalTxp_val ;
    }

    fin >> symbol >> eqSign;
    if (symbol != "Cooperativity" || eqSign != "Factor:") return RET_ERROR;

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

    fin.close();
    return fin ? RET_ERROR : 0;
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

    #if ACCESSIBILITY
    if ( acc_scale < ExprPar::min_acc_scale * ( 1.0 + ExprPar::delta ) ) acc_scale *= 2.0;
    if ( acc_scale > ExprPar::max_acc_scale * ( 1.0 - ExprPar::delta ) ) acc_scale /= 2.0;
    if ( acc_base < ExprPar::min_acc_base * ( 1.0 + ExprPar::delta ) ) acc_base *= 2.0;
    if ( acc_base > ExprPar::max_acc_base * ( 1.0 - ExprPar::delta ) ) acc_base /= 2.0;
    #endif // ACCESSIBILITY


}

void ExprPar::constrain_parameters()
{
    // constrain binding paramters
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) maxBindingWts[i] = ExprPar::min_weight;
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) maxBindingWts[i] = ExprPar::max_weight;
    }

    // constrain the interaction matrix 
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) ) { 
                factorIntMat( i, j ) = ExprPar::min_interaction; 
                factorIntMat( j, i ) = factorIntMat( i, j ); 
            }
            if ( factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
                factorIntMat( i, j ) = ExprPar::max_interaction;
                factorIntMat( j, i ) = factorIntMat( i, j ); 
            }
        }
    }
    
    // constrain transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC ) {
            if ( txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta ) txpEffects[i] = ExprPar::min_effect_Logistic;
            if ( txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta ) txpEffects[i] = ExprPar::max_effect_Logistic;
        } else {
            if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) txpEffects[i] = ExprPar::min_effect_Thermo;
            if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) txpEffects[i] = ExprPar::max_effect_Thermo;
        }
        
    }

    // constrain the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) ) repEffects[i] = ExprPar::min_repression;
            if ( repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) ) repEffects[i] = ExprPar::max_repression;
        }
    }

    // constrain the basl transcription
    for( int _i = 0; _i < basalTxps.size(); _i ++ )
    {
    	if ( modelOption == LOGISTIC ) {
        	if ( basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta ) basalTxps[ _i ] = ExprPar::min_basal_Logistic;
        	if ( basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta ) basalTxps[ _i ] = ExprPar::max_basal_Logistic;
    	} else {
        	if ( basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) ) basalTxps[ _i ] = ExprPar::min_basal_Thermo;
        	if ( basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) ) basalTxps[ _i ] = ExprPar::max_basal_Thermo;
    	}
    }

    #if ACCESSIBILITY
    if ( acc_scale < ExprPar::min_acc_scale * ( 1.0 + ExprPar::delta ) ) acc_scale = ExprPar::min_acc_scale;
    if ( acc_scale > ExprPar::max_acc_scale * ( 1.0 - ExprPar::delta ) ) acc_scale = ExprPar::max_acc_scale;
    if ( acc_base < ExprPar::min_acc_base * ( 1.0 + ExprPar::delta ) ) acc_base = ExprPar::min_acc_base;
    if ( acc_base > ExprPar::max_acc_base * ( 1.0 - ExprPar::delta ) ) acc_base = ExprPar::max_acc_base;
    #endif // ACCESSIBILITY

}

ModelType ExprPar::modelOption = DIRECT;
SearchType ExprPar::searchOption = CONSTRAINED;
int ExprPar::estBindingOption = 1;  // 1. estimate binding parameters; 0. not estimate binding parameters
 
// Parameter limits
double ExprPar::default_acc_scale = 1.0;
double ExprPar::default_acc_base = 1;
double ExprPar::default_weight = 1.0;
double ExprPar::default_interaction = 1.0;
double ExprPar::default_effect_Logistic = 0.0;
double ExprPar::default_effect_Thermo = 1.0;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -5.0;
double ExprPar::default_basal_Thermo = 0.01;
double ExprPar::min_acc_scale = 0.001;
double ExprPar::min_acc_base = 0.001;
double ExprPar::max_acc_scale = 10.0;
double ExprPar::max_acc_base = 3;
double ExprPar::min_weight = 0.0001;		
double ExprPar::max_weight = 5000;//500;		
double ExprPar::min_interaction = 0.001;	
double ExprPar::max_interaction = 500;
double ExprPar::min_effect_Logistic = -5;	
double ExprPar::max_effect_Logistic = 5;
// double ExprPar::min_effect_Direct = 0.01;
double ExprPar::min_effect_Thermo = 0.0001;	
double ExprPar::max_effect_Thermo = 1000;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Logistic = -9.0;	
double ExprPar::max_basal_Logistic = -1.0;
double ExprPar::min_basal_Thermo = 1.0E-5;	
double ExprPar::max_basal_Thermo = 1.0;
double ExprPar::delta = 0.0001;


bool ExprPar::one_qbtm_per_crm = ONE_QBTM;
bool ExprFunc::one_qbtm_per_crm = ONE_QBTM;

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, const ExprPar& _par, const vector< Sequence >& _seqs ) : motifs( _motifs ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ), par( _par ), seqs( _seqs )
{
    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors ); 
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
}


// Returns the efficiency Z_ON/Z_OFF
double ExprFunc::predictExpr_scalefree(int length, const vector< double >& factorConcs, int seq_num )
{

   int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

     #if TOSCA  
  // timeval start, end;
  //gettimeofday(&start, 0);
    double Z_off = 0;
    double Z_on = 0;
    compPartFunc_seq_interfactor(Z_on, Z_off, seq_num, factorConcs);
  //gettimeofday(&end, 0);
  //cout <<end.tv_usec-start.tv_usec << endl;
    #else

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] * factorConcs[sites[i].factorIdx] / ( 1.0 + bindingWts[i] * factorConcs[sites[i].factorIdx] );
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

    double Z_off = 0;
    double Z_on = 0;

 
    Z_off = compPartFuncOff(factorConcs);
    Z_on = compPartFuncOn(factorConcs);
    
    #endif //TOSCA


    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    return efficiency;
}

double ExprFunc::predictExpr( int length, const vector< double >& factorConcs, int seq_num )
{

   int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

     #if TOSCA  
  // timeval start, end;
  //gettimeofday(&start, 0);
    double Z_off = 0;
    double Z_on = 0;
    compPartFunc_seq_interfactor(Z_on, Z_off, seq_num, factorConcs);
  //gettimeofday(&end, 0);
  //cout <<end.tv_usec-start.tv_usec << endl;
    #else

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] * factorConcs[sites[i].factorIdx] / ( 1.0 + bindingWts[i] * factorConcs[sites[i].factorIdx] );
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

    double Z_off = 0;
    double Z_on = 0;

    Z_off = compPartFuncOff(factorConcs);
    Z_on = compPartFuncOn(factorConcs);
    #endif //TOSCA


    // Test if partition functions are in the range of double and set output accordingly
    if(std::isnan(Z_off) || std::isinf(Z_off)) return 0;
    if(std::isnan(Z_on) || std::isinf(Z_on)) return 1;

    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    double promoterOcc = efficiency * par.basalTxps[ promoter_number ] / ( 1.0 + efficiency * par.basalTxps[ promoter_number ] );
    return promoterOcc;
}

double ExprFunc::predictExpr( int length, const vector< double >& factorConcs, int seq_num, std::ofstream& fout )
{

    int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

     #if TOSCA  
 // timeval start, end;
  //gettimeofday(&start, 0);
    double Z_off = 0;
    double Z_on = 0;
    compPartFunc_seq_interfactor(Z_on, Z_off, seq_num, factorConcs);
  //gettimeofday(&end, 0);
  //cout <<end.tv_usec-start.tv_usec << endl;
  
   #else

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] * factorConcs[sites[i].factorIdx] / ( 1.0 + bindingWts[i] * factorConcs[sites[i].factorIdx] );
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

    double Z_off = 0;
    double Z_on = 0;


    Z_off = compPartFuncOff(factorConcs);
    Z_on = compPartFuncOn(factorConcs);
    #endif //TOSCA

    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    double promoterOcc = efficiency * par.basalTxps[ promoter_number ] / ( 1.0 + efficiency * par.basalTxps[ promoter_number ] );
    fout << Z_on << "\t" << Z_off << "\t" << par.basalTxps [ promoter_number ] << endl;
    return promoterOcc;
}

double ExprFunc::predictExpr_scanning_mode( int length, const vector< double >& factorConcs, int seq_num )
{

    int promoter_number = seq_num;
	if( !one_qbtm_per_crm )
		promoter_number = 0;	// Only one promoter strength

    int n = sites.size() - 1;

    vector<double> p_bound;

    compProb_scanning_mode( factorConcs, p_bound );

    // Parse the sequence into fragments
    const int nFragments = max(static_cast<int>( (length - 300)/ 50 ),1);
    vector<int> sites_begin (nFragments+1, 1);    
    vector<int> sites_end (nFragments+1, n);    


    int fragment_counter_start = 0;
    int fragment_counter_end = 0;
    for ( int idx_site = 1; idx_site <= n; idx_site++ ) {
	int position = sites[idx_site].start;
	if( ( position > 50 * fragment_counter_start ) &&  ( fragment_counter_start < nFragments ) ) {
		sites_begin[fragment_counter_start] = idx_site; 
		fragment_counter_start++ ;	
	} 
	if( position > 50 * fragment_counter_end + 300) {
		sites_end[fragment_counter_end] = idx_site; 
		fragment_counter_end++ ;	
	}
    }

    double promoterOcc = 0.0;
    for(int idx_fragagment = 0; idx_fragagment < max(fragment_counter_end,1); idx_fragagment++){
	double weight_complex = 1;
	for(int idx = sites_begin[idx_fragagment]; idx < sites_end[idx_fragagment]; idx++){
		weight_complex *= exp( - p_bound[idx] * par.txpEffects[ sites[idx].factorIdx ] );
	}
	double p_tmp = weight_complex * par.basalTxps[ promoter_number ] / ( 1.0 + par.basalTxps[ promoter_number ] * weight_complex );
	promoterOcc = promoterOcc * (1-p_tmp) + p_tmp; 
    }

    
    return promoterOcc;

}

void ExprFunc::compProb_scanning_mode(const vector< double >& factorConcs, vector< double >& p_bound)
{

    // initialization
    int n = sites.size();

    p_bound.clear();
    p_bound.resize(n);

    for ( int i = 1; i < n; i++ ) {
	double weight_competition = 0;
	double weight_cooperativity = 1;
        for ( int j = max(i - 30,0); j < min(i + 30 , n); j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) weight_competition += bindingWts[ j ] * factorConcs[sites[ j ].factorIdx];
		else weight_cooperativity += (compFactorInt( sites[ i ], sites[ j ] )-1) *  bindingWts[ j ] * factorConcs[sites[ j ].factorIdx];
                }
	double weight = bindingWts[ i ] * factorConcs[sites[ i ].factorIdx];
        p_bound[i] = weight * weight_cooperativity / (1.0 + weight * weight_cooperativity + weight_competition);
    }
 
}


ModelType ExprFunc::modelOption = DIRECT;

double ExprFunc::compPartFuncOff(const vector< double >& factorConcs) const
{
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) return compPartFuncOffChrMod( factorConcs ); 

    int n = sites.size() - 1;

    // initialization
    vector <double> Z (n + 1, 0.0);
    Z[0] = 1.0;
    vector <double> Zt (n + 1, 0.0);
    Zt[0] = 1.0;

    // recurrence 
    for ( int i = 1; i <= n; i++ ) {
	double sum = Zt[boundaries[i]]; 
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];	
        }
        Z[i] =   bindingWts[ i ] * factorConcs[sites[ i ].factorIdx] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }       
    return Zt[n];
}

// TOSCA with interfactor cooperativity (very slow)
int ExprFunc::compPartFunc_seq_interfactor(double &result_Z_on, double &result_Z_off, int _seq_num, const vector< double >& factorConcs) const
{

    // The distance index
    int dmax = max( coopDistThr, repressionDistThr ) + 1;

    // The binder index
    int tmax = motifs.size();	
    vector< int > motif_length( tmax ); 
    for(  int t = 0; t < tmax; t++  ) { motif_length[t] = motifs[ t ].length(); } 
   
    // The Sequenz index
    int n = seqs[_seq_num].size();

//    vector<vector<vector<double> > > Z_off (n + 1,vector <vector <double> > ( tmax , vector <double>(dmax + 1,0.0)));
//    vector<vector<vector<double> > >Z_on  (n + 1,vector <vector <double> > ( tmax , vector <double>(dmax + 1,0.0))); 
    double* Z_off = new double[n  * tmax * dmax + tmax * dmax + dmax]();
    double* Z_on = new double[n  * tmax * dmax + tmax * dmax + dmax]();

    // recurrence 
    for (int i = 1; i < n; i++ ) {

//	#pragma omp parallel for
	for ( int t = 0; t < tmax; t++ ) {
		   int idx_t = (i - motif_length[t]);
		   if(idx_t < 1) continue;	   

		   double sum_on = 1.0;
		   double sum_off = 1.0;
		   for (int d = 1; d < dmax; d++ ){
			Z_on[i * tmax * dmax + t * dmax + d ] = Z_on[(i - 1) * tmax * dmax + t * dmax + d - 1];
			Z_off[i * tmax * dmax + t * dmax + d] = Z_off[(i - 1) * tmax * dmax + t * dmax + d - 1];
		   } 

		   for (int t_alt = 0; t_alt < tmax; t_alt++ ){ 
			double int_factor = min(par.factorIntMat(t,t_alt), 1.0);
				for ( int d = motif_length[t]; d < dmax+1; d++ ) {
					sum_on  += int_factor * Z_on[idx_t * tmax * dmax + t_alt * dmax + d] ;
					sum_off += int_factor * Z_off[idx_t * tmax * dmax + t_alt * dmax + d] ;
				}
		   }

      		   Sequence elem( seqs[_seq_num], idx_t, motif_length[t] , true );
	           Sequence elem_comp( elem, 0, motif_length[t] , false );
		   double e1 =  exp( -motifs[ t ].energy( elem ) );
		   double e2 =  exp( -motifs[ t ].energy( elem_comp ) );
        	   double binding_weight = par.maxBindingWts[ t ] * factorConcs[t] * ( e1+e2 );

	           Z_on[i * tmax * dmax + t * dmax ] = par.txpEffects[ t ] * binding_weight * sum_on;
	           Z_off[i * tmax * dmax + t * dmax ] =  binding_weight * sum_off;
 
		   Z_on[i * tmax * dmax + t * dmax + dmax] = Z_on[(i-1) * tmax * dmax + t * dmax + dmax] + Z_on[(i-1) * tmax * dmax + t * dmax + dmax - 1];
		   Z_off[i * tmax * dmax + t * dmax + dmax] = Z_off[(i-1) * tmax * dmax + t * dmax + dmax] + Z_off[(i-1) * tmax * dmax + t * dmax + dmax - 1];
	}
    }
    
    
    result_Z_on = 1;
    result_Z_off = 1;
    for( int t = 0; t < tmax; t++ ) { 
	for(int d = 0; d < dmax+1; d++ ) {
		result_Z_on += Z_on[n * tmax * dmax + t * dmax + d];
		result_Z_off += Z_off[n * tmax * dmax + t * dmax + d];
	}        
    }
 
    delete Z_off;
    delete Z_on;

    return 1;
}

// TOSCA without interfactor cooperativity (fast version)
int ExprFunc::compPartFunc_seq(double &result_Z_on, double &result_Z_off, int _seq_num, const vector< double >& factorConcs) const
{

    // The distance index
    int dmax = max( coopDistThr, repressionDistThr );

    // The binder index
    int tmax = motifs.size();	
    int n = seqs[_seq_num].size();


    result_Z_on = 0;
    result_Z_off = 0;

    // recurrence 
    int idx = 0;

    for ( int t = 0; t < tmax; t++ ) {
   // initialization of the partition sums with value 1
//    vector<vector<double> > Z_off (dmax + 1,vector <double>(n,1.0));
//    vector<vector<double> > Z_on  (dmax + 1,vector <double>(n,1.0)); 

    double Z_off [dmax + 1][n];
    double Z_on [dmax + 1][n];

    for(int d = 0; d < dmax + 1; d++){
	Z_on[d][0] = 0;
	Z_off[d][0] = 0;
    }

    int motif_length_tmp =  motifs[ t ].length();
   
    for ( int i = 1; i < n - motif_length_tmp; i++ ) {
	idx = i ;


		   
		   double sum_on = 1;
		   double sum_off = 1;
		   for (int d = 1; d < motif_length_tmp && d < i ; d++ ){
			Z_on[d][idx] = Z_on[d-1][idx-1];
			Z_off[d][idx] = Z_off[d-1][idx-1];
		   } 
		   for ( int d = motif_length_tmp; d < dmax && d < i; d++ ) {
			Z_on[d][idx] = Z_on[d-1][idx-1];
			Z_off[d][idx] = Z_off[d-1][idx-1];

			sum_on  += compFactorInt( t, t, d ) * Z_on[d - motif_length_tmp][idx - motif_length_tmp] ;
			sum_off += compFactorInt( t, t, d ) * Z_off[d - motif_length_tmp][idx - motif_length_tmp] ;
			
		   }

      		   Sequence elem( seqs[_seq_num], i, motif_length_tmp , true );
	           Sequence elem_comp( seqs[_seq_num], i, motif_length_tmp , false );
        	   double binding_weight = par.maxBindingWts[ t ] * factorConcs[t] * (exp( -motifs[ t ].energy( elem ) ) + exp( -motifs[ t ].energy( elem_comp ) ));

	           Z_on[0][idx] = par.txpEffects[ t ] * binding_weight * sum_on;
	           Z_off[0][idx] =  binding_weight * sum_off;
 
		   Z_on[dmax][idx] = Z_on[dmax][idx-1] + Z_on[dmax-1][idx-1];
		   Z_off[dmax][idx] = Z_off[dmax][idx-1] + Z_off[dmax-1][idx-1];
	}
	for(int d = 0; d < dmax + 1; d++ ) {
		result_Z_on += Z_on[d][idx];
		result_Z_off += Z_off[d][idx];
	}        
    }
    
    return 1;
}

double ExprFunc::compPartFuncOffChrMod(const vector< double >& factorConcs) const
{
    int n = sites.size()- 1;

    // initialization
    double Z0 [n + 1];
    Z0[0] = 1.0;
    double Z1 [n + 1];
    Z1[0] = 1.0;
    double Zt [n + 1];
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
        Z0[i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * par.repEffects[ sites[i].factorIdx ] * sum1; 
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n]; 
}

double ExprFunc::compPartFuncOn(const vector< double >& factorConcs) const
{	
    if ( modelOption == DIRECT ) return compPartFuncOnDirect(factorConcs);
    if ( modelOption == QUENCHING ) return compPartFuncOnQuenching(factorConcs);
    if ( modelOption == CHRMOD_UNLIMITED) return compPartFuncOnChrMod_Unlimited(factorConcs); 
    if ( modelOption == CHRMOD_LIMITED ) return compPartFuncOnChrMod_Limited(factorConcs);
}

double ExprFunc::compPartFuncOnDirect(const vector< double >& factorConcs) const
{
   int n = sites.size() - 1;
    
    // initialization
    vector <double> Z (n + 1, 0.0);
    Z[0] = 1.0;
    vector <double> Zt (n + 1, 0.0);
    Zt[0] = 1.0;
	

    // recurrence 
    for ( int i = 1; i <= n; i++ ) {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];	
        }
        Z[i] =  bindingWts[ i ] * factorConcs[sites[ i ].factorIdx] * par.txpEffects[ sites[i].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }
     return Zt[n];
}

double ExprFunc::compPartFuncOnQuenching(const vector< double >& factorConcs) const
{
    int n = sites.size() - 1;
    int N0 = maxContact;

    vector < vector <double> > Z0 (N0 + 1, vector <double> (n + 1, 0.0));
    vector < vector <double> > Z1 (N0 + 1, vector <double> (n + 1, 0.0));

    // k = 0
    for ( int i = 0; i <= n; i++ ) {
        double sum1 = 1, sum0 = 0;
        for ( int j = 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            bool R = testRepression( sites[j], sites[i] );
            double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1[0][j] + Z0[0][j] );
            sum1 += ( 1 - R )* term;
            sum0 += R * term;
        }
        Z1[0][i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * sum1;
        Z0[0][i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * sum0;
    }
    
    // k >= 1
    for ( int k = 1; k <= N0; k++ ) {
        for ( int i = 0; i <= n; i++ ) {
            if ( i < k ) {
                Z1[k][i] = 0; 
                Z1[k][i] = 0;
                continue;
            }
            double sum1 = 0, sum0 = 0;
            for ( int j = 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                bool R = testRepression( sites[j], sites[i] );
                double effect = actIndicators[sites[j].factorIdx] * ( 1 - testRepression( sites[i], sites[j] ) ) * Z1[k - 1][j] * par.txpEffects[sites[j].factorIdx];
                double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1[k][j] + Z0[k][j] + effect );
                sum1 += ( 1 - R )* term;
                sum0 += R * term;
            }
            Z1[k][i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * sum1;
            Z0[k][i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * sum0;
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
            double term = Z1[k][i] + Z0[k][i];
            Z_on += term;
        }	
        for ( int k = 0; k <= N0 - 1; k++ ) {
            Z_on += actIndicators[sites[i].factorIdx] * Z1[k][i] * par.txpEffects[sites[i].factorIdx];
        }
    }
    return Z_on;
}

double ExprFunc::compPartFuncOnChrMod_Unlimited(const vector< double >& factorConcs) const
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
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * factorConcs[sites[ i ].factorIdx] * par.repEffects[ sites[i].factorIdx ] * sum1; 
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];     
}

double ExprFunc::compPartFuncOnChrMod_Limited(const vector< double >& factorConcs) const
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
            if ( actIndicators[sites[i].factorIdx] ) Z0[i][k] += k > 0 ? bindingWts[i] * factorConcs[sites[ i ].factorIdx] * par.txpEffects[sites[i].factorIdx] * sum0A : 0;
            if ( repIndicators[ sites[i].factorIdx ] ) Z1[i][k] = bindingWts[i] * factorConcs[sites[ i ].factorIdx]  * par.repEffects[ sites[i].factorIdx ] * sum1; 
            else Z1[i][k] = 0;
            Zt[i][k] = Z0[i][k] + Z1[i][k] + Zt[i - 1][k];  
//             cout << "i = " << i << " k = " << k << " Z0 = " << Z0[i][k] << " Z1 = " << Z1[i][k] << " Zt = " << Zt[i][k] << endl;
        }
    }

    // the partition function
//     cout << "Zt[n] = " << Zt[n] << endl;
    return sum( Zt[n] );         
}

// The interaction function for direct model
double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{

// 	assert( !siteOverlap( a, b, motifs ) );
    //if(a.factorIdx != b.factorIdx)	return 1.0;	// Only TF of the same type interact 	
    double maxInt = par.factorIntMat( a.factorIdx, b.factorIdx );
    unsigned dist = abs( a.start - b.start );
    assert( dist >= 0 );

    //assert(  modelOption == DIRECT  );	// For now only Direct model implemented
    #if FactorIntFunc
    double spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 ); // Range 1.0 <-> (maxint + 1.0)
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

// The interaction function for TOSCA
double ExprFunc::compFactorInt( int t_1, int t_2, int _dist  ) const
{

    double maxInt = par.factorIntMat( t_1, t_2 );
    int dist = _dist;
    assert( dist >= 0 );

    //assert(  modelOption == DIRECT  );	// For now only Direct model implemented
    #if FactorIntFunc
    double spacingTerm = ( dist < coopDistThr ? maxInt *  (1 - float(dist/coopDistThr) ) : 0.0 );	// Range 0 <-> maxint 
    #else
    double spacingTerm = ( dist < coopDistThr ? maxInt : 0.0 );
    #endif // FactorIntFunc

    #if ORIENTATION
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

ExprPredictor::ExprPredictor( const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector <string>& _seqNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts, const vector< Sequence >& _seqs  ) : seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ) ,indicator_bool ( _indicator_bool ), axis_start ( _axis_start ), axis_end( _axis_end ), axis_wts( _axis_wts ), seqs(_seqs)
{

    motifNames = _motifNames;
    coopMat = _coopMat;
    seqNames = _seqNames;

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
    if ( objOption == SSE ||  objOption == NORM_CORR || objOption == PGP) return comp_SSE_NormCorr_PGP( par );	
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par ); 
}

int ExprPredictor::train( const ExprPar& par_init )
{   
    par_model = par_init;	// Initialise the model parameter
    par_curr = par_init;	// The working parameter, which get saved in case of an emergancy
    signal(SIGINT, catch_signal);

/*  
    cout << "*** Diagnostic printing BEFORE adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_curr );
    cout << endl;
    cout << "Objective function value: " << obj_model << endl;
    cout << "*******************************************" << endl << endl;

   if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
*/
    obj_model = objFunc( par_model ); 
/*    
    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << obj_model << endl;
    cout << "*******************************************" << endl << endl;
*/

    if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ) { par_model.constrain_parameters(); 
								       par_model.adjust(); }
    if ( nAlternations == 0 ) return 0;
    
    // alternate between two different methods
    ExprPar par_result;
    double obj_result;

    for ( int i = 0; i < nAlternations; i++ ) {
	cout << "Minimisation step " << i+1 << " of " << nAlternations << endl; 
	//objOption = SSE;
	cout << "Simplex minimisation: " << endl; 
        simplex_minimize( par_result, obj_result );
	//simulated_annealing( par_result, obj_result );
	cout << endl;
	//save result
        //par_model = par_result; 
	//save_param();	



	//objOption = NORM_CORR;
	//cout << "Simplex minimisation Norm_Corr: " << endl; 
        //simplex_minimize( par_result, obj_result );	
	//cout << endl;

	// save result
        par_model = par_result;
        par_model.adjust();
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
    targetExprs.resize( nConds() );
    // create site representation of the target sequence
//     SiteVec targetSites;
//     SeqAnnotator ann( motifs, energyThrs );	
//     ann.annot( targetSeq, targetSites );
            
    // predict the expression

    ExprFunc* func = createExprFunc( par_model );
    func->set_sites(targetSites);
    int n = targetSites.size();
    vector< double > _bindingWts (n,0.0);
    vector< int > _boundaries (n,0);


    // Determin the boundaries for func
    _boundaries[0] = 0;
    int range = max(coopDistThr, repressionDistThr );
    for ( int k = 1; k < n; k++ ) {
    	int l; 
	for ( l=0; l < k; l++ ) {
	    if ( ( targetSites[k].start - targetSites[l].start ) <= range ) {break;} 
	}
    _boundaries[k] = l ;
    }	
    func->set_boundaries(_boundaries);
    
    // compute the Boltzman weights of binding for all sites for func
    _bindingWts[0] = 1.0;
    for ( int k = 1; k < n; k++ ) {
	double access_tmp = 1.0;
	#if ACCESSIBILITY
	access_tmp = exp( - par_model.acc_scale * ( 1 - targetSites[k].accessibility) );
	#endif //ACCESSIBILITY
        _bindingWts[k] = access_tmp * par_model.maxBindingWts[ targetSites[k].factorIdx ] * targetSites[k].wtRatio ;	

    }
    func->set_bindingWts(_bindingWts); 
	
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr( targetSeqLength, concs, seq_num );	
        targetExprs[j] = predicted ;
    }
    
    delete func;

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

int ExprPredictor::nAlternations = 6;
int ExprPredictor::nRandStarts = 5;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
double ExprPredictor::min_delta_f_NormCorr = 1.0E-8;
double ExprPredictor::min_delta_f_PGP = 1.0E-8;
int ExprPredictor::nSimplexIters = 20000;
int ExprPredictor::nGradientIters = 50;
bool ExprPredictor::one_qbtm_per_crm = ONE_QBTM;

// Initialise static members as empty
ExprPar ExprPredictor::par_curr;
IntMatrix ExprPredictor::coopMat = IntMatrix();
vector <string> ExprPredictor::motifNames = vector <string>();
vector <string> ExprPredictor::seqNames = vector <string>();

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

    // print Accessibility parameters
	cout << par.acc_scale << "\t" << par.acc_base << "\t";

    cout << endl;
}

ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{
    return new ExprFunc( motifs, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, par, seqs );
}


int indices_of_crm_in_gene[] = {
	5, 11, 17, 23, 29
};

double ExprPredictor::comp_SSE_NormCorr_PGP( const ExprPar& par ) 
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
    // error of each sequence
    double squaredErr = 0;
    double correlation = 0;
    double pgp_score = 0;
    #if TIMER 
    timeval start, end;
    gettimeofday(&start, 0);
    #endif // TIMER

    for ( int i = 0; i < nSeqs(); i++ ) {

        vector< double > predictedExprs (nConds(), -1);
        vector< double > observedExprs (nConds(), 1);
        vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
        

	
        // Initiate the sites for func
	func->set_sites(seqSites[ i ]);

	int n = seqSites[i].size();
	vector< double > _bindingWts (n,0.0);
	vector< int > _boundaries (n,0);
	
	// Determin the boundaries for func
	_boundaries[0] = 0;
	int range = max(coopDistThr, repressionDistThr );
        for ( int k = 1; k < n; k++ ) {
		int l;	       	 
		for ( l=0; l < k; l++ ) {
		    if ( ( seqSites[i][k].start - seqSites[i][l].start ) <= range ) {break;} 
		}
	        _boundaries[k] = l;
	}	
        func->set_boundaries(_boundaries);
	// compute the Boltzman weights of binding for all sites for func
        _bindingWts[0] = 1.0;
        for ( int k = 1; k < n; k++ ) {
		double access_tmp = 1.0;
		#if ACCESSIBILITY
		access_tmp = exp(  - par.acc_scale * (1- seqSites[i][k].accessibility) ) ;
		#endif //ACCESSIBILITY
		_bindingWts[k] = access_tmp * par.maxBindingWts[ seqSites[i][k].factorIdx ] * seqSites[i][k].wtRatio ;	
        }
	func->set_bindingWts(_bindingWts); 


       #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < nConds(); j++ ) {		
         	concs[j] = factorExprData.getCol( j );
      		predictedExprs[j] = func->predictExpr( seqLengths[ i ], concs[j], i );	//_scanning_mode
                // observed expression for the i-th sequence at the j-th condition
                observedExprs[j] = exprData( i, j );
	        //weight += observedExprs[j];				// The weight of an enhancer is proportional to its expression width
        }

	//weight = ( 100 - weight ) / 100.0;			// Transformation of the weight from 0 - 100 to 1 - 0

        double beta = 1.0;
        squaredErr +=  least_square( predictedExprs, observedExprs, beta );
	correlation  +=  corr( predictedExprs, observedExprs ); 
        pgp_score += pgp( predictedExprs, observedExprs, beta );
    }	

    #if TIMER 
    gettimeofday(&end, 0);
    cout << "Time " << (end.tv_sec-start.tv_sec)+1e-6*(end.tv_usec-start.tv_usec) << endl;
    #endif // TIMER

    delete func;
    obj_norm_corr = correlation / nSeqs();
    obj_sse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
    obj_pgp = pgp_score / nSeqs();
    if (objOption == SSE)	return obj_sse;
    if (objOption == NORM_CORR)	return -obj_norm_corr;
    if (objOption == PGP)	return -obj_pgp;
    return 0;
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

	int n = seqSites[i].size();
	vector< double > _bindingWts (n,0.0);
	vector< int > _boundaries (n,0.0);
	
	// Determin the boundaries for func
	_boundaries[0] = 0;
	int range = max(coopDistThr, repressionDistThr );
        for ( int k = 1; k < n; k++ ) {
	       	int l; 
		for ( l = k - 1; l >= 1; l-- ) {
		    if ( ( seqSites[i][k].start - seqSites[i][l].start ) > range ) break; 
		}
	    _boundaries[k] = l ;
	}	
        func->set_boundaries(_boundaries);
    
	// compute the Boltzman weights of binding for all sites for func
        _bindingWts[0] = 1.0;
        for ( int k = 1; k < n; k++ ) {
		double access_tmp = 1.0;
		#if ACCESSIBILITY
		access_tmp = exp(  - par_model.acc_scale * (1- seqSites[i][k].accessibility) ) ;
		#endif //ACCESSIBILITY
		_bindingWts[k] = access_tmp * par_model.maxBindingWts[ seqSites[i][k].factorIdx ] * seqSites[i][k].wtRatio ;
        }
	func->set_bindingWts(_bindingWts); 

       // #pragma omp parallel for schedule(dynamic)
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
		for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
			if( i == indices_of_crm_in_gene[ _i ] ){
				gene_crm_fout << i << "\t" << j << "\t";
            			//predicted = func->predictExpr_scanning_mode( seqLengths[ i ], concs, i, gene_crm_fout );
            			predicted = func->predictExpr( seqLengths[ i ], concs, i );
				break;
			}
		}	
		if ( predicted < 0 ){
            		predicted = func->predictExpr( seqLengths[ i ], concs, i );
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

        vector< double > predictedExprs (nConds(), -1);
        vector< double > observedExprs (nConds(), 1);
        vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
        

	
        // Initiate the sites for func
	func->set_sites(seqSites[ i ]);

	int n = seqSites[i].size();
	vector< double > _bindingWts (n,0.0);
	vector< int > _boundaries (n,0);
	
	// Determin the boundaries for func
	_boundaries[0] = 0;
	int range = max(coopDistThr, repressionDistThr );
        for ( int k = 1; k < n; k++ ) {
		int l;	       	 
		for ( l=0; l < k; l++ ) {
		    if ( ( seqSites[i][k].start - seqSites[i][l].start ) <= range ) {break;} 
		}
	        _boundaries[k] = l;
	}	
        func->set_boundaries(_boundaries);
	// compute the Boltzman weights of binding for all sites for func
        _bindingWts[0] = 1.0;
        for ( int k = 1; k < n; k++ ) {
 		double access_tmp = 1.0;
		#if ACCESSIBILITY
		access_tmp =  exp( - par_model.acc_scale * ( 1 - seqSites[i][k].accessibility ) ) ;
		#endif //ACCESSIBILITY
		_bindingWts[k] = access_tmp * par_model.maxBindingWts[ seqSites[i][k].factorIdx ] * seqSites[i][k].wtRatio ;
        }
	func->set_bindingWts(_bindingWts); 

       // #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < nConds(); j++ ) {

			
            	concs[j] = factorExprData.getCol( j );
		//for( int _i = 0; _i < sizeof ( indices_of_crm_in_gene ) / sizeof ( int ); _i++ ){
		//	if( i == indices_of_crm_in_gene[ _i ] ){
		//		gene_crm_fout << i << "\t" << j << "\t";
           	//		predictedExprs[j] = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs[j], i, gene_crm_fout );
		//		break;
		//	}
		//}	
		//if ( predictedExprs[j] < 0 ){
            		predictedExprs[j] = func->predictExpr( seqLengths[ i ], concs[j], i );
		//}

            // observed expression for the i-th sequence at the j-th condition
            observedExprs[j] = exprData( i, j );
	    //weight += observedExprs[j];				// The weight of an enhancer is proportional to its expression width
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs ); 
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
    gsl_vector_set_all( ss, 3 );
            
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
	#if FILE_OUTPUT
	if(iter % 1000 == 0){
    		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_norm_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	#ifdef SHORT_OUTPUT
	if(iter % SHORT_OUTPUT == 0){
   		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_norm_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	cout << "======================================" << endl;
	cout << "======================================" << endl;
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
	cout << "======================================" << endl;
	par_curr.print( cout, motifNames, seqNames, coopMat);
	cout << "======================================" << endl;
	cout << "======================================" << endl << endl;
	#endif // SHORT_OUTPUT
	#endif // FILE_OUTPUT
	// Save parameters
	if( iter % 1000 == 0 ) save_param();


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
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr; // Chose Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm 
		
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
        	printf( "\r %zu \t f() = %8.5f", iter, s->f );
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
	par_curr.print( cout, motifNames, seqNames, coopMat);
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

int ExprPredictor::simulated_annealing( ExprPar& par_result, double& obj_result ) 
{

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_siman_params_t siman_params = {100, 10000, 0.1, 0.2, 0.1, 1.1, 2.0e-7};

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
    // set the initial values to be searched 
    gsl_vector* x = vector2gsl( pars ); 
    gsl_siman_solve(r, x, (void*)this, gsl_obj_f, siman_stepper,siman_print, NULL, NULL, NULL, sizeof(*x), siman_params);

    // free the minimiser
    gsl_vector_free( x );    
    gsl_rng_free (r);
    
    return 0;
}

void siman_print(gsl_vector * xp)
{
}

void siman_stepper(const gsl_rng * r, gsl_vector* v, double step_size)
{
    int idx = int(gsl_rng_uniform(r) * v->size);
    double u = gsl_rng_uniform(r) * 2 * step_size - step_size + gsl_vector_get( v, idx );
    gsl_vector_set( v , idx , u ) ; 
}

// function to save parameters to file
int ExprPredictor::save_param()
{
	ofstream fparam_sm( "param.save" );
	par_curr.print( fparam_sm, motifNames, seqNames, getCoopMat() );
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
    predictor->par_curr = par;
    // call the ExprPredictor object to evaluate the objective function 
    double obj = predictor->objFunc( par );
    return obj;
}

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
