#include </usr/local/include/gsl/gsl_math.h>
#include </usr/local/include/gsl/gsl_multimin.h>

#include "cmaes.h"
#include "ExprPredictor.h"
#include "param.h"
#include <sys/time.h>
#include <omp.h>
#include <unistd.h>
#include <set>

ExprPar::ExprPar( int _nFactors, int _nSeqs ) : factorIntMat(), factorSynMat(), factorSkewMat()
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );  

    factorSynMat.setDimensions( _nFactors, _nFactors );
    factorSynMat.setAll( ExprPar::default_synergy );   

    factorSkewMat.setDimensions( _nFactors, _nFactors );
    factorSkewMat.setAll( ExprPar::default_synergy );       

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
}
	
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const Matrix& _factorSynMat, const Matrix& _factorSkewMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, int _nSeqs, double _acc_scale) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), factorSynMat( _factorSynMat ), factorSkewMat( _factorSkewMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), nSeqs( _nSeqs  ), acc_scale(_acc_scale)
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() ); 	
    if ( !factorSynMat.isEmpty() ) assert( factorSynMat.nRows() == maxBindingWts.size() && factorSynMat.isSquare() ); 	
    assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
    	if ( one_qbtm_per_crm ){
    		assert( basalTxps.size() == nSeqs );
	}
	else{
		assert( basalTxps.size() == 1 );
	}
}

ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const IntMatrix& SynMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs ) : factorIntMat(), factorSynMat(), factorSkewMat()
{	
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( SynMat.isSquare() && SynMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
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

    // set the synergy matrix
    factorSynMat.setDimensions( _nFactors, _nFactors );
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( SynMat( i, j ) ) {
                double synergy = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_synergy ), log( max_synergy ) ) ) : exp( pars[counter++] );
                factorSynMat( i, j ) = synergy;  
            }
            else factorSynMat( i, j ) = ExprPar::default_synergy;
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {
            factorSynMat( i, j ) = factorSynMat( j, i );
        }
    }       

    // set the skew matrix
    factorSkewMat.setDimensions( _nFactors, _nFactors );
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
                double skew = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_skew, max_skew) : pars[counter++];
                factorSkewMat( i, j ) = skew;  
            }
            else factorSkewMat( i, j ) = ExprPar::default_skew;
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {
            factorSkewMat( i, j ) = factorSkewMat( j, i );
        }
    }  

    // set the transcriptional effects
    for ( int i = 0; i < _nFactors; i++ ) {
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
	} else {
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
    #endif // ACCESSIBILITY

}

double ExprPar::parameter_L2_norm() const
{
	double L2_weights = 0;
	double L2_effects = 0;

    for ( int i = 0; i < nFactors(); i++ ) {

		L2_weights += pow( maxBindingWts[ i ] / max_weight,2 )/nFactors();
		if(txpEffects[i] >= 1)
			L2_effects += pow( (txpEffects[i] -1) / ( max_effect_Thermo -1 ), 2) / nFactors();
		else
			L2_effects += pow( min_effect_Thermo / txpEffects[i], 2) / nFactors();
	}


	return L2_weights + L2_effects;
}


double ExprPar::parameter_L2_norm_interactions() const
{
	double L2_norm_coop = 0;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
			#if NEGATIVE_COOP
			if(factorIntMat( i, j ) >= 1)
				L2_norm_coop += pow((factorIntMat(i, j) - 1)/ max_interaction, 2);
			else
				L2_norm_coop +=  pow(min_interaction / factorIntMat( i, j ), 2);
			#else
			L2_norm_coop += pow(factorIntMat(i, j)/ max_interaction, 2);
			#endif // NEGATIVE_COOP
        }	
	}

	L2_norm_coop /= nFactors();

	double L2_norm_syn = 0;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
			if(factorSynMat( i, j ) >= 1)
				L2_norm_syn += pow((factorSynMat(i, j) - 1)/ max_synergy, 2);
			else
				L2_norm_syn +=  pow(min_synergy / factorSynMat( i, j ), 2);
		}
	}

	L2_norm_syn /= nFactors();

	return L2_norm_coop + L2_norm_syn;
}

double ExprPar::parameter_L1_norm_skew() const
{
	double L1_norm_skew = 0;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
			L1_norm_skew += fabs(factorSkewMat(i, j));
		}
	}
    return L1_norm_skew;
}

double ExprPar::parameter_L1_norm() const
{
	double L1_weights = 0;
	double L1_effects = 0;

    for ( int i = 0; i < nFactors(); i++ ) {

		L1_weights +=  maxBindingWts[ i ] / max_weight / nFactors();
		if(txpEffects[i] >= 1)
			L1_effects +=  (txpEffects[i] -1) / ( max_effect_Thermo -1 ) / nFactors();
		else
			L1_effects +=  min_effect_Thermo / txpEffects[i] / nFactors();

	}

	return L1_weights + L1_effects;
}


double ExprPar::parameter_L1_norm_interactions() const
{
	double L1_norm_coop = 0;
	int Nparameter = 0;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
			#if NEGATIVE_COOP
			if(factorIntMat( i, j ) >= 1)
				L1_norm_coop +=  (factorIntMat(i, j) - 1) / max_interaction;
			else
				L1_norm_coop +=  min_interaction / factorIntMat( i, j );
			#else
			L1_norm_coop += factorIntMat(i, j)/ max_interaction;
			#endif // NEGATIVE_COOP
        }	
	}
	L1_norm_coop /= nFactors();

	double L1_norm_syn = 0;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
			if(factorSynMat( i, j ) >= 1)
				L1_norm_syn += (factorSynMat(i, j) - 1)/ max_synergy;
			else
				L1_norm_syn +=  min_synergy / factorSynMat( i, j );
		}
	}
	L1_norm_syn /= nFactors();

	return L1_norm_coop + L1_norm_syn;
}


void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const IntMatrix& SynMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( SynMat.isSquare() && SynMat.nRows() == nFactors() );  
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

    // write the synergy matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( SynMat( i, j ) ) {
                    double synergy = searchOption == CONSTRAINED ? infty_transform( log( factorSynMat( i, j ) ), log( min_synergy ), log( max_synergy ) ) : log( factorSynMat( i, j ) ); 
                   pars.push_back( synergy );
                }
            }
        }	       
    }

    #if ORIENTATION
    // write the skew matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( coopMat( i, j ) ) {
                    double skew = searchOption == CONSTRAINED ? infty_transform( factorSkewMat( i, j ), min_skew, max_skew) : factorSkewMat( i, j ); 
                   pars.push_back( skew );
                }
            }
        }	       
    }
    #endif // ORIENTATION
	
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
    #endif // ACCESSIBILITY

}

void ExprPar::print( ostream& os, const vector< string >& motifNames, const vector< string >& seqNames, const IntMatrix& coopMat, const IntMatrix& SynMat) const
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
    os << "Accessibility = " << acc_scale << endl;
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
    // print the interaction interactions
    os << "Cooperativity Factor:"  << endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << motifNames[i] << "\t" << motifNames[j] << "\t" << factorIntMat( i, j ) << "\t" << factorSkewMat(i, j) << endl;			
        }
    }
    os << "Synergy Factor:"  << endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( SynMat( i, j ) ) os << motifNames[i] << "\t" << motifNames[j] << "\t" << factorSynMat( i, j ) << endl;			
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
    fin >> symbol >> eqSign >> value;
    if (symbol != "Accessibility" || eqSign != "=") return RET_ERROR;
    acc_scale = atof( value.c_str() );
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
	string line;
    if (symbol != "Cooperativity" || eqSign != "Factor:") return RET_ERROR;
    // read the cooperative interactions
    string factor1, factor2;
    double coopVal, skewVal;
	getline(fin, line);
	while( getline(fin, line) ){
		if( line == "Synergy Factor:") break;
	    istringstream iss(line);
		iss >> factor1 >> factor2 >> coopVal >> skewVal;
		if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorIntMat( idx1, idx2 ) = coopVal;
        factorIntMat( idx2, idx1 ) = coopVal;
        factorSkewMat( idx1, idx2 ) = skewVal;
        factorSkewMat( idx2, idx1 ) = skewVal;
    }

    // read the synergy interactions
    double SynVal;
	while (getline(fin, line)){
	    istringstream iss(line);
		iss >> factor1 >> factor2 >> SynVal;
        if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorSynMat( idx1, idx2 ) = SynVal;
        factorSynMat( idx2, idx1 ) = SynVal;
    }


    fin.close();
    return fin ? RET_ERROR : 0;
}

int ExprPar::load( const string& file, const vector <string>& seqNames, const vector <string>& motifNames )
{
    // open the file
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl;	exit( 1 ); } 


    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    // read the factor information
	string name_tmp;
    for ( int i = 0; i < nFactors(); i++ ) {
        fin >> name_tmp;
        int idx = factorIdxMap[name_tmp];
		fin >> maxBindingWts[idx] >> txpEffects[idx]; 
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) fin >> repEffects[idx];
    }
    string symbol, eqSign, value;
    #if ACCESSIBILITY
    fin >> symbol >> eqSign >> value;
    if (symbol != "Accessibility" || eqSign != "=") return RET_ERROR;
    acc_scale = atof( value.c_str() );
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
	string line;
    if (symbol != "Cooperativity" || eqSign != "Factor:") return RET_ERROR;
    // read the cooperative interactions
    string factor1, factor2;
    double coopVal, skewVal;
	getline(fin, line);
	while( getline(fin, line) ){
		if( line == "Synergy Factor:") break;
	    istringstream iss(line);
		iss >> factor1 >> factor2 >> coopVal >> skewVal;
		if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorIntMat( idx1, idx2 ) = coopVal;
        factorIntMat( idx2, idx1 ) = coopVal;
        factorSkewMat( idx1, idx2 ) = skewVal;
        factorSkewMat( idx2, idx1 ) = skewVal;
    }

    // read the synergy interactions
    double SynVal;
	while (getline(fin, line)){
	    istringstream iss(line);
		iss >> factor1 >> factor2 >> SynVal;
        if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorSynMat( idx1, idx2 ) = SynVal;
        factorSynMat( idx2, idx1 ) = SynVal;
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

    // adjust the synergy matrix 
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( factorSynMat( i, j ) < ExprPar::min_synergy * ( 1.0 + ExprPar::delta ) ) { 
                factorSynMat( i, j ) *= 2.0; 
                factorSynMat( j, i ) = factorSynMat( i, j ); 
            }
            if ( factorSynMat( i, j ) > ExprPar::max_synergy * ( 1.0 - ExprPar::delta ) ) {
                factorSynMat( i, j ) /= 2.0;
                factorSynMat( j, i ) = factorSynMat( i, j ); 
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
            if ( repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) ) repEffects[i] *= 2.0;
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
    
    // constrain the synergy matrix 
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( factorSynMat( i, j ) < ExprPar::min_synergy * ( 1.0 + ExprPar::delta ) ) { 
                factorSynMat( i, j ) = ExprPar::min_synergy; 
                factorSynMat( j, i ) = factorSynMat( i, j ); 
            }
            if ( factorSynMat( i, j ) > ExprPar::max_synergy * ( 1.0 - ExprPar::delta ) ) {
                factorSynMat( i, j ) = ExprPar::max_synergy;
                factorSynMat( j, i ) = factorSynMat( i, j ); 
            }
        }
    }

    // constrain the skew matrix 
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( factorSkewMat( i, j ) < ExprPar::min_skew * ( 1.0 + ExprPar::delta ) ) { 
                factorSkewMat( i, j ) = ExprPar::min_skew; 
                factorSkewMat( j, i ) = factorSkewMat( i, j ); 
            }
            if ( factorSkewMat( i, j ) > ExprPar::max_skew * ( 1.0 - ExprPar::delta ) ) {
                factorSkewMat( i, j ) = ExprPar::max_skew;
                factorSkewMat( j, i ) = factorSkewMat( i, j ); 
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
    #endif // ACCESSIBILITY

}

ModelType ExprPar::modelOption = DIRECT;
SearchType ExprPar::searchOption = CONSTRAINED;
int ExprPar::estBindingOption = 1;  // 1. estimate binding parameters; 0. not estimate binding parameters
 
// training parameter limits / range

double ExprPar::delta = 0.0001;
double ExprPar::default_acc_scale = 1;
double ExprPar::par_penalty = 1;
double ExprPar::interaction_penalty = 1;
double ExprPar::default_weight = 1.0;
double ExprPar::default_interaction = 1.0;
double ExprPar::default_synergy = 1.0;
double ExprPar::default_skew = 0.0;
double ExprPar::default_effect_Logistic = 0.0;
double ExprPar::default_effect_Thermo = 1.0;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -5.0;
double ExprPar::default_basal_Thermo = 0.01;
double ExprPar::min_acc_scale = 0.001;
double ExprPar::max_acc_scale = 100.0;
double ExprPar::min_basal_Logistic = -9.0;	
double ExprPar::max_basal_Logistic = -1.0;
// double ExprPar::min_effect_Direct = 0.01;
double ExprPar::min_effect_Logistic = -5;	
double ExprPar::max_effect_Logistic = 5;


#if	PARAMETER_SPACE == 0	// small
double ExprPar::min_weight = 0.001;		
double ExprPar::max_weight = 500;	
double ExprPar::min_interaction = 0.0001;	
double ExprPar::max_interaction = 500;
double ExprPar::min_synergy = 0.1;	
double ExprPar::max_synergy = 10;
double ExprPar::min_effect_Thermo = 0.01;	
double ExprPar::max_effect_Thermo = 500;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Thermo = 1.0E-5;	
double ExprPar::max_basal_Thermo = 0.05;

#elif PARAMETER_SPACE == 1	// medium
double ExprPar::min_weight = 0.001;		
double ExprPar::max_weight = 1000;		
double ExprPar::min_interaction = 0.0001;	
double ExprPar::max_interaction = 10000;
double ExprPar::min_synergy = 0.05;	
double ExprPar::max_synergy = 50;
double ExprPar::min_effect_Thermo = 0.001;	
double ExprPar::max_effect_Thermo = 1000;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Thermo = 1.0E-6;	
double ExprPar::max_basal_Thermo = 0.1;

#elif PARAMETER_SPACE == 2	// large
double ExprPar::min_weight = 0.001;		
double ExprPar::max_weight = 1000;		
double ExprPar::min_interaction = 0.001;	
double ExprPar::max_interaction = 1000;
double ExprPar::min_synergy = 0.001;	
double ExprPar::max_synergy = 1000;
double ExprPar::min_skew = -10;	
double ExprPar::max_skew = 10;
double ExprPar::min_effect_Thermo = 0.0001;	
double ExprPar::max_effect_Thermo = 10000;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Thermo = 1.0E-7;	
double ExprPar::max_basal_Thermo = 0.1;
#endif // PARAMETER_SPACE

bool ExprPar::one_qbtm_per_crm = ONE_QBTM;
bool ExprFunc::one_qbtm_per_crm = ONE_QBTM;

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, int _SynDistThr, const ExprPar& _par, const vector< Sequence >& _seqs ) : motifs( _motifs ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ), SynDistThr( _SynDistThr ), par( _par ), seqs( _seqs )
{
    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors ); 
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
}

void ExprFunc::set_sites( SiteVec _sites )
{
	sites = _sites;
    int n = sites.size();
    vector< double > _bindingWts (n,0.0);
    vector< int > _boundaries (n,0);

    // Determin the boundaries
    _boundaries[0] = 0;
    int range = max(coopDistThr, repressionDistThr);
    for (int k = 1; k < n; k++) {
    	int l; 
		for (l = 1; l < k; l++) {
	    	if ((sites[k].start - sites[l].start) <= range or siteOverlap( sites[k], sites[l], motifs )) {break;} 
		}
    	_boundaries[k] = l - 1;
    }	
    set_boundaries(_boundaries);
    
    // compute the Boltzman weights of binding for all sites for func
    _bindingWts[0] = 1.0;
    for ( int k = 1; k < n; k++ ) {
	double access_tmp = 1.0;
	#if ACCESSIBILITY
	access_tmp = exp(  - par.acc_scale * ( 1 - sites[k].accessibility) );
	#endif //ACCESSIBILITY
        _bindingWts[k] = access_tmp * par.maxBindingWts[ sites[k].factorIdx ] * sites[k].wtRatio ;	

    }
   set_bindingWts(_bindingWts); 

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
    int dmax = max( coopDistThr, repressionDistThr) + 1;

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
    int dmax = max( coopDistThr, repressionDistThr);

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
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * factorConcs[sites[i].factorIdx] * par.repEffects[ sites[i].factorIdx ] * sum1; 
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
            	sum += compFactorSyn( sites[ i ], sites[ j ] ) * compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
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
	if(par.factorIntMat(a.factorIdx, b.factorIdx) == 1.0)
		return 1.0;
 
    double maxInt = par.factorIntMat(a.factorIdx, b.factorIdx);
   unsigned dist = abs( a.start - b.start );
    assert( dist >= 0 );
	double spacingTerm = 1;
	if(FactorIntOption == BINARY) 
		#if NEGATIVE_COOP 
	    spacingTerm = ( dist < coopDistThr ? maxInt : 1.0 );
		#else
	    spacingTerm = ( dist < coopDistThr ? maxInt + 1 : 1.0 );
		#endif //NEGATIVE_COOP  
	else if(FactorIntOption == LINEAR){
		double d = float(dist)/float(coopDistThr);
		#if NEGATIVE_COOP 
	    spacingTerm = ( dist < coopDistThr ? maxInt  * (1 - d) + d : 1.0 );
		#else
	    spacingTerm = ( dist < coopDistThr ? maxInt  * (1 - d) + 1 : 1.0 );
		#endif //NEGATIVE_COOP  
	}
	else if(FactorIntOption == GAUSSIAN){
		double d = float(dist)/float(coopDistThr);
		#if NEGATIVE_COOP 
	    spacingTerm = ( dist < coopDistThr ? (maxInt - 1) * exp( - 4.5 * ( d * d ) ) + 1: 1.0 );	// Sigma is one third of coopDistThr
		#else
	    spacingTerm = ( dist < coopDistThr ? maxInt * exp( - 4.5 * ( d * d ) ) + 1: 1.0 );	// Sigma is one third of coopDistThr
		#endif //NEGATIVE_COOP  
	}

    #if ORIENTATION
    double x = 3 * float(dist)/float(coopDistThr); // Sigma is one third of coopDistThr
    double orientationTerm = 1; 
    if( a.strand == b.strand )
        orientationTerm = 1 + erf(par.factorSkewMat(a.factorIdx, b.factorIdx) * x);
    else
        orientationTerm = 1 + erf( - par.factorSkewMat(a.factorIdx, b.factorIdx) * x);
    return spacingTerm * orientationTerm;
    #else
    return spacingTerm;
    #endif //ORIENTATION
}

// The interaction function for TOSCA
double ExprFunc::compFactorInt( int t_1, int t_2, int _dist ) const
{
    double maxInt = par.factorIntMat( t_1, t_2 );
    int dist = _dist;
    assert( dist >= 0 );
    assert( dist >= 0 );
	double spacingTerm = 1;
	if(FactorIntOption == BINARY) 
	    spacingTerm = ( dist < coopDistThr ? maxInt + 1: 1.0 );
	else if(FactorIntOption == LINEAR){
		double d = float(dist)/float(coopDistThr);
	    spacingTerm = ( dist < coopDistThr ? maxInt * (1 - d) + 1 : 1.0 );
	}
	else if(FactorIntOption == GAUSSIAN){
		double d = float(dist)/float(coopDistThr);
	    spacingTerm = ( dist < coopDistThr ? maxInt * exp( - 4.5 * ( d * d ) ) + 1: 1.0 );	// Sigma is one third of coopDistThr
	}
    return spacingTerm;
}



// The synergy function for direct model
double ExprFunc::compFactorSyn( const Site& a, const Site& b ) const
{
	if(par.factorSynMat(a.factorIdx, b.factorIdx) == 1.0)
		return 1.0;

    double maxSyn = par.factorSynMat(a.factorIdx, b.factorIdx);

    unsigned dist = abs( a.start - b.start );
    assert( dist >= 0 );
	double spacingTerm = 1;
	double d = float(dist)/float(SynDistThr);
    spacingTerm = ( dist < SynDistThr ? (maxSyn - 1) * exp( - 4.5 * ( d * d ) ) + 1 : 1.0 );

    return spacingTerm;
}

bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
// 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}
/*
ExprPredictor::ExprPredictor(){

    vector< SiteVec >& seqSites;		
    const vector< int > _seqLengths(1,0); 
	seqLengths = &_seqLengths;       
    Matrix exprData;		
    vector< Motif >& motifs;		
    Matrix& factorExprData;		
    const vector< bool > _actIndicators(1, true);   
	actIndicators = &_actIndicators;
    const vector< bool > _repIndicators(1, true);   
	repIndicators = &_repIndicators;
    IntMatrix _repressionMat;    
}
*/
ExprPredictor::ExprPredictor( const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const IntMatrix& _coopMat, const IntMatrix& _SynMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, int _repressionDistThr, int _coopDistThr, int _SynDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector <string>& _seqNames, const vector< Sequence >& _seqs  ) : seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), coopDistThr( _coopDistThr ), SynDistThr( _SynDistThr ), indicator_bool ( _indicator_bool ), seqs(_seqs)
{

    motifNames = _motifNames;
    coopMat = _coopMat;
    SynMat = _SynMat;
    seqNames = _seqNames;

    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
    assert( coopMat.isSquare() && coopMat.isSymmetric() && coopMat.nRows() == nFactors() );
    assert( SynMat.isSquare() && SynMat.isSymmetric() && SynMat.nRows() == nFactors() );
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

    ExprFunc* func = createExprFunc( par );
    double squaredErr = 0;
    double correlation = 0;
    double pgp_score = 0;
    #if TIMER 
    timeval start, end;
    gettimeofday(&start, 0);
    #endif // TIMER
	if(one_qbtm_per_crm){        
		for ( int i = 0; i < nSeqs(); i++ ) {
			// Initiate the sites for func
			func->set_sites(seqSites[ i ]);        
			vector< double > predictedEfficiency (nConds(), -1);
		    vector< double > observedExprs (nConds(), 1);
		    vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
		    
			//#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < nConds(); j++ ) {		
				concs[j] = factorExprData.getCol( j );
				predictedEfficiency[j] = func->predictExpr_scalefree( seqLengths[ i ], concs[j], i );	
				observedExprs[j] = exprData( i, j );									// observed expression for the i-th sequence at the j-th condition
			}
			correlation += train_btr(predictedEfficiency, observedExprs, i);

			double beta = 1.0;
			//squaredErr += least_square( predictedExprs, observedExprs, beta );
			//correlation += corr( predictedExprs, observedExprs ); 
			//pgp_score += pgp( predictedExprs, observedExprs, beta );
		}	
	}
	else{
		for ( int i = 0; i < nSeqs(); i++ ) {
			// Initiate the sites for func
			func->set_sites(seqSites[ i ]);        
			vector< double > predictedExprs (nConds(), -1);
		    vector< double > observedExprs (nConds(), 1);
		    vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
		    
			//#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < nConds(); j++ ) {		
				concs[j] = factorExprData.getCol( j );
				predictedExprs[j] = func->predictExpr( seqLengths[ i ], concs[j], i );	
				observedExprs[j] = exprData( i, j );									// observed expression for the i-th sequence at the j-th condition
			}

			double beta = 1.0;
			squaredErr += least_square( predictedExprs, observedExprs, beta );
			correlation += corr( predictedExprs, observedExprs ); 
			pgp_score += pgp( predictedExprs, observedExprs, beta );
		}	
	}
    #if TIMER 
    gettimeofday(&end, 0);
    cout << "Time " << (end.tv_sec-start.tv_sec)+1e-6*(end.tv_usec-start.tv_usec) << endl;
    #endif // TIMER

    delete func;
    obj_corr = correlation / nSeqs();
    obj_sse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
    obj_pgp = pgp_score / nSeqs();
	double penalty = 0;
	if (PenaltyOption == L1) 

		penalty = par.par_penalty * par.parameter_L1_norm() + par.interaction_penalty * par.parameter_L1_norm_interactions();
	if (PenaltyOption == L2) 
		penalty = par.par_penalty * par.parameter_L2_norm() + par.interaction_penalty * par.parameter_L2_norm_interactions();

    if (objOption == SSE)	return obj_sse - penalty;
    else if (objOption == CORR)	return -obj_corr + penalty;
    else if (objOption == PGP)	return -obj_pgp + penalty;
    return 0;
}


double ExprPredictor::objFunc( const ExprPar& par, int crm ) 
{
    ExprFunc* func = createExprFunc( par );
    double squaredErr = 0;
    double correlation = 0;
    double pgp_score = 0;

   	// Initiate the sites for func
	func->set_sites(seqSites[ crm ]);        
	vector< double > predictedExprs (nConds(), -1);
    vector< double > observedExprs (nConds(), 1);
    vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
        
	//#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < nConds(); j++ ) {		
		concs[j] = factorExprData.getCol( j );
		predictedExprs[j] = func->predictExpr(seqLengths[ crm ], concs[j], crm);	
		observedExprs[j] = exprData(crm, j);									// observed expression for the i-th sequence at the j-th condition
	}

	double beta = 1.0;
	squaredErr += least_square( predictedExprs, observedExprs, beta );
	correlation += corr( predictedExprs, observedExprs ); 
	pgp_score += pgp( predictedExprs, observedExprs, beta );

    delete func;
    obj_corr = correlation;
    obj_sse = sqrt( squaredErr / ( nConds() ) ); 
    obj_pgp = pgp_score;
	double penalty = 0;
	if (PenaltyOption == L1) 
		penalty = par.par_penalty * par.parameter_L1_norm() + par.interaction_penalty * par.parameter_L1_norm_interactions();
	if (PenaltyOption == L2) 
		penalty = par.par_penalty * par.parameter_L2_norm() + par.interaction_penalty * par.parameter_L2_norm_interactions();

    if (objOption == SSE)	return obj_sse - penalty;
    else if (objOption == CORR)	return -obj_corr + penalty;
    else if (objOption == PGP)	return -obj_pgp + penalty;
    return 0;
}

int ExprPredictor::train( const ExprPar& par_init )
{   
    par_model = par_init;	// Initialise the model parameter
    par_curr = par_init;	// The working parameter, which get saved in case of an emergancy
    signal(SIGINT, catch_signal);	// Initialise signal catching
    obj_model = objFunc( par_model );

    if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ){ 
		par_model.constrain_parameters(); 
	    }//par_model.adjust(); }
    if ( nAlternations == 0 ) return 0;
    ExprPar par_result;
    double obj_result;

    for ( int i = 0; i < nAlternations; i++ ) {
		if(optimizationOption == CMAES){
			double tolerance = 1e-4;
			cout << "CMA-ES minimization (sigma = " << cmaes_sigma << "):" << endl; 
			cmaes_minimize(par_result, obj_result, cmaes_sigma, tolerance);
		}
		else if(optimizationOption == Simplex){
			cout << "Simplex minimization " << i + 1 << " of " << nAlternations << ":" << endl; 
			simplex_minimize(par_result, obj_result);
		}
		else if(optimizationOption == BFGS){
			cout << "Gradient minimization step " << i + 1 << " of " << nAlternations << ":" << endl; 
			gradient_minimize(par_result, obj_result);
		}
		if(one_qbtm_per_crm){
			par_result.basalTxps = par_model.basalTxps;
		}
		if(obj_result <= obj_model){ 
			par_model = par_result;
			obj_model = obj_result;
			save_param();
			}

    }
	
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

int ExprPredictor::predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num )
{
    targetExprs.clear();
    targetExprs.resize( nConds() );
    ExprFunc* func = createExprFunc( par_model );
    func->set_sites(targetSites);
	
	vector< double > predictedEfficiency (nConds(), -1);
    vector< double > observedExprs = exprData.getRow( seq_num );
    vector < vector < double > > concs (nConds(), vector <double> (factorExprData.nRows(), 0) );
	if(one_qbtm_per_crm){        
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < nConds(); j++ ) {		
			concs[j] = factorExprData.getCol( j );
			predictedEfficiency[j] = func->predictExpr_scalefree( targetSeqLength, concs[j], seq_num );	
		}

		train_btr(predictedEfficiency, observedExprs, seq_num);

		#pragma omp parallel for schedule(dynamic)
		for ( int j = 0; j < nConds(); j++ ) {
			double efficiency = par_model.basalTxps[seq_num] * predictedEfficiency[j];
		    double predicted = efficiency / (1 + efficiency);	
		    targetExprs[j] = predicted;
		}
	}
	else{
		#pragma omp parallel for schedule(dynamic)
		for ( int j = 0; j < nConds(); j++ ) {
			concs[j] = factorExprData.getCol( j );
			double predicted = func->predictExpr( targetSeqLength, concs[j], seq_num );       
		    targetExprs[j] = predicted;
		}
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
OptimizationType ExprPredictor::optimizationOption = OPTIMIZATION_ALGORITHM;
PenaltyType ExprPredictor::PenaltyOption = PARAMETER_PENALTY;
FactorIntType ExprFunc::FactorIntOption = FactorIntFunc;

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
double ExprPredictor::min_delta_f_PGP = 1.0E-8;
int ExprPredictor::nSimplexIters = 20;
int ExprPredictor::nCMAESIters = 10000;
int ExprPredictor::nGradientIters = 5000;
double ExprPredictor::cmaes_sigma = 0.1;

bool ExprPredictor::one_qbtm_per_crm = ONE_QBTM;

// Initialise static members as empty
ExprPar ExprPredictor::par_curr;
IntMatrix ExprPredictor::coopMat = IntMatrix();
IntMatrix ExprPredictor::SynMat = IntMatrix();
vector <string> ExprPredictor::motifNames = vector <string>();
vector <string> ExprPredictor::seqNames = vector <string>();

// Initialise a global pointer that can smuggle ExprPredictor into FitFunc (I know! I hate it, too.)
void* global_pointer;

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

    // sample the interaction matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                double rand_synergy = exp( gsl_ran_flat( rng, log( ExprPar::min_synergy ), log( ExprPar::max_synergy ) ) );
                if ( SynMat( i, j ) ) par.factorSynMat( i, j ) = rand_synergy;
            }
        }
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = i + 1; j < nFactors(); j++ ) {
                par.factorSynMat( i, j ) = par.factorSynMat( j, i );
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

    // test the synergy matrix
    if ( modelOption != LOGISTIC ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( par.factorSynMat( i, j ) < ExprPar::min_synergy * ( 1.0 + ExprPar::delta ) || par.factorSynMat( i, j ) > ExprPar::max_synergy * ( 1.0 - ExprPar::delta ) )
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

   // print the synergy matrix
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( SynMat( i, j ) ) cout << par.factorSynMat( i, j ) << "\t";
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
	cout << par.acc_scale << "\t";

    cout << endl;
}

ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{

    return new ExprFunc( motifs, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, SynDistThr, par, seqs );
}


int indices_of_crm_in_gene[] = {
	5, 11, 17, 23, 29
};

double ExprPredictor::comp_impact( const ExprPar& par, int tf ) 
{
	// Calculate the objecttive function with the factor tf basicly deleted
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
	par_deleted.maxBindingWts[tf] = ExprPar::min_weight;
	par_deleted.txpEffects[tf] = 1.0;
	double obj_deleted = objFunc(par_deleted);
	double obj_full = objFunc(par_full);	
	double impact = (obj_full - obj_deleted) / obj_full;

	if (objOption == SSE) return -impact;	// SSE gets minimised
	else return impact;
}

double ExprPredictor::comp_impact_coop( const ExprPar& par, int tf ) 
{
	// Calculate the objecttive function without cooperativity between tf1 and tf2
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;    
    # if NEGATIVE_COOP
	vector<double > zero_vector (nFactors(),1);
    #else
	vector<double > zero_vector (nFactors(),ExprPar::min_interaction);
    #endif //NEGATIVE_COOP
	par_deleted.factorIntMat.setRow( tf, zero_vector );
	par_deleted.factorIntMat.setCol( tf, zero_vector );
	double obj_deleted = objFunc(par_deleted);
	double obj_full = objFunc(par_full);	
	double impact = (obj_full - obj_deleted);	

	if (objOption == SSE)	return	impact;	// SSE gets minimised
	else return -impact;
}

double ExprPredictor::comp_impact_acc( const ExprPar& par ) 
{
	// Calculate the objecttive function without accessibility
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
	par_deleted.acc_scale = ExprPar::min_acc_scale;
	double obj_deleted = objFunc(par_deleted);
	double obj_full = objFunc(par_full);	
	double impact = (obj_full - obj_deleted) / obj_full;	

	if (objOption == SSE)	return	-impact;	// SSE gets minimised
	else return impact;
}

double ExprPredictor::comp_impact( const ExprPar& par, int tf, int crm ) 
{
	// Calculate the objecttive function with the factor tf basicly deleted
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
	par_deleted.maxBindingWts[tf] = ExprPar::min_weight;
	par_deleted.txpEffects[tf] = 1.0;
	double obj_deleted = objFunc(par_deleted, crm);
	double obj_full = objFunc(par_full, crm);	
	double impact = obj_full - obj_deleted;			

	if (objOption == SSE)	return	impact;	// SSE gets minimised
	else return -impact;
}

double ExprPredictor::comp_impact_coop( const ExprPar& par, int tf, int crm ) 
{
	// Calculate the objecttive function without cooperativity between tf1 and tf2
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
    # if NEGATIVE_COOP
	vector<double > zero_vector (nFactors(),1);
    #else
	vector<double > zero_vector (nFactors(),ExprPar::min_interaction);
    #endif //NEGATIVE_COOP
	par_deleted.factorIntMat.setRow( tf, zero_vector );
	par_deleted.factorIntMat.setCol( tf, zero_vector );
	double obj_deleted = objFunc(par_deleted, crm);
	double obj_full = objFunc(par_full, crm);	
	double impact = obj_full - obj_deleted;

	if (objOption == SSE)	return	impact;	// SSE gets minimised
	else return -impact;	
}


double ExprPredictor::comp_impact_coop_pair( const ExprPar& par, int tf1, int tf2 ) 
{
	// Calculate the objecttive function without cooperativity between tf1 and tf2
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
	#if NEGATIVE_COOP 
	par_deleted.factorIntMat( tf1, tf2 ) = 1;
	#else
	par_deleted.factorIntMat( tf1, tf2 ) = ExprPar::min_interaction;
	#endif //NEGATIVE_COOP  
	double obj_deleted = objFunc(par_deleted);
	double obj_full = objFunc(par_full);	
	double impact = obj_full - obj_deleted;

	if (objOption == SSE)	return	impact;	// SSE gets minimised
	else return -impact;	
}


double ExprPredictor::comp_impact_synergy_pair( const ExprPar& par, int tf1, int tf2 ) 
{
	// Calculate the objecttive function without cooperativity between tf1 and tf2
	ExprPar par_deleted = par;
	ExprPar par_full = par;
	par_full.par_penalty = 0;
	par_deleted.par_penalty = 0;
	#if NEGATIVE_COOP 
	par_deleted.factorSynMat( tf1, tf2 ) = 1;
	#else
	par_deleted.factorSynMat( tf1, tf2 ) = ExprPar::min_interaction;
	#endif //NEGATIVE_COOP  
	double obj_deleted = objFunc(par_deleted);
	double obj_full = objFunc(par_full);	
	double impact = obj_full - obj_deleted;

	if (objOption == SSE)	return	impact;	// SSE gets minimised
	else return -impact;	
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
	int range = max(coopDistThr, repressionDistThr);
        for ( int k = 1; k < n; k++ ) {
	       	int l; 
			for ( l = k - 1; l >= 0; l-- ) {
			    if ( ( seqSites[i][k].start - seqSites[i][l].start ) > range or siteOverlap(seqSites[i][k], seqSites[i][l], motifs ) ) break; 
			}
	   		_boundaries[k] = l ;
	}	
        func->set_boundaries(_boundaries);
    
	// compute the Boltzman weights of binding for all sites for func
        _bindingWts[0] = 1.0;
        for ( int k = 1; k < n; k++ ) {
		double access_tmp = 1.0;
		#if ACCESSIBILITY
		access_tmp = exp( -par_model.acc_scale * (1- seqSites[i][k].accessibility) ) ;
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
    //#pragma omp parallel for schedule(dynamic)
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
    // extract initial parameters
    vector < double > pars;
    par_model.getFreePars( pars, coopMat, SynMat, actIndicators, repIndicators ); 
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

	
        par_curr = ExprPar ( pars, coopMat, SynMat, actIndicators, repIndicators, nSeqs() );

	//Hassan end
	// check if the current values of parameters are valid
        //the following line should be uncommented if you remove all the changes by Hassan
	//ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, SynMat actIndicators, repIndicators );
        
	
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
    		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	#ifdef SHORT_OUTPUT
	if(iter % SHORT_OUTPUT == 0){
   		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	cout << "======================================" << endl;
	cout << "======================================" << endl;
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
	cout << "======================================" << endl;
	par_curr.print( cout, motifNames, seqNames, coopMat, SynMat);
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

        par_result = ExprPar ( pars, coopMat, SynMat, actIndicators, repIndicators, nSeqs() );
	//Hassan end
	//uncomment the following line if you remove all the changes by Hassan
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, SynMat, actIndicators, repIndicators );
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
    par_model.getFreePars( pars, coopMat, SynMat, actIndicators, repIndicators ); 
            
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
 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
//    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr; // Chose Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm 
		
    // create the minimizer
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = 2, tol = 0.1;
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

        par_curr = ExprPar ( pars, coopMat, SynMat, actIndicators, repIndicators, nSeqs() );
	
	//Hassan end
        //ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, SynMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;
    

        // check for stopping condition
//        double f_curr = s->f;
//        double delta_f = abs( f_curr - f_prev ); 
//		cout << delta_f << endl;   
//     	if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
//        if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
//        if ( objOption == PGP && delta_f < min_delta_f_PGP ) break;
        
        status = gsl_multimin_test_gradient( s->gradient, 5e-4 );
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

         // print the current parameter and function values
	#if FILE_OUTPUT
	if(iter % 1000 == 0){
    		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	#ifdef SHORT_OUTPUT
	if(iter % SHORT_OUTPUT == 0){
   		printf( "\r %zu \t SSE = %8.5f \t Corr = %8.5f \t PGP = %8.5f", iter, obj_sse, obj_corr, obj_pgp);
		fflush(stdout);
	}	
	#else
	cout << "======================================" << endl;
	cout << "======================================" << endl;
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
	cout << "======================================" << endl;
	par_curr.print( cout, motifNames, seqNames, coopMat, SynMat);
	cout << "======================================" << endl;
	cout << "======================================" << endl << endl;
	#endif // SHORT_OUTPUT
	#endif // FILE_OUTPUT
	// Save parameters
	if( iter % 1000 == 0 ) save_param();

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

    par_result = ExprPar ( pars, coopMat, SynMat, actIndicators, repIndicators, nSeqs() );
	//Hassan end
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, SynMat, actIndicators, repIndicators );
    obj_result = s->f;
    
    // free the minimizer
    gsl_vector_free( x );    
    gsl_multimin_fdfminimizer_free( s );
    
    return 0;
}

libcmaes::ProgressFunc<libcmaes::CMAParameters<>,libcmaes::CMASolutions> select_time = [](const libcmaes::CMAParameters<> &cmaparams, const libcmaes::CMASolutions &cmasols)
{
	#if FILE_OUTPUT
	if (cmasols.niter() % 100 == 0){
		double current_score = - cmasols.best_candidate().get_fvalue();
		double best_score = - cmasols.get_best_seen_candidate().get_fvalue();		 	
		printf( "\r %i \t current score = %8.5f \t best score = %8.5f", cmasols.niter(), current_score, best_score);
		fflush(stdout);
	}
	#else
	if (cmasols.niter() % 1 == 0){
		double current_score = - cmasols.best_candidate().get_fvalue();
		double best_score = - cmasols.get_best_seen_candidate().get_fvalue();		 	
		printf( "\r %i \t current score = %8.5f \t best score = %8.5f", cmasols.niter(), current_score, best_score);
		fflush(stdout);
	}
	#endif // FILE_OUTPUT
  return 0;
};


libcmaes::FitFunc obj_func_wrapper = [](const double *x, const int N)
{

    ExprPredictor* global_predictor = (ExprPredictor*)global_pointer;

	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;

	for( int index = 0; index < global_predictor->indicator_bool.size(); index ++ ){
		if( global_predictor->indicator_bool[index]){
			all_pars.push_back( x[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( global_predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}

    ExprPar par_tmp(all_pars, global_predictor->getCoopMat(), global_predictor->getSynMat(), global_predictor->getActIndicators(), global_predictor->getRepIndicators(), global_predictor -> nSeqs());
	global_predictor->par_curr = par_tmp;
    double obj = global_predictor->objFunc(par_tmp);
    return obj;
};


int ExprPredictor::cmaes_minimize(ExprPar& par_result, double& obj_result, double sigma, double tolerance) 
{
    // extract initial parameters
    vector < double > pars;
    par_model.getFreePars(pars, coopMat, SynMat, actIndicators, repIndicators); 
        
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

	libcmaes::CMAParameters<> cmaparams(pars, sigma);
	cout << "Number offspring: " << cmaparams.lambda() << endl;
	cmaparams.set_ftolerance(tolerance);	
	cmaparams.set_algo(aCMAES);
	cmaparams.set_max_iter(nCMAESIters);
	cmaparams.set_mt_feval(true);
	global_pointer = (void*)this;

    #if MONITOR_PARAMS
	string fname;
	for(int idx = 1; idx < 1000; idx ++){
		fname = "monitor_params_"+to_string(idx)+".dat";
		if(!fexists(fname))
			break;
	}

	cmaparams.set_fplot(fname);
    #endif //MONITOR_PARAMS
	libcmaes::CMASolutions cmasols = libcmaes::cmaes<>(obj_func_wrapper, cmaparams, select_time);    


	libcmaes::Candidate best_candidate = cmasols.get_best_seen_candidate();
	free_pars = best_candidate.get_x();
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

    par_result = ExprPar(pars, coopMat, SynMat, actIndicators, repIndicators, nSeqs());
    obj_result = best_candidate.get_fvalue();	

    return 0;
}

double ExprPredictor::train_btr(vector< double >& predictedEfficiency, vector< double >& observedExprs, int i) 
{
	double btr = par_model.basalTxps[i];
	double btr_log = log(btr);
	double step_size = 1;
	double delta_corr = -corr_gradient( predictedEfficiency, observedExprs, btr);
	double delta_log_corr = delta_corr * btr;

	int Nsteps = 0;
	while (abs(delta_log_corr * step_size) > 1e-5 && Nsteps < 5000){
		Nsteps += 1;
		btr_log -= step_size * delta_log_corr;
		btr = exp(btr_log);
		delta_corr = -corr_gradient( predictedEfficiency, observedExprs, btr);
		delta_log_corr = delta_corr * btr;
		// Optionally use a cooling process
		//step_size *= 0.99;
		if(btr < ExprPar::min_basal_Thermo) {
			btr = ExprPar::min_basal_Thermo;
			break;
		} 
		if(btr > ExprPar::max_basal_Thermo) {
			btr = ExprPar::max_basal_Thermo;
			break;
		} 
	}
	par_model.basalTxps[i] = btr;
	par_curr.basalTxps[i] = btr;
	double corr = corr_scalefree( predictedEfficiency, observedExprs, btr);
	return corr;
}

// function to save parameters to file
int ExprPredictor::save_param()
{
	ofstream fparam_sm( "param.save" );
	par_curr.print( fparam_sm, motifNames, seqNames, getCoopMat(), getSynMat());
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
		if( predictor->indicator_bool[index]){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getSynMat(), predictor->getActIndicators(), predictor->getRepIndicators(), predictor -> nSeqs() );
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getSynMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
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
