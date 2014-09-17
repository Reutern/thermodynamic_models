/*****************************************************
* Train and test the expression model
* Input: sequence, expression, motif, factor expr, cooperativity rules, 
*   activator list, repression rules
* Output: trained model and expression of training data
* File formats: 
* (1) Sequence: Fasta format
* (2) Expression: one line per sequence
*       <seq_name expr_1 ... expr_m>
* (3) Motif: Stubb format
* (4) Factor expression: one line per factor
*       <factor expr_1 ... expr_m>
* (5) Cooperativity: list of cooperative factor pairs 
* (6) Factor roles: the role field is either 1 (yes) or 0 (no)
*       <factor activator_role repressor_role>
* (7) Repression: list of <R A> where R represses A
* (8) Parameters: the format is:
*       <factor binding activation repression>
*       <factor1 factor2 coop>
*       <basal_transcription = x>
*     where repression is optional, and the coop. lines are optional.
* Note that (5), (6), (7) and (8) may be empty
******************************************************/
#include "ExprPredictor.h"
#include "OccPredictor.h"
#include <time.h>
#include <math.h> // for fmod
#include "param.h" 


int main( int argc, char* argv[] ) 
{

    // Set the timer
    long time_start = time(NULL);
  
    // command line processing
    string seqFile, test_seqFile, annFile, exprFile, test_exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, print_parFile, axis_wtFile;
    string outFile, occFile;     // output files
    int coopDistThr = 50;
    double factorIntSigma = 25.0;   // sigma parameter for the Gaussian interaction function
    int repressionDistThr = 50;
    int maxContact = 1;
	vector<double> eTF (8);

	eTF[0] = 0.7 ;
	eTF[1] = 0.7 ;
	eTF[2] = 0.7 ;
	eTF[3] = 0.7 ;
	eTF[4] = 0.7 ;
	eTF[5] = 0.7 ;
	eTF[6] = 0.7 ;
	eTF[7] = 0.7 ;

	string free_fix_indicator_filename;
	ExprPredictor::one_qbtm_per_crm = ONE_QBTM;
	ExprPar::one_qbtm_per_crm = ONE_QBTM;
	ExprFunc::one_qbtm_per_crm = ONE_QBTM;

    ExprPredictor::nAlternations = 0;
    for ( int i = 1; i < argc; i++ ) {
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
        else if ( !strcmp( "-ts", argv[ i ] ) )
            test_seqFile = argv[ ++i ];
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];      
        else if ( !strcmp( "-e", argv[ i ] ) )
            exprFile = argv[ ++i ];            
        else if ( !strcmp( "-te", argv[ i ] ) )
            test_exprFile = argv[ ++i ];            
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
        else if ( !strcmp( "-f", argv[ i ] ) )
            factorExprFile = argv[ ++i ];    
        else if ( !strcmp( "-o", argv[ i ] ) )
            ExprPredictor::modelOption = getModelOption( argv[++i] );
        else if ( !strcmp( "-c", argv[ i ] ) )
            coopFile = argv[ ++i ];
        else if ( !strcmp( "-i", argv[ i ] ) )
            factorInfoFile = argv[ ++i ];            
        else if ( !strcmp( "-r", argv[ i ] ) )
            repressionFile = argv[ ++i ];  
        else if ( !strcmp( "-oo", argv[ i ] ) )
            ExprPredictor::objOption = getObjOption( argv[++i] );    
        else if ( !strcmp( "-mc", argv[i] ) )
            maxContact = atoi( argv[++i] );
        else if ( !strcmp( "-fo", argv[i] ) )
            outFile = argv[++i];    
        else if ( !strcmp( "-oc", argv[i] ) )
            occFile = argv[++i];    
        else if ( !strcmp( "-p", argv[i] ) )
            parFile = argv[++i]; 
        else if ( !strcmp( "-pp", argv[i] ) )
            print_parFile = argv[++i]; 
	else if ( !strcmp( "-wt", argv[ i ]) )
		axis_wtFile = argv[ ++ i ];
        else if ( !strcmp( "-ct", argv[i] ) )
            coopDistThr = atof( argv[++i] ); 
        else if ( !strcmp( "-sigma", argv[i] ) )
            factorIntSigma = atof( argv[++i] );    
        else if ( !strcmp( "-rt", argv[i] ) )
            repressionDistThr = atof( argv[++i] );    
        else if ( !strcmp( "-na", argv[i] ) )
            ExprPredictor::nAlternations = atoi( argv[++i] );    
        else if ( !strcmp( "-ff", argv[i] ) )
            free_fix_indicator_filename = argv[++i];    
        else if ( !strcmp( "-oq", argv[i] ) ){
            	ExprPredictor::one_qbtm_per_crm = ONE_QBTM;    
		ExprPar::one_qbtm_per_crm = ONE_QBTM;
		ExprFunc::one_qbtm_per_crm = ONE_QBTM;
	}
        else if ( !strcmp( "-et", argv[i] ) ) {}
           // eTF = atof( argv[ ++i ] );    
    }
    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -ts test_seqFile -e exprFile -te test_exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -pp print_parFile -rt repressionDistThr -na nAlternations -ct coopDistThr -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }           


//     bool readSites = false;     // whether read sites (if true) or read sequences 
    
    // additional control parameters
    double gcContent = 0.5;
    FactorIntType intOption = BINARY;     // type of interaction function
    ExprPar::searchOption = CONSTRAINED;      // search option: unconstrained; constrained. 
    ExprPar::estBindingOption = 1;

    ExprPredictor::nRandStarts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 20000;
    ExprPredictor::nGradientIters = 2000;

    int rval;
    vector< vector< double > > data;    // buffer for reading matrix data
    vector< string > labels;    // buffer for reading the labels of matrix data
    string factor1, factor2;    // buffer for reading factor pairs

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();

    // read the expression data
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
  //  assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ){
    	if( labels[ i ] != seqNames[ i ] ) cout << labels[i] << seqNames[i] << endl;
	//assert( labels[i] == seqNames[i] );
    }
    Matrix exprData( data ); 
    int nConds = exprData.nCols();
    
    // read the motifs
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    rval = readMotifs( motifFile, background, motifs, motifNames ); 
    assert( rval != RET_ERROR );
    int nFactors = motifs.size();

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
     
    // read the factor expression data
    labels.clear();
    data.clear();
    rval = readMatrix( factorExprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nFactors && condNames.size() == nConds );
    for ( int i = 0; i < nFactors; i++ ) assert( labels[i] == motifNames[i] );
    Matrix factorExprData( data ); 
    assert( factorExprData.nCols() == nConds ); 
	//initialize the energy threshold factors
	vector < double > energyThrFactors;
	energyThrFactors.clear( );
	for ( int index = 0; index < nFactors; index++ ){
		//cout << motifNames[index] << " " << eTF[index] << endl;  // Print eTF for all TF
		energyThrFactors.push_back( eTF[index] );
	}
    // site representation of the sequences
	
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrFactors );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqs; i++ ) {
		//cout << "Annotated sites for CRM: " << seqNames[i] << endl;
            	ann.annot( seqs[ i ], seqSites[ i ] );
    		seqSites[i].insert( seqSites[i].begin(), Site() );  // start with a pseudo-site at position 0 
            	seqLengths[i] = seqs[i].size();

        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }


    // read the cooperativity matrix 
    int num_of_coop_pairs = 0;
    IntMatrix coopMat( nFactors, nFactors, false );
    if ( !coopFile.empty() ) {
        ifstream fcoop( coopFile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the cooperativity file " << coopFile << endl;
            exit( 1 );
        }  
        while ( fcoop >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            if( coopMat( idx1, idx2 ) == false && coopMat( idx2, idx1 ) == false ){
	    	coopMat( idx1, idx2 ) = true;
	    	coopMat( idx2, idx1 ) = true;
	    	num_of_coop_pairs ++;
	    }
        }        
    } 

    // read the roles of factors
    vector< bool > actIndicators( nFactors, false );
    vector< bool > repIndicators( nFactors, false );
    if ( !factorInfoFile.empty() ) {
        ifstream finfo( factorInfoFile.c_str() );
        if ( !finfo ) {
            cerr << "Cannot open the factor information file " << factorInfoFile << endl;
            exit( 1 );
        }      
        string name;
        int i = 0, actRole, repRole;
        while ( finfo >> name >> actRole >> repRole ) {
            assert( name == motifNames[i] );
            if ( actRole ) actIndicators[i] = true;
            if ( repRole ) repIndicators[i] = true;
            i++;
        }
    }
    
    // read the repression matrix 
    IntMatrix repressionMat( nFactors, nFactors, false );
    if ( !repressionFile.empty() ) {
        ifstream frepr( repressionFile.c_str() );
        if ( !frepr ) {
            cerr << "Cannot open the repression file " << repressionFile << endl;
            exit( 1 );
        }        
        while ( frepr >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            repressionMat( idx1, idx2 ) = true;
        }        
    }

	// read the axis wt file
	vector < int > axis_start;
	vector < int > axis_end;
	vector < double > axis_wts;

	axis_start.clear();
	axis_end.clear();
	axis_wts.clear();

	if( !axis_wtFile.empty() ){
		ifstream axis_wtInfo ( axis_wtFile.c_str() );
		if( !axis_wtInfo ){
			cerr << "Cannot open the axis weight information file " << axis_wtFile << endl;
			exit( 1 );
		}
		int temp1, temp2; 
		double temp3;
		double temp3_sum = 0;
		while( axis_wtInfo >> temp1 >> temp2 >> temp3 ){
			axis_start.push_back( temp1 );
			axis_end.push_back( temp2 );
			axis_wts.push_back( temp3 );
			temp3_sum += temp3;
		}
		assert( !( temp3_sum > 100 ) && !( temp3_sum < 100 ));
	}
	else{
		axis_start.push_back( 0 );
		axis_end.push_back( condNames.size() - 1 );
		axis_wts.push_back( 100 );
	}
	
	vector <bool> indicator_bool;
	indicator_bool.clear();
	if( !free_fix_indicator_filename.empty() ){
		ifstream free_fix_indicator_file ( free_fix_indicator_filename.c_str() );
		while( !free_fix_indicator_file.eof( ) ){
			int indicator_var;
			free_fix_indicator_file >> indicator_var;
			assert ( indicator_var == 0 || indicator_var == 1 );
			indicator_bool.push_back( indicator_var );
		}
	}
	else{
		//for binding weights, coop pairs and transcriptional effects
		for( int index = 0; index < nFactors + num_of_coop_pairs + nFactors; index++ ){
			indicator_bool.push_back( true );
		}
		if( ExprPredictor::one_qbtm_per_crm ){
			for( int index = 0; index < nSeqs; index++ ){
				indicator_bool.push_back( true );
			}
		}
		else{
			indicator_bool.push_back( true );
		}
	}

    // print the parameters for running the analysis
    cout << "Program settings: " << endl; 
    cout << "Model = " << getModelOptionStr( ExprPredictor::modelOption ) << endl;
    if ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_LIMITED ) {
        cout << "Maximum_Contact = " << maxContact << endl;
    }
    if ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_LIMITED || ExprPredictor::modelOption == CHRMOD_UNLIMITED ) {
        cout << "Repression_Distance_Threshold = " << repressionDistThr << endl;
    }
    cout << "Objective_Function = " << getObjOptionStr( ExprPredictor::objOption ) << endl;
    if ( !coopFile.empty() ) {
        cout << "Interaction_Model = " << getIntOptionStr( intOption ) << endl;
        cout << "Interaction_Distance_Threshold = " << coopDistThr << endl;
        if ( intOption == GAUSSIAN ) cout << "Sigma = " << factorIntSigma << endl;
    }
    cout << "Search_Option = " << getSearchOptionStr( ExprPar::searchOption ) << endl;
    cout << "Num_Random_Starts = " << ExprPredictor::nRandStarts << endl;
    cout << "Energy Threshold Factor = " << eTF << endl;
    cout << "Pseudo count = " << PSEUDO_COUNT << endl;
    cout << endl;
    #if PRINT_STATISTICS
    cout << "Statistics: " << endl; 
    cout << "Factors "<< nFactors << "\t " << "Sequences " << nSeqs <<  endl;
    cout << motifNames[0] << " \t " << motifNames[1] << " \t " << motifNames[2] << " \t " << motifNames[3] << " \t " << motifNames[4] << " \t " << motifNames[5] << " \t " << motifNames[6] << " \t " << motifNames[7] << " \t " << "Sum \t Length \t Name" << endl;
    double average_number = 0;
    for(int seqs_idx = 0; seqs_idx < nSeqs; seqs_idx++){
	average_number += seqSites[seqs_idx].size()/nSeqs;
	int sites_count[] = {0,0,0,0,0,0,0,0};
        double weight_count[] = {0,0,0,0,0,0,0,0};
	for( int idx = 1; idx < seqSites[seqs_idx].size() ; idx++ ){
			sites_count[seqSites[seqs_idx][idx].factorIdx]++;
			weight_count[seqSites[seqs_idx][idx].factorIdx] = weight_count[seqSites[seqs_idx][idx].factorIdx] + seqSites[seqs_idx][idx].wtRatio;
		}
        for( int l = 0; l < nFactors; l++){
	        cout <<  round(100*  weight_count[l]) / 100.0 << " \t "; }
	cout << seqSites[seqs_idx].size() - 1 << " \t " << seqNames[seqs_idx] << " \t " << seqLengths[seqs_idx] <<  endl;}
    cout << endl; 
    cout << average_number << endl;
	
    #endif // PRINT_STATISTICS


    #if SAVE_ENERGIES
    cout << "Save site energies" << endl;
    ofstream site_energies;
    site_energies.open ("../data/site_energies.txt");
    for(int tf = 0; tf < nFactors; tf++){
	site_energies << motifNames[tf] << "\t" << factorExprData.getRow(tf) << endl;
	vector <double> weight(nConds);
        for(int n=0; n < nConds; n++){
	weight[n] = 0;
		for(int seqs_idx = 0; seqs_idx < nSeqs; seqs_idx++){
			for( int idx = 1; idx <= seqSites[seqs_idx].size() ; idx++ ){
				if(seqSites[seqs_idx][idx].factorIdx != tf)  continue;
				weight[n] += seqSites[seqs_idx][idx].wtRatio * exprData.getRow(seqs_idx)[n];
			}
		}
	}
	site_energies << motifNames[tf] << "\t" << weight << endl;
    }
    // Write site energies to site_energies.txt	
    site_energies.close();
    cout << "Site energies saved" << endl;
    #endif //SAVE_ENERGIES


 // read the initial parameter values
    ExprPar par_init( nFactors, nSeqs );
    if ( !parFile.empty() ) {
        rval = par_init.load( parFile );
        if ( rval == RET_ERROR ) {
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        } 
    }

    // Initialise the predictor class
    ExprPredictor* predictor = new ExprPredictor( seqSites, seqLengths, exprData, motifs, factorExprData, coopMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, indicator_bool, motifNames, axis_start, axis_end, axis_wts, seqs );
   

    // random number generator
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)

    // model fitting
    predictor->train( par_init, rng );
    gsl_rng_free( rng );

    // print the training results
    ExprPar par = predictor->getPar();

    // print the parameters
    ofstream pout( print_parFile.c_str() );
    if ( !pout ) {
        cout << "Estimated values of parameters:" << endl;
        par.print( cout, motifNames, coopMat );
    }
    else {
        par.print( pout, motifNames, coopMat );
    }
    cout << "Performance = " << setprecision( 5 ) << ( ExprPredictor::objOption == SSE ? predictor->getObj() : -predictor->getObj() ) << endl;

    // print the predictions
    ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }
    fout << "Rows\t" << condNames << endl;

    // print the impact
    cout << "Impact of the TF:" << endl; 
    for(int tf = 0; tf < nFactors; tf++ ){

	double totalweight = 0;
	for(int seqs_idx = 0; seqs_idx < nSeqs; seqs_idx++){
		for( int idx = 1; idx <= seqSites[seqs_idx].size(); idx++ ){
			if(seqSites[seqs_idx][idx].factorIdx != tf)  continue;
			totalweight += seqSites[seqs_idx][idx].wtRatio;
		}
	}
	// effect * binding energy * total weight
	double impact = max(par.txpEffects[tf] , 1/par.txpEffects[tf]) * par.maxBindingWts[tf] * totalweight;
	cout << motifNames[tf] << "\t" << impact << endl;

    }

    #if CALCULATE_OCCUPANCY

    ofstream Occout( occFile.c_str() );

    for ( int seq_idx = 0; seq_idx < nSeqs; seq_idx++ ) {

	printf( "\rOccupancy prediction for sequence: %i ", seq_idx+1);
	fflush(stdout);

    	// Initialise the occupancy predictor
    	OccPredictor* occpred = new OccPredictor( seqSites[seq_idx], motifs, factorExprData, coopMat, coopDistThr, par_init );
	
	vector < vector <double> > tf_occupancy (nConds, vector <double> (nFactors, 0.0));

	    for(int position_idx = 0; position_idx < 100; position_idx++){ 
		vector <int> tf_number (nFactors, 0);
	    	for(int site_idx = 1; site_idx < seqSites[seq_idx].size(); site_idx++){
			tf_number[ seqSites[seq_idx][site_idx].factorIdx ] += 1;
			tf_occupancy[position_idx][ seqSites[seq_idx][site_idx].factorIdx ] += occpred -> predictOcc(site_idx,position_idx);
		}
	
		for(int tf_idx = 0; tf_idx < nFactors; tf_idx++ )
			tf_occupancy[position_idx][tf_idx] = tf_occupancy[position_idx][tf_idx] / tf_number[tf_idx];

	    }

	// First line is the name of the enhancer and the position indices
        Occout << seqNames[seq_idx] ; 
	for ( int position_idx = 0; position_idx < nConds; position_idx++ )
		Occout << "\t" << position_idx;
	Occout << endl;

	// The following lines are the TF names and their occupancies
	for ( int tf_idx = 0; tf_idx < nFactors; tf_idx++ ){
		Occout << motifNames[tf_idx] ; 
       		for ( int position_idx = 0; position_idx < nConds; position_idx++ )
			Occout << "\t" << tf_occupancy[position_idx][tf_idx] ;       // occupancy predictions
	        Occout << endl;
		}

    	delete occpred;		// clean up the predictor

    }
    cout << endl;	
    #endif // CALCULATE_OCCUPANCY

    #if CROSS_VALIDATION

    // read the sequences
    vector< Sequence > test_seqs;
    vector< string > test_seqNames;
    rval = readSequences( test_seqFile, test_seqs, test_seqNames );
    assert( rval != RET_ERROR );
    int test_nSeqs = test_seqs.size();

    // read the expression data
    rval = readMatrix( test_exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
  //  assert( labels.size() == nSeqs );
    for ( int i = 0; i < test_nSeqs; i++ ){
    	if( labels[ i ] != test_seqNames[ i ] ) cout << labels[i] << test_seqNames[i] << endl;
	//assert( labels[i] == seqNames[i] );
    }
    Matrix test_exprData( data ); 

    // site representation of the sequences
	
    vector< SiteVec > test_seqSites( test_nSeqs );
    vector< int > test_seqLengths( test_nSeqs );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < test_nSeqs; i++ ) {
		//cout << "Annotated sites for CRM: " << seqNames[i] << endl;
            	ann.annot( test_seqs[ i ], test_seqSites[ i ] );
            	test_seqLengths[i] = test_seqs[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, test_seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( test_seqs[i], test_seqSites[i] );
            test_seqLengths[i] = test_seqs[i].size();
        }
    }

    for ( int i = 0; i < test_nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor->predict( test_seqSites[i], test_seqLengths[i], targetExprs, i );
        vector< double > observedExprs = test_exprData.getRow( i );
        
        // print the results
        fout << test_seqNames[i] << "\t" << observedExprs << endl;      // observations
        fout << test_seqNames[i];

        for ( int j = 0; j < nConds; j++ ) fout << "\t" << targetExprs[j];       // predictions
        fout << endl;
    }

    #else
    for ( int i = 0; i < nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor->predict( seqSites[i], seqLengths[i], targetExprs, i );
        vector< double > observedExprs = exprData.getRow( i );
        
        // error
        double beta = 1.0; 
        double error =  least_square( targetExprs, observedExprs, beta ) / nConds;

        // print the results
        fout << seqNames[i] << "\t" << observedExprs << endl;      // observations
        fout << seqNames[i]; 

        for ( int j = 0; j < nConds; j++ )
		fout << "\t" << beta * targetExprs[j];       // predictions
        fout << endl;

	#if !SHORT_OUTPUT
        // print the agreement bewtween predictions and observations
        cout << seqNames[i] << "\t" << beta << "\t"; 
        if ( ExprPredictor::objOption == SSE ) 
            cout << error << endl;       
	else if ( ExprPredictor::objOption == SSE_V ) 
            cout << least_square_variance( targetExprs, observedExprs, beta, 1.0 ) << endl;
        else if ( ExprPredictor::objOption == CORR ) 
            cout << corr( targetExprs, observedExprs ) << endl;
        else if ( ExprPredictor::objOption == CROSS_CORR )
            cout << ExprPredictor::exprSimCrossCorr( observedExprs, targetExprs ) << endl; 
        else if ( ExprPredictor::objOption == NORM_CORR )
            cout << norm_corr( observedExprs, targetExprs ) << endl; 
	#endif // SHORT_OUTPUT
    }
    #endif // CROSS_VALIDATION

    delete predictor;

    int dmm = 0;
    int dss = 0;
    long time_end = time(NULL) - time_start; // gives the time elapsed since time_start in seconds
    dss =time_end % 60; // the remainder is seconds to be displayed
    int minutes= time_end / 60;  // the total minutes in float
    dmm= minutes % 60;  // the remainder are minutes to be displayed
    int dhh = minutes / 60 ; // the total hours in float


    string string_dmm;
    string string_dss;
    if(dmm < 10)    
	string_dmm = ":0";
    else
	string_dmm = ":";
    if(dss < 10)    
	string_dss = ":0";		
    else
	string_dss = ":";

    cout << "Runtime: " << dhh << string_dmm << dmm << string_dss << dss << endl;    


    return 0;	
}

