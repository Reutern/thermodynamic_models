#include "ExprPredictor.h"
#include "OccPredictor.h"
#include <fstream>
#include <time.h>
#include <math.h> // for fmod
#include "param.h" 

int main( int argc, char* argv[] ) 
{

    // Set the timer
    long time_start = time(NULL);
  
    // command line processing
    string seqFile, test_seqFile, accFile, test_accFile, annFile, exprFile, test_exprFile, motifFile, factorExprFile, coopFile, SynFile, factorInfoFile, repressionFile, parFile, print_parFile;
    string outFile, occFile, impactFile;     // output files
    int coopDistThr = 150;
    int SynDistThr = 50;
    double factorIntSigma = 25.0;   // sigma parameter for the Gaussian interaction function
    int repressionDistThr = 0;
    int maxContact = 1;
	double hyperparameter = 0.1;
	vector<double> eTF (20,0.6);


	string free_fix_indicator_filename;
	ExprPredictor::one_qbtm_per_crm = ONE_QBTM;
	ExprPar::one_qbtm_per_crm = ONE_QBTM;
	ExprFunc::one_qbtm_per_crm = ONE_QBTM;

    for ( int i = 1; i < argc; i++ ) {
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
        else if ( !strcmp( "-ts", argv[ i ] ) )
            test_seqFile = argv[ ++i ];
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];      
        else if ( !strcmp( "-acc", argv[ i ] ) ){
            accFile = argv[ ++i ];}
        else if ( !strcmp( "-tacc", argv[ i ] ) )
            test_accFile = argv[ ++i ];      
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
        else if ( !strcmp( "-sy", argv[ i ] ) )
            SynFile = argv[ ++i ];
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
        else if ( !strcmp( "-ct", argv[i] ) )
            coopDistThr = atof( argv[++i] ); 
        else if ( !strcmp( "-syt", argv[i] ) )
            SynDistThr = atof( argv[++i] ); 
        else if ( !strcmp( "-sigma", argv[i] ) )
            factorIntSigma = atof( argv[++i] );    
        else if ( !strcmp( "-rt", argv[i] ) )
            repressionDistThr = atof( argv[++i] );    
        else if ( !strcmp( "-na", argv[i] ) )
            ExprPredictor::nAlternations = atoi( argv[++i] );    
        else if ( !strcmp( "-ff", argv[i] ) )
            free_fix_indicator_filename = argv[++i];
        else if ( !strcmp( "-if", argv[i] ) )
            impactFile = argv[++i];    
        else if ( !strcmp( "-hy", argv[i] ) )
            hyperparameter = atof( argv[++i] );    
        else if ( !strcmp( "-oq", argv[i] ) ){
           	ExprPredictor::one_qbtm_per_crm = true;    
			ExprPar::one_qbtm_per_crm = true;
			ExprFunc::one_qbtm_per_crm = true;
		}
        else if ( !strcmp( "-et", argv[i] ) ) {}
           // eTF = atof( argv[ ++i ] );    
    }
    if (seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || 
	   ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) 
	   &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
    	cerr << "Usage: " << argv[ 0 ] << " -s seqFile -ts test_seqFile -e exprFile -te test_exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile  -c SynFile -i	factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -pp print_parFile -rt repressionDistThr -na nAlternations -ct coopDistThr -ct SynDistThr -hy hyperparameter -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }           

    // additional control parameters
    double gcContent = 0.5;
    FactorIntType FactorIntOption = FactorIntFunc;     // type of interaction function
    ExprPar::searchOption = CONSTRAINED;      // search option: unconstrained; constrained. 
    ExprPar::estBindingOption = 1;
	ExprPar::default_par_penalty = hyperparameter;

    ExprPredictor::nRandStarts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_PGP = 1.0E-10;
    ExprPredictor::nSimplexIters = 2000;
    ExprPredictor::nCMAESIters = 10000;
    ExprPredictor::nGradientIters = 200;

    int rval;
    vector< vector< double > > data;    // buffer for reading matrix data
    vector< string > labels;    // buffer for reading the labels of matrix data
    string factor1, factor2;    // buffer for reading factor pairs

    // read the sequences
	cout << "Read sequence" << endl;
    vector< Sequence > seqs;
    vector< string > seqNames;

    #if ACCESSIBILITY
    rval = readSequences( seqFile, accFile, seqs, seqNames );
    #else
    rval = readSequences( seqFile, seqs, seqNames );
    #endif // ACCESSIBILITY

    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();

    // read the expression data
	cout << "Read expression" << endl;
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    for ( int i = 0; i < nSeqs; i++ ){
    	if( labels[ i ] != seqNames[ i ] ) cout << labels[i] << seqNames[i] <<  endl;
    }
    Matrix exprData( data ); 
    int nConds = exprData.nCols();

    // read the motifs
	cout << "Read motifs" << endl;
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
	cout << "Read TF concentrations" << endl;
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
		energyThrFactors.push_back( eTF[index] );
	}

    // site representation of the sequences
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrFactors );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqs; i++ ) {
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

/*
	//Delete every cic site that is not in 15n distance to a second cic site
	
	for ( int seq_idx = 0; seq_idx < nSeqs; seq_idx++ ) {
		SiteVec sites_tmp = seqSites[seq_idx];
		for (int site_idx = seqLengths[seq_idx]; site_idx > 0; site_idx--){		
			if( seqSites[seq_idx][site_idx].factorIdx != 7 ) continue;
			bool delete_site = true;
			int position_1 = seqSites[seq_idx][site_idx].start;
	 		for (int site_idx_2 = 0; site_idx_2 < seqLengths[seq_idx]; site_idx_2++){
				if( seqSites[seq_idx][site_idx_2].factorIdx != 7 ) continue;	
				int position_2 = seqSites[seq_idx][site_idx_2].start;	
				if (abs(position_1 - position_2) == 15) {
					cout << position_1 << " " << position_2 << endl;
					delete_site = false;
					break;
				}
			}
			if (delete_site == true)
				sites_tmp.erase(sites_tmp.begin() + site_idx);		
		}
		seqSites[seq_idx] = sites_tmp;
	}
*/
    // read the cooperativity matrix 
    int num_of_coop_pairs = 0;
    IntMatrix coopMat( nFactors, nFactors, false );
    if (!coopFile.empty()) {
        ifstream fcoop( coopFile.c_str() );
        if (!fcoop) {
            cerr << "Cannot open the cooperativity file " << coopFile << endl;
            exit( 1 );
        }  
        while (fcoop >> factor1 >> factor2) {
			if(factor1 == "all"){
				IntMatrix fullMat( nFactors, nFactors, true );
				coopMat = fullMat;
				num_of_coop_pairs = (nFactors + 1) * nFactors / 2;
			}
			else{
		        assert( factorIdxMap.count(factor1) && factorIdxMap.count(factor2) );
		        int idx1 = factorIdxMap[factor1];
		        int idx2 = factorIdxMap[factor2];
		        if(coopMat(idx1, idx2) == false && coopMat(idx2, idx1) == false ){
					coopMat(idx1, idx2) = true;
					coopMat(idx2, idx1) = true;
					num_of_coop_pairs ++;
				}
	    	}
        }        
    } 

    // read the synergy matrix 
    int num_of_Syn_pairs = 0;
    IntMatrix SynMat( nFactors, nFactors, false );
    if (!SynFile.empty()) {
        ifstream fSyn( SynFile.c_str() );
        if (!fSyn) {
            cerr << "Cannot open the Synergy file " << SynFile << endl;
            exit( 1 );
        }  
        while (fSyn >> factor1 >> factor2) {
			if(factor1 == "all"){
				IntMatrix fullMat( nFactors, nFactors, true );
				SynMat = fullMat;
				num_of_Syn_pairs = (nFactors + 1) * nFactors / 2;
			}
			else{
		        assert( factorIdxMap.count(factor1) && factorIdxMap.count(factor2) );
		        int idx1 = factorIdxMap[factor1];
		        int idx2 = factorIdxMap[factor2];
		        if(SynMat(idx1, idx2) == false && SynMat(idx2, idx1) == false ){
					SynMat(idx1, idx2) = true;
					SynMat(idx2, idx1) = true;
					num_of_Syn_pairs ++;
				}
	    	}
        }        
    } 

    // read the roles of factors
    vector< bool > actIndicators(nFactors, false);
    vector< bool > repIndicators(nFactors, false);
    if (!factorInfoFile.empty()) {
        ifstream finfo(factorInfoFile.c_str());
        if (!finfo) {
            cerr << "Cannot open the factor information file " << factorInfoFile << endl;
            exit(1);
        }      
        string name;
        int i = 0, actRole, repRole;
        while (finfo >> name >> actRole >> repRole) {
            assert(name == motifNames[i]);
            if (actRole) actIndicators[i] = true;
            if (repRole) repIndicators[i] = true;
            i++;
        }
    }

    // read the repression matrix 
    IntMatrix repressionMat(nFactors, nFactors, false);
    if (!repressionFile.empty()) {
        ifstream frepr(repressionFile.c_str());
        if (!frepr) {
            cerr << "Cannot open the repression file " << repressionFile << endl;
            exit(1);
        }        
        while (frepr >> factor1 >> factor2) {
            assert(factorIdxMap.count(factor1) && factorIdxMap.count(factor2));
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            repressionMat(idx1, idx2) = true;
        }        
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
		for( int index = 0; index < nFactors + num_of_coop_pairs + nFactors ; index++ ){
			indicator_bool.push_back( true );
		}
		for( int index = 0; index < nFactors + num_of_Syn_pairs + nFactors ; index++ ){
			indicator_bool.push_back( true );
		}
		if( ExprPredictor::one_qbtm_per_crm ){
			for( int index = 0; index < nSeqs; index++ ){
				indicator_bool.push_back( false );
			}
		}
		else{
			indicator_bool.push_back( true );
		}
		#if ACCESSIBILITY
		indicator_bool.push_back( true ); 	// acc_scale
		#endif //ACCESSIBILITY
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
    cout << "Penalty_Function = " << getPenaltyOptionStr( ExprPredictor::PenaltyOption ) << "\tlambda = " << hyperparameter << endl;
    if ( !coopFile.empty() ) {
        cout << "Interaction_Model = " << getIntOptionStr( FactorIntOption ) << endl;
        cout << "Interaction_Distance_Threshold = " << coopDistThr << endl;
        if ( FactorIntOption == GAUSSIAN ) cout << "Sigma = " << factorIntSigma << endl;
    }
    cout << "Search_Option = " << getSearchOptionStr( ExprPar::searchOption ) << endl;
    cout << "Num_Random_Starts = " << ExprPredictor::nRandStarts << endl;
    cout << "Energy Threshold Factor = " << eTF << endl;
    cout << "Pseudo count = " << PSEUDO_COUNT << endl;
    cout << endl;
	#if PRINT_STATISTICS
    ofstream weights_file;
    weights_file.open ("../data/weights_stark.txt");
    cout << "Statistics: " << endl; 
    cout << "Factors "<< nFactors << "\t " << "Sequences " << nSeqs <<  endl;
    for(int motif_idx = 0; motif_idx < nFactors; motif_idx++){
   		weights_file << motifNames[ motif_idx ] << " \t ";}
	weights_file << "\n";
    cout << "Sum\t Name\t Length" << endl;
    double average_number = 0;
    for(int seqs_idx = 0; seqs_idx < nSeqs; seqs_idx++){
		weights_file << seqNames[seqs_idx] << " \t ";
		average_number += seqSites[seqs_idx].size()/nSeqs;
		vector<double> weight_count (nFactors,0);
		for( int idx = 1; idx < seqSites[seqs_idx].size() ; idx++ ){
				#if ACCESSIBILITY
				double weight_tmp = (1-seqSites[seqs_idx][idx].accessibility )* seqSites[seqs_idx][idx].wtRatio;
				#else
				double weight_tmp = seqSites[seqs_idx][idx].wtRatio;
				#endif // ACCESSIBILITY
				weight_count[seqSites[seqs_idx][idx].factorIdx] = weight_count[seqSites[seqs_idx][idx].factorIdx] + weight_tmp;
		}
		for( int l = 0; l < nFactors; l++){
			cout << weight_count[l]  << " \t "; }	
		cout << seqSites[seqs_idx].size() - 1 << " \t " << seqNames[seqs_idx] << " \t " << seqLengths[seqs_idx] <<  endl;
	}

    cout << average_number << endl;
    #endif // PRINT_STATISTICS

    #if SAVE_ENERGIES
    cout << "Save site energies" << endl;
    ofstream site_energies;
    site_energies.open ("../data/sites/site_weights.txt");
    for(int seqs_idx = 0; seqs_idx < nSeqs; seqs_idx++){
		site_energies << seqNames[seqs_idx] << "\t Nsites=" << seqSites[seqs_idx].size()-1 << "\t len=" << seqLengths[seqs_idx] << endl;
        for(int sites_idx = 1; sites_idx < seqSites[seqs_idx].size(); sites_idx++){
			int idx = seqSites[seqs_idx][sites_idx].factorIdx;	
			#if ACCESSIBILITY
			double weight_tmp = (1-seqSites[seqs_idx][sites_idx].accessibility )* seqSites[seqs_idx][sites_idx].wtRatio;
			#else
			double weight_tmp = seqSites[seqs_idx][sites_idx].wtRatio;
			#endif // ACCESSIBILITY
			site_energies << seqSites[seqs_idx][sites_idx].start << "\t"
						  << seqSites[seqs_idx][sites_idx].start + static_cast<int> (motifs[idx].length()) << "\t" 
						  << seqSites[seqs_idx][sites_idx].strand << "\t"
						  << motifNames[idx] << "\t" 
						  << weight_tmp << "\t" << endl;
		}
    }
    // Write site energies to site_energies.txt	
    site_energies.close();
    cout << "Site energies saved" << endl;
    #endif //SAVE_ENERGIES

    // read the initial parameter values
    ExprPar par_init( nFactors, nSeqs );
    if ( !parFile.empty() ) {
        rval = par_init.load( parFile, seqNames, motifNames);
        if ( rval == RET_ERROR ) {
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        } 
    }
    // Initialise the predictor class
    ExprPredictor* predictor = new ExprPredictor( seqSites, seqLengths, exprData, motifs, factorExprData, coopMat, SynMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, SynDistThr, indicator_bool, motifNames, seqNames, seqs );

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
        par.print( cout, motifNames, seqNames, coopMat, SynMat );
    }
    else {
        par.print( pout, motifNames, seqNames, coopMat, SynMat );
    }
    cout << "Performance = " << setprecision( 5 ) << ( ExprPredictor::objOption == SSE ? predictor->getObj() : -predictor->getObj() ) << endl;

    // print the predictions
    ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }
    fout << "Rows\t" << condNames << endl;

	#if ACCESSIBILITY
	double impact_acc = predictor->comp_impact_acc(par);
	cout << "Impact accessibility: " << impact_acc << endl;
	
	#endif // ACCESSIBILITY


    #if CALCULATE_OCCUPANCY
    ofstream Occout( occFile.c_str() );
    for ( int seq_idx = 0; seq_idx < nSeqs; seq_idx++ ) {
		printf( "\rOccupancy prediction for sequence: %i ", seq_idx+1);
		fflush(stdout);

    	// Initialise the occupancy predictor
    	OccPredictor* occpred = new OccPredictor( seqSites[seq_idx], motifs, factorExprData, coopMat, SynMat, coopDistThr, SynDistThr, par_init );
		Occout << seqNames[seq_idx] << "\t" << seqSites[seq_idx].size() - 1 << "\t" << seqLengths[seq_idx] << "\t" << "strand" << endl ; 

    	for(int site_idx = 1; site_idx < seqSites[seq_idx].size(); site_idx++){
			double occ_avg = 0;
			for(int position_idx = 0; position_idx < 100; position_idx++){ 
				occ_avg += occpred -> predictOcc(site_idx, position_idx) / 100.0 ;
			}
			Occout << seqSites[seq_idx][site_idx].start << "\t" << seqSites[seq_idx][site_idx].factorIdx <<  "\t" << occ_avg <<  "\t" << seqSites[seq_idx][site_idx].strand << endl;
		}
    	delete occpred;		// clean up the predictor
    }

    cout << endl;
    #endif // CALCULATE_OCCUPANCY
    #if CROSS_VALIDATION
    // read the sequences
    vector< Sequence > test_seqs;
    vector< string > test_seqNames;
    #if ACCESSIBILITY
    rval = readSequences( test_seqFile, test_accFile, test_seqs, test_seqNames );
    #else
    rval = readSequences( test_seqFile, test_seqs, test_seqNames );
    #endif // ACCESSIBILITY
    assert( rval != RET_ERROR );
    int test_nSeqs = test_seqs.size();

    // read the expression data
    vector< vector< double > > test_data;
    rval = readMatrix( test_exprFile, labels, condNames, test_data );
    assert( rval != RET_ERROR );
    for ( int i = 0; i < test_nSeqs; i++ ){
    	if( labels[ i ] != test_seqNames[ i ] ) cout << labels[i] << test_seqNames[i] << endl;
    }

    Matrix test_exprData( test_data ); 
    // site representation of the sequences	
    vector< SiteVec > test_seqSites( test_nSeqs );
    vector< int > test_seqLengths( test_nSeqs );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < test_nSeqs; i++ ) {
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

    // Additional parameter training for the qbtm of the test enhancers
    if(ExprPredictor::one_qbtm_per_crm){
		vector <double> basalTxps_modified;
		basalTxps_modified.clear();
		for( int index = 0; index < test_nSeqs; index++ ){
			basalTxps_modified.push_back( 0.1 );
		}
		// Modify parameters
		par.basalTxps = basalTxps_modified;
    }
	// New predictor 
	ExprPredictor* predictor_CV = new ExprPredictor( test_seqSites, test_seqLengths, test_exprData, motifs, factorExprData, coopMat, SynMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, coopDistThr, SynDistThr, indicator_bool, motifNames, test_seqNames, test_seqs );
	predictor_CV->setPar(par);


    // random number generator
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)

    double squaredErr = 0;
    double correlation = 0;
    for ( int i = 0; i < test_nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor_CV->predict( test_seqSites[i], test_seqLengths[i], targetExprs, i );
        vector< double > observedExprs = test_exprData.getRow( i );
        // print the results
        fout << test_seqNames[i] << "\t" << observedExprs << endl;      // observations
        fout << test_seqNames[i];

        for ( int j = 0; j < nConds; j++ ) fout << "\t" << targetExprs[j];       // predictions
        fout << endl;

        double beta = 1.0;
        squaredErr +=  least_square( targetExprs, observedExprs, beta );
		correlation += corr( targetExprs, observedExprs ); 
    }	

    double obj_corr = correlation / test_nSeqs;
    double obj_sse = sqrt( squaredErr / ( test_nSeqs * nConds ) ); 
    cout << "Performance on test set: SSE = " << obj_sse << "\t" << "Corr = " << obj_corr << endl;

	#if PRINT_IMPACT
    cout << "Save impact" << endl;
	ExprPar par_impact = predictor_CV->getPar();
    ofstream impact_stream;
    impact_stream.open (impactFile);
	impact_stream << "CRM";
    for(int tf = 0; tf < nFactors; tf++ ){
		impact_stream << "\t" << motifNames[tf];}
	impact_stream << endl;
    for(int seqs_idx = 0; seqs_idx < test_nSeqs; seqs_idx++){
		impact_stream << test_seqNames[seqs_idx];
        for(int tf = 0; tf < nFactors; tf++ ){
			double impact = predictor_CV->comp_impact(par_impact,tf,seqs_idx);
			impact_stream << "\t" << impact;
		}
		impact_stream << endl;
    }

	impact_stream << "TF/TF";
    for(int tf = 0; tf < nFactors; tf++ ){
		impact_stream << "\t" << motifNames[tf];}
		impact_stream << endl;
    	for(int tf_1 = 0; tf_1 < nFactors; tf_1++){
			impact_stream << motifNames[tf_1];
			for(int tf_2 = 0; tf_2 < nFactors; tf_2++){
				double impact = 0;			
				if(tf_2 >= tf_1 and coopMat(tf_1, tf_2) == 1 )				
					impact = predictor_CV->comp_impact_coop_pair(par_impact, tf_1, tf_2);
				impact_stream << "\t" << impact;
			}
			impact_stream << endl;
    	}
    impact_stream.close();
    cout << "Impact saved" << endl;

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
	double log_impact = log10( abs( impact ) ); 

	double impact_new = predictor_CV->comp_impact(par,tf) * 100;
	double impact_coop = predictor_CV->comp_impact_coop(par,tf) *100;
	printf ("%s \t %4.2f \t %4.1f %% \t %4.1f %% \n", motifNames[tf].c_str(), log_impact, impact_new, impact_coop);
	//cout << motifNames[tf] << "\t" << log_impact << "\t" << impact_new << "\t" << impact_coop << endl;
    }
	#endif // PRINT_IMPACT


	delete predictor_CV;
    #else // CROSS_VALIDATION
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
        else if ( ExprPredictor::objOption == CORR ) 
            cout << corr( targetExprs, observedExprs ) << endl;
	#endif // SHORT_OUTPUT
    }
    #endif // CROSS_VALIDATION

    delete predictor;

    int dmm = 0;
    int dss = 0;
    long time_end = time(NULL) - time_start; // gives the time elapsed since time_start in seconds
    dss =time_end % 60; // the remainder in seconds to be displayed
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
