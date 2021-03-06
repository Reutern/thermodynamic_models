#include <sstream>
#include "SeqAnnotator.h"
#include "param.h" 

vector<double> energy_balance {0,
0,
0,
0.630233,
0,
0.459308,
1.26417,
0,
2.43632,
1.26638,
0.185422,
0.681974,
0,
1.16636,
1.79246,
2.27374};


bool isNt( int a )
{
    if ( a < 0 || a > 3 ) return false;
    else return true;	
}

int complement( int a )
{
    assert( a >= 0 && a < ALPHABET_SIZE );
            
    if ( a == 0 ) return 3;
    if ( a == 1 ) return 2;
    if ( a == 2 ) return 1;
    if ( a == 3 ) return 0;	
    if ( a == MISSING ) return MISSING;
    if ( a == GAP ) return GAP;	
}

int symbolToInt( char c )
{
    char upper = toupper( c );
    for ( int i = 0; i < ALPHABET_SIZE; i++ ) {
        if ( ALPHABET[ i ] == upper ) return i;	
    }
    
    return -1;
}

char strand2char( bool strand )
{
    if ( strand ) return '+';
    else return '-';	
}

bool char2strand( char c )
{
    assert( c == '+' || c == '-' );
    
    if ( c == '+' ) return true;
    else return false;
}

vector< double > createNtDistr( double gcContent )
{
    assert( gcContent >= 0 && gcContent <= 1.0 );
    
    vector< double > freqs( 4 );
    freqs[0] = ( 1.0 - gcContent ) / 2.0;
    freqs[1] = gcContent / 2.0;
    freqs[2] = freqs[1];
    freqs[3] = freqs[0];

    return freqs;
}

Sequence::Sequence( const string& str )
{
    nts.resize(str.size());
    accessibility.resize(str.size());
    for ( int i = 0; i < str.size(); i++ ) {
        int nt = symbolToInt( str[ i ] );	// could be a NNN or gap
        if ( nt >= 0 && nt < ALPHABET_SIZE ) {
            nts[i] = nt ;
	    accessibility[i] = 1.0;
        } else {
            cerr << "Illegal symbol: " << nt << " in " << str << endl;
            exit( 0 );	
        }       
    }
}

Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
    assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );	

    nts.resize(length);
    accessibility.resize(length);
    for ( int i = 0; i < length; i++ ) {
        if ( strand ) {	nts[i] = other[ start + i ] ; accessibility[i] = other.get_accessibility(i); }
        else { nts[i] = complement( other[ start + length - 1 - i ] ) ; accessibility[i] = other.get_accessibility( start + length - 1 - i ); }
    }	
}

int Sequence::push_back( int nt, double ac )
{	
    assert( nt >= 0 && nt < ALPHABET_SIZE );
    nts.push_back( nt );
    accessibility.push_back( ac );
    return 0;
}

int Sequence::push_back( const Sequence& elem )
{
    for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ], elem.accessibility[ i ] );	
    return 0;
}

Sequence Sequence::compRevCompl() const
{
    return Sequence( *this, 0, size(), false );	
}

void Sequence::getNtCounts( vector< int >& counts ) const
{
    counts.clear();
    for ( int i = 0; i < NBASES; i++ ) {
        counts.push_back( 0 );
    }
    
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;	
    }
}

bool Sequence::containsMissing() const
{
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] == MISSING ) return true;	
    }	
    
    return false;
}

int Sequence::load( const string& file, const string& accfile, string& name, int format )
{
    vector< Sequence > seqs;
    vector< string > names;
    int rval = readSequences( file, accfile, seqs, names, format );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( seqs[ 0 ] );
    name = names[ 0 ];
    return rval;
}

int Sequence::load( const string& file,const string& accfile, int format )
{
    string name;
    int rval = load( file, accfile, name, format );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Sequence& seq )
{
    // output the nts
    for ( int i = 0; i < seq.size(); i++ ) {
        os << ALPHABET[ seq[ i ] ];	
    }	
                    
    return os;
}

int readSequences( const string& file, const string& accfile, vector< Sequence >& seqs, vector< string >& names, int format )
{
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }
    seqs.clear();
    names.clear();

    // 	open the files
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

    ifstream facc( accfile.c_str() );
    if ( !facc ) { cerr << "Cannot open" << accfile << endl; exit( 1 ); }

    string line;
    string line_acc;
    Sequence seq;
    
    // read sequences: FASTA format
    if ( format == FASTA ) {
		bool rlin = static_cast<bool> (getline( fin, line ));

        while ( rlin ) {
	        // add the sequence and start a new sequence if the line starts with >
            //cout << line << endl;
            if ( line[ 0 ] == '>' ) { 
                if ( seq.size() ) {
                    seqs.push_back( seq );
                    seq.clear();	
                }
                        
                stringstream ss( line.substr( 1 ) );
                string name; 
                ss >> name;
                names.push_back( name );
				rlin = static_cast<bool>(getline( fin, line ));
            } else { 
				string line_tmp;
			    do{
            		if ( line[ 0 ] == '>' ) { break; }
					line_tmp += line;
					rlin = static_cast<bool> (getline( fin, line ));
                 } while( rlin ); 

		         // check if the line contains content
		         int start = line_tmp.find_first_not_of( " \t\r" );
		         int last = line_tmp.find_last_not_of( " \t\r" );
		         if ( start == string::npos || last == string::npos ) continue;  

                 // append the sequence
		double ac = 1.0;
	    getline( facc, line_acc );
        for ( int i = start; i < last; i++ ) {
                int nt = symbolToInt( line_tmp[ i ] );	// could be a NNN or gap
		        getline( facc, line_acc );
			    //ac = std::stod (line_acc);	// read in according accessibility
 	            std::size_t pos = line_acc.find("\t");      // position of next gap in line_acc
		    	line_acc = line_acc.substr(pos+1);   	// Truncate line_acc for next read-out
			    ac = std::stod (line_acc);	// read in according accessibility
    	        if ( nt >= 0 && nt < ALPHABET_SIZE ) {   
                	seq.push_back( nt, ac ); 	
      	        } 
				else {          
         	    	//cerr << "Illegal symbol: " << nt << " in " << file << endl;
         	        return RET_ERROR;	
      	     	}           
         	   } 
   	    getline( facc, line_acc );
        	}			
    	}
        // add the last sequence
        if( seq.size() ) seqs.push_back( seq );
                        
        return 0;
    }	
}

int readSequences( const string& file, const string& accfile, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, accfile, seqs, names, format );	
    return rval;
}

int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }
    seqs.clear();
    names.clear();
     
    // 	open the files
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

    string line;
    Sequence seq;
    
    // read sequences: FASTA format
    if ( format == FASTA ) {
        while ( getline( fin, line ) ) {
	    

            // add the sequence and start a new sequence if the line starts with >
            //cout << line << endl;
            if ( line[ 0 ] == '>' ) { 	
                if ( seq.size() ) {
                    seqs.push_back( seq );
                    seq.clear();	
                }
                        
                stringstream ss( line.substr( 1 ) );
                string name; 
                ss >> name;
                names.push_back( name );
            } else { 
                // check if the line contains content
                int start = line.find_first_not_of( " \t\r" );
                int last = line.find_last_not_of( " \t\r" );
                if ( start == string::npos || last == string::npos ) continue;
                        

                for ( int i = start; i <= last; i++ ) {
                    int nt = symbolToInt( line[ i ] );	// could be a NNN or gap
                    if ( nt >= 0 && nt < ALPHABET_SIZE ) {
                        seq.push_back( nt, 1.0 ); 	// seq.push_back( nt, ac );
                    } else {
                        //cerr << "Illegal symbol: " << nt << " in " << file << endl;
                        return RET_ERROR;	
                    } 
                }
            }			
        }
            
        // add the last sequence
        if( seq.size() ) seqs.push_back( seq );
                        
        return 0;
    }	
}

int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );	
    return rval;
}


int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
    assert( seqs.size() == names.size() );
    
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }
            
    ofstream fout( file.c_str() );
    
    if ( format == FASTA ) {
        for ( int i = 0; i < seqs.size(); i++ ) {
            fout << ">" << names[ i ] << endl;
            fout << seqs[ i ] << endl;
        }
    }
    
    return 0;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
    // default name: integer starting from 1
    vector< string > names;
    for ( int i = 0; i < seqs.size(); i++ ) {
        char buffer[ 10 ];
        sprintf( buffer, "%i", i );
        names.push_back( string( buffer ) );	
    }	
    
    // print
    return writeSequences( file, seqs, names, format );
}

Matrix compWtmx( const Matrix& countMatrix, double pseudoCount )
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    
    int l = countMatrix.nRows();		// l: the length of motif		
    Matrix pwm( l, 4 );

//     // the sum of each position/column should be a const. (number of sequences)
//     double n = 0;		// number of sites used in the count matrix
//     for ( int j = 0; j < 4; j++ ) {
//         n += countMatrix( 0, j );
//     }
//     for ( int i = 1; i < l; i++ ) {
//         double count = 0;
//         for ( int j = 0; j < 4; j++ ) {
//             count += countMatrix( i, j );
//         }
//         if ( count != n ) { cout << "count matrix incorrect" << endl; exit( 1 ); }
//     }
    
    // the multinomial distribution at each column
    for ( int i = 0; i < l; i++ ) {
        double n = 0;       // total counts at this position 
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );	// pseudoCount ensures a minimum possibility of a nucleotide at every position 
        }	
    }
    return pwm;		
}


Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 )
{
    assert( background.size() == 4 );	
    
    init();
}

Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 )
{
    assert( background.size() == 4 );
    
    pwm = compWtmx( countMatrix, pseudoCount );
    init();
}

double Motif::LLR( const Sequence& elem ) const
{
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;
    
    double result = 0;
    for ( int i = 0; i < l; i++ ) {
        result += LLRMat( i, elem[ i ] ); 	
    }
    
    return result;
}

double Motif::energy( const Sequence& elem ) const
{
	return ( -LLR( elem ) + maxLLR );	
}

void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
    assert( rng != NULL );
    
    int l = pwm.nRows();
    Sequence sampleElem;
    for ( int i = 0; i < l; i++ ) {
        // nt. distribution at position i
        vector< double > distr = pwm.getRow( i );
        
        // sample nt. from this distribution	
        int nt = sampleMul( rng, distr );
        sampleElem.push_back( nt, 1.0 );
    }		
    
    if ( strand == 0 ) elem = sampleElem.compRevCompl();
    else elem = sampleElem;
}

int Motif::load( const string& file, const vector< double >& background, string& name )
{
    vector< Motif > motifs;
    vector< string > names;
    int rval = readMotifs( file, background, motifs, names );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( motifs[ 0 ] );
    name = names[ 0 ];
    return rval;				
}

int Motif::load( const string& file, const vector< double >& background )
{
    string name;
    int rval = load( file, background, name );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Motif& motif )
{
    os << motif.pwm;
    
    return os;
}

void Motif::init()
{
    int l = pwm.nRows();
    
    // compute the LLR matrix
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat( i, j ) = log( pwm( i, j ) / background[ j ] );
        }
    }
    
    // the strongest site
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ) , b_max );
        maxSite.push_back( b_max, 1.0 );	
    }
    
    // compute the LLR of the strongest site
    maxLLR = 0;
    for ( int i = 0; i < l; i++ ) {
        maxLLR += LLRMat( i, maxSite[ i ] );	
    }
}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
    // 	open the file
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }
    motifs.clear(); 
    names.clear();

    string line;
    
    // read the motifs
    do {
        getline( fin, line );
        
        if ( line[ 0 ] != '>' ) continue;
        
        // read the names, length and pseudocount
        int MAX_SIZE = 100;
        char lineStr[ MAX_SIZE ];
        strcpy( lineStr, ( line.substr( 1 ) ).c_str() );
        char *name, *lengthStr, *pseudoCountStr;
        name = strtok( lineStr, " \t" );
        lengthStr = strtok( NULL, " \t" );
        pseudoCountStr = strtok( NULL, " \t PSEUDO_COUNT ");
        int length;
        double pseudoCount;
        if ( lengthStr ) length = atoi( lengthStr );
        else { return RET_ERROR; }
        if ( pseudoCountStr ) pseudoCount = atof( pseudoCountStr ); 
        else pseudoCount = PSEUDO_COUNT;
        
	// read the count matrix
        Matrix countMat( length, NBASES );
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < NBASES; ++j ) {
                fin >> countMat( i, j );
            }	
        }
        
        // create the motif
        names.push_back( string( name ) );
        motifs.push_back( Motif( countMat, pseudoCount, background ) );	
    } while ( !fin.eof() );
                                    
    return 0;
}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
    vector< string > names;
    return readMotifs( file, background, motifs, names );	
}

ostream& operator<<( ostream& os, const Site& site )
{
    char strandChar = site.strand ? '+' : '-';
    os << site.start + 1 << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;
    
    return os;
}


int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names)
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    sites.clear();
    names.clear();

    SiteVec currVec;
    int nrecords = 0;       // number of ">" read so far
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;

        if ( line.substr( 0, 1 ) == ">" ) {
            stringstream ss( line.substr( 1 ) );
            string name; 
            ss >> name;
            names.push_back( name );
            nrecords++;
            if ( nrecords > 1 ) {
                sites.push_back( currVec );
                currVec.clear();
            }
        } else {
            int start;
            char strandChar;
            string factor;
            double energy = 0;
            stringstream ss( line );
            ss >> start >> strandChar >> factor >> energy; 
            bool strand = strandChar == '+' ? 1 : 0;
            map<string, int>::const_iterator iter = factorIdxMap.find( factor );
            currVec.push_back( Site( start, strand, iter->second , energy, 1.0 ) );		// TODO: read in accessibility
        }
    }

    sites.push_back( currVec );

    return 0;
}

int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites)
{
    vector< string > names;
    return readSites( file, factorIdxMap, sites, names);
}


int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
	//cout << "start annotation:" << endl;
    	sites.clear();
    
    	// scan the sequence for the sites of all motifs
    	for ( int i = 0; i < seq.size(); i++ ) {
        // test for each motif
		    for ( int k = 0; k < motifs.size(); k++ ) {
                int l = motifs[ k ].length();
            	if ( i + l > seq.size() ) continue;
               	double energy;
			    double access;
                
                // positive strand
                Sequence elem( seq, i, l, 1 );
			    access = elem.get_accessibility();
                energy = motifs[ k ].energy( elem );
			    if ( energy <=  energyThrFactors[ k ] * motifs[ k ].getMaxLLR() ) {		// accessibility is incorporated in the decision if something is a bninding site
				    //cout << "Energy Diff for motif: " << k << " = " << energy << "\t";
				    //cout << elem << Site( i, 1, k,  energy ) << endl;	      
                    #if APPLY_ENERGY_BALANCE  
                    energy = energy - energy_balance[k]; 
                    #endif //APPLY_ENERGY_BALANCE
		            sites.push_back( Site( i, 1, k, energy, access ) );
                }	
                
                // negative strand
                Sequence rcElem( seq, i, l, 0 );
                energy = motifs[ k ].energy( rcElem );
			    access = rcElem.get_accessibility();
			    if ( energy <= energyThrFactors[ k ] * motifs[k].getMaxLLR() ) {		// accessibility is incorporated in the decision if something is a binding site
				    //cout << "Energy Diff for motif: " << k << " = " << energy << "\t";
		     	   	//cout << rcElem << Site( i, 0, k,  energy ) << endl;		
                    #if APPLY_ENERGY_BALANCE  
                    energy = energy - energy_balance[k]; 
                    #endif //APPLY_ENERGY_BALANCE	
                	sites.push_back( Site( i, 0, k, energy, access ) );
                }				
        	}	
    	}
    
 	//cout << "end annotation" << endl;
    	return sites.size();
}


int SeqAnnotator::compEnergy( const Sequence& seq, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ ) {
        Sequence elem( seq, sites[i].start, motifs[sites[i].factorIdx].length(), sites[i].strand );
        sites[i].energy = motifs[sites[i].factorIdx].energy( elem );
        sites[i].wtRatio = exp( -sites[i].energy );
    }

    return 0;
}

int SeqAnnotator::printEnergy( ostream& os, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ ) {
	os << sites[i].factorIdx << "\t" << sites[i].wtRatio << "\t" << sites[i].energy << endl; 
    }

    return 0;
}

