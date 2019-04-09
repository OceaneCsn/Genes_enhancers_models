#include <algorithm>
#include <sstream>
#include <libgen.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <limits>

#include "SSA/SSA.h"
#include "GenomeInfo.h"
#include "matrixFile.hpp"
#include "seqGenerator.hpp"
#include "utils.hpp"

#include <math.h>


static void
display_usage( char const * const argv_0 ) {
    std::cout << "USAGE: " << argv_0 << "[OPTIONS] PREFIX_INDEX MATRIX [MATRIX...]" << std::endl ;
    std::cout << "OPTION: "                                                         << std::endl ;
    std::cout << "\t-r, --reverse-complement"                                       << std::endl ;
    std::cout << "\t\tsearch on reverse complement strand"                          << std::endl ;
    std::cout << "\t-h, --homology N"                                               << std::endl ;
    std::cout << "\t\thomology percentage (default 90)"                             << std::endl ;
    std::cout << "\t-t, --threshold "                                               << std::endl ;
    std::cout << "\t\tfor number of potential occurences(default inf)"              << std::endl ;
}


static void
display_relatif_name( GenomeInfo &gInfo, unsigned int const * const poccs, unsigned long int const nb_poccs,
					  const char sequence[], uint8_t const len_sequence, const float score_sequence, const float ratio_sequence ) {
	if (score_sequence < 0) { return; }

	unsigned long int i=0;
	for ( ; i < nb_poccs; ++i) {
			try {
				pair<const unsigned long int, const unsigned long int> rname_and_abs_pos = gInfo.getRelativePositionOfSequence( poccs[ i ], len_sequence ) ;
				std::cout << sequence << '\t' << score_sequence << '\t' << ratio_sequence << std::endl ;
				std::cout << gInfo.getChrName( rname_and_abs_pos.first ) <<  ":" << rname_and_abs_pos.second + 1 << std::endl ;
				break;
			}
			catch (SequenceMismatches & e) {
				continue ;
			}
	}

	for ( ++i; i < nb_poccs; ++i) {
		try {
			pair<const unsigned long int, const unsigned long int> rname_and_abs_pos = gInfo.getRelativePositionOfSequence( poccs[ i ], len_sequence ) ;
			std::cout << gInfo.getChrName( rname_and_abs_pos.first ) <<  ":" << rname_and_abs_pos.second + 1 << std::endl ;
		}
		catch (SequenceMismatches & e) {
			continue ;
		}
	}
}

float
compute_sequence_score( const DnaMatrix pwm, char const * const sequence, const unsigned int sequence_size ) {
	float score = 0 ;
	for (unsigned int i = 0; i < sequence_size; ++i) {
		switch ( sequence[i] ) {
			case 'A': score += pwm[i][0]; break ;
			case 'C': score += pwm[i][1]; break ;
			case 'G': score += pwm[i][2]; break ;
			case 'T': score += pwm[i][3]; break ;
		}
	}
	return score ;
}

static float min( float a, float b )
{
#define MIN(a,b) (((a)<(b))?(a):(b))
	return MIN(a,b) ;
#undef MIN
}

static float max( float a, float b )
{
#define MAX(a,b) (((a)>(b))?(a):(b))
	return MAX(a,b) ;
#undef MAX
}

int
main (int argc, char * const * argv)
{
    int find_reverse_complement = 0 ;
    float homology_percent = 90 ;
    double nb_mot_max = std::numeric_limits<double>::infinity();

    while (1) {
    	int option_index = 0;
        static struct option long_options[] = {
            {"reverse-complement", no_argument, &find_reverse_complement, 'r' },
            {"homology",     required_argument, NULL,                     'h' },
            {"threshold",     required_argument, NULL,                     't' },
            { 0,                   0,           0,                        0   }
        };

        int c = getopt_long( argc, argv, "h:r:t", long_options, &option_index );
        if (c == -1)
            break;

        switch (c) {
			case 'r':
				find_reverse_complement = 1 ;
				break ;

			case 'h':
				homology_percent = atof( optarg ) ;
				break ;

            		case 't':
				nb_mot_max = atof( optarg ) ;
				break ;

			case '?':
				/* getopt_long already printed an error message. */
				std::exit( EXIT_FAILURE ) ;
				break ;
        }
    }
    // Si pas assez d'arguments, lance l'usage
    if (argc - optind < 2 ) {
    	display_usage( argv[0]) ;
        return 1;
    }

    // Chargement de l'index et du GenomeInfo pour convertir les positons absoluer en relatif (ie. 33655 -> chr5:654)
    TFMindex *index;
    load_index( (char *)argv[optind], (void **)&index );
    GenomeInfo gInfo( (std::string( argv[optind++] ) + ".conf").c_str() ) ;

    // get the nucleotide frequencies in the genome from the .conf file
    float background[4] ;
    background[0] = 0.25 ;  // A
    background[1] = 0.25 ;  // C
    background[2] = 0.25 ;  // G
    background[3] = 0.25 ;  // T

    float origin_percent = homology_percent;
    while (optind < argc) {
		std::ifstream ifs( argv[ optind ], std::ios::in | std::ios::binary );
		if (!ifs)
		{
			cerr << "Error reading matrices '" << optind << "'" << std::endl   ;
			exit(1);
		}

		MotifMatrix matrix( ifs ) ;

		DnaMatrix pwm = matrix.get_pwm( background ) ;
		//nb_mot_max = 25000000;
		float threshold;
		float score_max;
		float score_min;
        homology_percent = origin_percent;
		bool done = false;
		while(!done){
		
  		threshold = get_threshold_by_percent( pwm, matrix.get_length(), homology_percent ) ;
      		score_max = get_score_max(pwm, matrix.get_length());
  		score_min = get_score_min(pwm, matrix.get_length());
  
  		threshold = max(threshold, 0);
  		SeqGenerator *sg = new CpSeqGenerator( pwm, matrix.get_length(), threshold ) ;
  		int count = 0;
  		char gs[matrix.get_length() + 1] ;
  		while(sg->has_next()){
  		  sg->get_next_sequence( gs ) ;
  		  count++;
  		  if(count > nb_mot_max){
  		    break;
  		  }
		}
		if(count < nb_mot_max){
		  done = true;
		}
		else{
		  homology_percent +=2.5;
		}
		}
		
		std::cout << ">" << matrix.get_name() << std::endl ;
		SeqGenerator *seq_generator = new CpSeqGenerator( pwm, matrix.get_length(), threshold ) ;
		std::cerr << matrix.get_name() << " " << homology_percent << std::endl;
		// Search and display
		char generated_sequence[matrix.get_length() + 1] ;
		while ( seq_generator->has_next()  ) {
			seq_generator->get_next_sequence( generated_sequence ) ;

			// Récupérer les positions absolues des sequences.
			unsigned int *poccs_sens = NULL ;
			unsigned long int const nb_sens = index->locate( (unsigned char *)generated_sequence, matrix.get_length(), &poccs_sens, 0 ) ;
			float sequence_score = compute_sequence_score( pwm, generated_sequence, matrix.get_length() ) ;

			float sequence_ratio = min((sequence_score - score_min) / (score_max - score_min), 1) ;

			std::sort( poccs_sens, poccs_sens + nb_sens ) ;
			display_relatif_name( gInfo, poccs_sens, nb_sens, generated_sequence, matrix.get_length(), sequence_score, sequence_ratio ) ;

			if ( find_reverse_complement ) {
				unsigned int *poccs_reverse = NULL ;
				reverse_complement( generated_sequence ) ;
				unsigned long int const nb_reverse = index->locate( (unsigned char *)generated_sequence, matrix.get_length(), &poccs_reverse, 0 ) ;

				std::sort( poccs_reverse, poccs_reverse + nb_reverse ) ;
				display_relatif_name( gInfo, poccs_reverse, nb_reverse, generated_sequence, matrix.get_length(), sequence_score, sequence_ratio ) ;
			}
		}

		delete seq_generator ;
		delete pwm ;
		++optind ;
	}
	return 0;
}
