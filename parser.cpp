/* parser.cpp
 *
 * Copyright (C) 2012 Ido Tal <idotal@ieee.org> and Alexander Vardy <avardy@ucsd.edu>
 *
 * This file is part of polarList.
 *
 */

#include "parser.h"
#include <argp.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "myfloat.h"

using namespace std;

const char *argp_program_version = "polar list 1.0";
const char *argp_program_bug_address =  "<idotal@ieee.org>";

enum {
    OPT_BSC = 256,
    OPT_BEC,
    OPT_AWGN,
    OPT_SIGMA_SQUARE,
    OPT_SIGMA,
    OPT_SNR_INFORMATION_BIT,
    OPT_SNR_CODEWORD_BIT,
    /*    OPT_PRINT_FIRST_MISDECODED_INDEX,
        OPT_PRINT_WORD_ERROR_RATE,
        OPT_PRINT_WORD_ERROR_RATE_CRC,
        OPT_PRINT_BIT_ERROR_RATE,
        OPT_PRINT_BIT_ERROR_RATE_CRC,
        OPT_PRINT_LIST_BIT_ERROR_RATE,
        OPT_PRINT_ML_WORD_ERROR_RATE,
        OPT_PRINT_ML_BIT_ERROR_RATE,
        OPT_PRINT_FIRST_MISDECODED_HISTOGRAM,
        OPT_PRINT_ERROR_RATE_AFTER_FIRST_MISDECODED_BIT,
    //    OPT_PRINT_AVERAGETREEWIDTHFOR_LISTERROR_HISTOGRAM,
    //    OPT_PRINT_AVERAGETREEWIDTHFOR_CORRECTLIST_HISTOGRAM,
        OPT_PRINT_MISDECODED_HISTOGRAM,
        OPT_PRINT_LIST_WORD_ERROR_RATE,*/
    OPT_PRINT_COMMAND_LINE,
    OPT_USE_MIN_SUM,
    OPT_CRC32,
    OPT_CRC24,
    OPT_CRC16,
    OPT_CRC15,
    OPT_CRC10,
    OPT_CRC8,
    OPT_RRC,
    OPT_SYSTEMATIC,
    OPT_DYNAMIC_FROZEN,
    OPT_DYNAMIC_FROZEN_SUFFIX_LENGTH,
};

void setCrcArguments(int crcBitsCount_rvalue, int &crcBitsCount_lvalue, struct argp_state *state)
{
    if ( crcBitsCount_lvalue != 0 )
    {
        printf("Error: Only one CRC value allowed");
        argp_usage(state);
    }

    crcBitsCount_lvalue = crcBitsCount_rvalue;

}

int calcKFromR(int r, int m)
{
    int i;
    unsigned int k;

    k = 0;

    for ( i = 0; i < (1<<m); i++ )
    {
            k++;
    }

    return k;
}

bool yesNoParser(const char *option, const char *val, struct argp_state *state, bool defaultValue = true)
{
    if (val == NULL)
    {
        return defaultValue;
    }

    if (strcmp(val,"yes")==0)
    {
        return true;
    }
    if (strcmp(val,"no")==0)
    {
        return false;
    }

    printf("Error: The value of %s should be either yes or no", option);
    argp_usage(state);
    return true; // keep the compiler happy

}

myfloat positiveFloatParser(const char *option, const char *val, struct argp_state *state, bool positive = true)
{
    double tempDouble = -1;


    if ( sscanf(val, "%lf", &tempDouble) != 1)
    {

        if ( positive == true )
        {
            printf("Error: %s  must be a positive float.\n", option);
        }
        else
        {
            printf("Error: %s  must be a float.\n", option);
        }

        argp_usage (state);
    }

    if ( positive == true && tempDouble <= 0 )
    {
        printf("Error: %s  must be a positive float.\n", option);
        argp_usage (state);
    }

    return tempDouble;
}

int positiveIntParser(const char *option, const char *val, struct argp_state *state)
{
    int tempInt = -1;

    if ( sscanf(val, "%d", &tempInt) != 1)
    {
        printf("Error: %s  must be a positive integer.\n", option);
        argp_usage (state);
    }

    if ( tempInt <= 0 )
    {
        printf("Error: %s  must be a positive integer.\n", option);
        argp_usage (state);
    }

    return tempInt;
}

unsigned long unsignedLongParser(const char *option, const char *val, struct argp_state *state)
{
    unsigned long tempInt = -1;

    if ( sscanf(val, "%lu", &tempInt) != 1)
    {
        printf("Error: %s  must be a positive integer.\n", option);
        argp_usage (state);
    }

    if ( tempInt == 0 )
    {
        printf("Error: %s  must be a positive integer.\n", option);
        argp_usage (state);
    }

    return tempInt;
}


static error_t parse_opt_filename (int key, char *arg, struct argp_state *state)
{
    struct Arguments *arguments = (Arguments *) state->input;
    ifstream infile;
    double tempDouble;


    if (arguments->pass == 0)
    {
        switch (key)
        {
        case 'i':
            arguments->infile = arg;
            arguments->pass = 1;
            state->next = 0; // parse the options again

            infile.open(arg);

            if (infile.fail())
            {
                printf("Error: the channel input filename %s does not exist.\n", arg);
                argp_usage (state);

            }


            infile.close();
            break;

        case 'm': // do nothing, since we've looked at this before
        	arguments->mu = positiveFloatParser("-p (--prob)", arg, state);
            break;


        case 'b':
        	int i;
        	int binInd;
        	i = binInd = 0;
        	while (arg[i]) {
        		if (arg[i]=='0')
        			binInd = binInd << 1;
        		else
        			binInd = (binInd << 1) + 1;
        		i++;

        	}
        	arguments->binInd = binInd;
        	arguments->binIndLength = i;
            break;

        case 'u':
        	arguments->degUp = string(arg);
            break;


        case ARGP_KEY_ARG:
            argp_usage(state);
            break;

        case ARGP_KEY_END:
            printf("Error: no -i (--input-filename)\n");
            argp_usage (state);
            break;

        default:
            break;

        }
    }
    else if (arguments->pass == 1)
    {
        switch (key)
        {
        case 'm': // do nothing, since we've looked at this before
               	arguments->mu = positiveFloatParser("-p (--prob)", arg, state);
                   break;

        case 'b':
        	int i;
        	int binInd;
        	i = binInd = 0;
        	while (arg[i]) {
        		if (arg[i]=='0')
        			binInd = binInd << 1;
        		else
        			binInd = (binInd << 1) + 1;
        		i++;

        	}
        	arguments->binInd = binInd;
        	arguments->binIndLength = i;
            break;

        case 'u':
        	arguments->degUp = string(arg);
            break;
 //       case 'k':
 //           arguments->pass = 2;
 //           state->next = 0; // parse the options again
 //
 //           arguments->k = positiveIntParser("-k (--code-dimension)", arg, state);
 //
 //           if (arguments->k > 1<<arguments->m)
 //           {
 //               printf("Error: -k (--code-dimension)  can not be larger than the number of channels: %d.\n", 1<<arguments->m);
 //               argp_usage (state);
 //           }
 //
 //           break;
 //
 //       case 'r':
 //           arguments->pass = 2;
 //           state->next = 0; // parse the options again
 //
 //           arguments->r = positiveIntParser("-r (--reed-muller-r)", arg, state);
 //
 //           if (arguments->r > arguments->m)
 //           {
 //               printf("Error: -r (--reed-muller-r) can not be larger than -m (--reed-muller-m): %d.\n", arguments->m);
 //               argp_usage (state);
 //           }

            // *for now* set k as if there is no crc. Fix this during the next pass
//            arguments->k = calcKFromR(arguments->r, arguments->m);
//            break;

        case ARGP_KEY_ARG:
            argp_usage(state);
            break;

//        case ARGP_KEY_END:
//            printf("Error: no -k (--code-dimension) option specified.\n");
//            argp_usage (state);
//            break;

        default:
            break;
        }
    }
    else if (arguments->pass == 2)
    {
        switch (key)
        {
        case 'i': // do nothing, since we've looked at this before
            break;

        case 'k': // do nothing, since we've looked at this before
            break;

        case 'r': // do nothing, since we've looked at this before
            break;

        case 'm': // do nothing, since we've looked at this before
        	arguments->mu = positiveFloatParser("-p (--prob)", arg, state);
            break;

        case 'b':
        	int i;
        	int binInd;
        	i = binInd = 0;
        	while (arg[i]) {
        		if (arg[i]=='0')
        			binInd = binInd << 1;
        		else
        			binInd = (binInd << 1) + 1;
        		i++;

        	}
        	arguments->binInd = binInd;
        	arguments->binIndLength = i;

        	break;

        case 'u':
        	arguments->degUp = string(arg);
            break;

        case OPT_BSC:
            arguments->bsc = true;
            break;

        case OPT_BEC:
            arguments->bec = true;
            break;

        case OPT_AWGN:
            arguments->awgn = true;
            break;

        case 'p':
            arguments->p = positiveFloatParser("-p (--prob)", arg, state);
            break;

        case OPT_SIGMA_SQUARE:
            arguments->sigmaSquare = positiveFloatParser("--sigma-square", arg, state);
            break;

        case OPT_SIGMA:
            tempDouble = positiveFloatParser("--sigma-square", arg, state, false);
            arguments->sigmaSquare = tempDouble*tempDouble;
            break;

//        case OPT_SNR_INFORMATION_BIT:
//            tempDouble = positiveFloatParser("--SNR-information", arg, state, false);
//            {
//                double rate;
//
//                rate = ((double) arguments->k)/(1<<arguments->m);
//
//
//
//                arguments->sigmaSquare = 1.0 / (2.0 * rate * pow(10, tempDouble/10) );
//            }
//            break;

        case OPT_SNR_CODEWORD_BIT:
            tempDouble = positiveFloatParser("--SNR-codeword", arg, state, false);
            arguments->sigmaSquare = 1.0 / (2.0 * pow(10, tempDouble/10) );
            break;



        case 'g':
            if ( arg == NULL ) // default
            {
                arguments->genieInterventionsAllowed = arguments->k;
            }
            else
            {
                arguments->genieInterventionsAllowed = positiveIntParser("-g (--genie-interventions)", arg, state);
            }
            break;

        case 's':
            arguments->seed = unsignedLongParser("-s (--seed)", arg, state);
            break;

        case 't':
            arguments->trials = positiveIntParser("-t (--trials)", arg, state);
            break;

        case 'o':
            arguments->outfile = arg;
            break;

        case 'd':
        	break;



        case OPT_PRINT_COMMAND_LINE:
            arguments->printCommandLine = yesNoParser("--print-command-line", arg, state);
            break;



        case ARGP_KEY_ARG:
            argp_usage(state);
            break;

        case ARGP_KEY_END:
            //TODO: check for more errors




            if (arguments->trials == -1)
            {
                printf("No -t (--trials) option specified.\n");
                argp_usage (state);
            }

            if (arguments->bec == false && arguments->awgn==false)
            {
                arguments->bsc = true;
            }

            {
                int channelOptions = 0;
                channelOptions += arguments->bec == true ? 1 : 0;
                channelOptions += arguments->bsc == true ? 1 : 0;
                channelOptions += arguments->awgn == true ? 1 : 0;

                if ( channelOptions > 1 )
                {
                    printf("Error: must select at most one channel type.\n");
                    argp_usage (state);

                }
            }

            if ( arguments->bec == true || arguments->bsc == true )
            {
                if (arguments->p == -1.0)
                {
                    printf("Error: -p (--prob) option not set.\n");
                    argp_usage (state);

                }

            }

            if ( arguments->awgn == true )
            {
                if (arguments->sigmaSquare == -1.0)
                {
                    printf("Error: must set one of --sigma, --sigma-square, --SNR-information, or --SNR-codeword.\n");
                    argp_usage (state);

                }

            }


            break;

        default:
            return ARGP_ERR_UNKNOWN;
        }

    }
    return 0;
}


void parse(Arguments &arguments, int argc, char **argv)
{
    /* Start of various structures used to parse input paramaters to program */
    const struct argp_option options[] =
    {
        {NULL, 0, NULL, 0, "Channel specification",0},
        {"input-filename-channel",   'i', "INPUT", 0, "The file INPUT is assumed to have |X| rows and |Y| columns.",0},
//      {"bsc",  OPT_BSC, NULL, 0, "Take the underlying channel as the BSC channel (default).",0},
//      {"bec",  OPT_BEC, NULL, 0, "Take the underlying channel as the BEC channel.",0},
//      {"awgn",  OPT_AWGN, NULL, 0, "Take the underlying channel as the binary input AWGN channel.",0},
//      {"prob",  'p', "P", 0, "Take the positive floating point number P to be the error probability of the bsc or bec channel.",0},
//      {"sigma-square", OPT_SIGMA_SQUARE, "SIGMA-SQUARE", 0, "Take the positive floating point number SIGMA_SQUARE to be the variance of the awgn channel.",0},
//      {"sigma", OPT_SIGMA, "SIGMA", 0, "Take the positive floating point number SIGMA to be the standard deviation of the  awgn channel.",0},
//      {"SNR-information", OPT_SNR_INFORMATION_BIT, "SNR", 0, "Take the positive floating point number SNR to be such that variance (sigma^2) of the  awgn channel equals 1 / (2 * R * 10^(SNR/10)), where R is the code rate. Put another way, SNR = 10 log_10 (1/(2*R*sigma^2)).",0},
//      {"EbN0", OPT_SNR_INFORMATION_BIT, "SNR", OPTION_ALIAS, "", 0},
//      {"SNR-codeword", OPT_SNR_CODEWORD_BIT, "SNR", 0, "Take the positive floating point number SNR to be such that the variance (sigma^2) of the  awgn channel equals 1 / (2 * 10^(SNR/10)).  Put another way, SNR = 10 log_10 (1/(2*sigma^2)).",0},
//      {"EsN0", OPT_SNR_CODEWORD_BIT, "SNR", OPTION_ALIAS, "", 0},
        {NULL, 0, NULL, 0, "Input distribution specification",0},
        {"input-filename-input-distribution",   'd', "INPUT", 0, "The file INPUT is assumed to have the 2^m Inputs probabilities. A line starting with # or * is treated as a comment.",0},

		{NULL, 0, NULL, 0, "Degrade/Upgrade Spec",0},
        {"deg-up",   'u', "INPUT", 0, "Specify Degrading merge [deg-up deg] or Upgrading merge [deg-up up] ",0},
		{"mu",   'm', "M", 0, "This option specifies the mu value for the Degrading/Upgrading merge algorithm",0},
//        {"m",   'm', "M", 0, "This option specifies the m value - number of iterations.",0},
		{"bin-ind",   'b', "B", 0, "This option specifies the mu value for the Degrading/Upgrading merge algorithm",0},

        {NULL, 0, NULL, 0, "Output",0},
        {"output-filename",  'o', "OUTFILE", 0, "Output to OUTFILE instead of to standard output.",0},
        {NULL, 0, NULL, 0, NULL ,0},
        {NULL, 0, NULL, 0, NULL ,0}
    };

    const char args_doc[] = "";

    const char doc[] = "polar list  -- A program to decode polar codes.\vExample usage: polarList --input-filename=channels.txt --trials=1000 -k 350 -l 4 -p 0.11 --crc10\nor: polarList -r 3 -m 8 --trials=1000 -l 4 --awgn --SNR-information=3 --print-misdecoded-histogram=no --print-first-misdecoded-histogram=no";

    const struct argp argp = {options, parse_opt_filename, args_doc, doc, NULL, NULL, NULL};
    /* End of various structures used to parse input paramaters to program */

    /* Set argument defaults */
    arguments.pass = 0;
    arguments.verbosity = -1;
    arguments.p = -1.0;
    arguments.sigmaSquare = -1.0;
    arguments.outfile = NULL;
    arguments.infile = NULL;
    arguments.awgn = false;
    arguments.trials = -1;
    arguments.seed = 0;
    arguments.k = -1;
    arguments.r = -1;
    arguments.genieInterventionsAllowed = 0;
    arguments.binIndLength = 0;

    arguments.printCommandLine = true;


    arguments.degUp = "deg";


    /* Parse command line parameters */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
}

