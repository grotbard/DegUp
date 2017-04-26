/* parser.h
 *
 * Copyright (C) 2012 Ido Tal <idotal@ieee.org> and Alexander Vardy <avardy@ucsd.edu>
 *
 * This file is part of polarList.
 *
 */

#ifndef PARSER
#define PARSER

#include <string>


struct Arguments;

void parse(Arguments &arguments, int argc, char **argv);

struct Arguments
{
    int pass;
    int verbosity;
    float p;
    float sigmaSquare;
    bool bsc;
    bool bec;
    bool awgn;
    char *outfile;            /* Argument for -o */
    char *infile;            /* Argument for -i */
    int trials;
    unsigned long seed;
    int k;
//    int m;
    double mu;
    int binInd; //0b101 for example
    int binIndLength;
    int r;
    int genieInterventionsAllowed;

    std::string degUp; //degrade or upgrade

    std::string commandLine;
    bool printCommandLine;

};

int calcHammingWeight(unsigned int i);

#endif // PARSER
