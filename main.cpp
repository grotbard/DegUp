#include <iostream>

#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <iostream>

#include "parser.h"
#include "Channel.h"

using namespace std;

// The Main function starts here.
int main(int argc, char **argv)
{
    int i;






    struct Arguments arguments;

    for ( i = 0; i < argc; i++ )
    {
        arguments.commandLine.append(argv[i]);
        if ( i != argc - 1  )
        {
            arguments.commandLine.append(" ");
        }
    }


    parse(arguments, argc, argv);
    try
    {
    	Channel channel = Channel(arguments.infile);
    	if (arguments.degUp == "deg") {
    		if (arguments.binIndLength == 0)
    			channel.degrade(arguments.mu);
    		else
    			channel.degradingMerge(arguments.mu, arguments.binInd, arguments.binIndLength);
    		channel.printChannel();
    		cout << "P_err of W_plus after degrading merge : " << channel.get_p_err() << endl;
    	    cout << "I(X;Y) of W_plus after degrading merge : " << channel.get_mutual_information() << endl;
    	}
    	else {
    		if (arguments.binIndLength == 0)
    			channel.upgrade(arguments.mu);
    		else
    			channel.upgradingMerge(arguments.mu, arguments.binInd, arguments.binIndLength);
    		channel.printChannel();
    		cout << "P_err of W_plus after upgrading merge : " << channel.get_p_err() << endl;
    		cout << "I(X;Y) of W_plus after upgrading merge : " << channel.get_mutual_information() << endl;

    	}



    }
    catch (const string& str_err)
    {
    	cout << str_err << endl;

    }





    cout << arguments.outfile << endl;
    return 0;
}
