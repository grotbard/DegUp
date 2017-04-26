/*
 * Channel.h
 *
 *  Created on: Feb 26, 2017
 *      Author: gal
 */

#ifndef CHANNEL_H_
#define CHANNEL_H_

#include "myfloat.h"
#include <vector>
#include <string>
using namespace std;

class Channel {
public:
	//Constructors
	Channel(const string& fullFilePathChannel, const string& fullFilePathInputDist ="");

	//Public Methods
	int getRowCount(); //x_abs
	int getColCount(); //y_abs
	void degrade(double mu);
	void upgrade(double mu);
	void plus(); //this <= W+
	void minus(); //this <= W-
	void degradingMerge(double mu, int bin_ind, int bin_ind_length = 0);
	void upgradingMerge(double mu, int bin_ind, int bin_ind_length = 0);
	myfloat get_p_err();
	myfloat get_mutual_information();
	void printChannel();






	//Public Variables
	vector< vector<myfloat> > probMat;
	vector<myfloat> inputDistribution;
	int x_abs;
	int y_abs;
	int p; //alphabet size


	//Static Methods:
	static int get_b_x(myfloat x, myfloat mu);
	static vector<myfloat> get_b_x_vec_upgrade(myfloat mu);

	// Destructor
	~Channel();

private:
	//Private Methods
	void cleanZeroColumns();
};

#endif /* CHANNEL_H_ */
