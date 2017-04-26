/*
 * Channel.cpp
 *
 *  Created on: Feb 26, 2017
 *      Author: gal
 */

#include "Channel.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

using namespace std;


Channel::Channel(const string& fullFilePathChannel, const string& fullFilePathInputDist) {
	x_abs = -1;
	y_abs = -1;
	this->p = p;
	string line;
	ifstream myfile(fullFilePathChannel.c_str());
	if (myfile.is_open()){
		int line_n = 0;
		while ( getline (myfile,line) )
	    {
			line_n++;
			stringstream ss(line);
			string item;
			vector<string> tokens;
			while (getline(ss, item, ',')) {
				tokens.push_back(item);
			}

			if (y_abs > -1) {
				// we already have read one line, so let's check that this line has also y_abs entries
				if (y_abs != tokens.size()){
					ostringstream ostr_err;
					ostr_err <<"Error: input file: " << fullFilePathChannel << " line #" << line_n << " should have " << y_abs << "entries";
					throw ostr_err.str();
				}
			}
			else {
				// setting number of columns for the first time
				y_abs = tokens.size();
			}
			vector<myfloat> vec_line;
			for (int i = 0; i < tokens.size(); i++)
			{
				vec_line.push_back(atof(tokens[i].c_str()));
			}
			probMat.push_back(vec_line);
	    }
	    myfile.close();
	    x_abs = line_n;
	    this->cleanZeroColumns();

	    if (fullFilePathInputDist == "")
	    {
	    	//no input distribution - assuming uniform
	    	for (int i = 0; i < x_abs; i++)
	    		inputDistribution.push_back(1/(double)x_abs);
	    }
	    else
	    {
	    	//read input (x) distribution from file - TODO Add this
	    }

	  }

	  else cout << "Unable to open file";

}

int Channel::getRowCount() {
	return x_abs;
}


int Channel::getColCount() {
	return y_abs;
}


void Channel::degrade(double mu) {
	this->cleanZeroColumns();
	myfloat beta = 1/exp(1)/log(2);
	int mu_hat = ceil(beta*mu);

	vector< vector<myfloat> > Q;
	//resize Q [x_abs, (2*mu_hat)^x_abs]
	int y_abs_Q = pow(2*mu_hat, x_abs);
	Q.resize(x_abs);
	for (int i = 0; i < x_abs; i++)
		Q[i].resize(y_abs_Q);
	for (int j = 0; j < y_abs; j++) {
		int b_x = 0;
		for (int i=0; i < x_abs; i++) {
			myfloat P_y = 0;
			for (int k=0; k < x_abs; k++)
				P_y += probMat[k][j]*inputDistribution[k];	//reminder: P(y) = sum_x(W(y|x)*P(x))

			b_x = b_x + (pow(2*mu_hat, i))*(get_b_x(probMat[i][j]*inputDistribution[i]/P_y, mu)-1);
		}
		for (int i=0; i < x_abs; i++)
			Q[i][b_x] += probMat[i][j];
	}
	//updating probMat and other internal vars:
	probMat = Q;
	y_abs = y_abs_Q;
	this->cleanZeroColumns();
}


void Channel::upgrade(double mu) {
	this->cleanZeroColumns();
	vector<myfloat> binningVec = get_b_x_vec_upgrade(mu);

	vector<int> yBinVec; //this vector will contain the binning index corresponding to the index number of y symbol.
	//resize Q [x_abs, (2*mu_hat)^x_abs]
	int M = binningVec.size();
	yBinVec.resize(y_abs);

	for (int j = 0; j < y_abs; j++) {
		int b_x = 0;
		for (int i=0; i < x_abs; i++) {
			myfloat P_y = 0;
			for (int k=0; k < x_abs; k++)
				P_y += probMat[k][j]*inputDistribution[k];	//reminder: P(y) = sum_x(W(y|x)*P(x))

			//finding binning index (l) per x:
			myfloat phiXy = probMat[i][j]*inputDistribution[i]/P_y;
			int binningInd = 0;
			for ( ; binningInd < binningVec.size(); binningInd++) {
				if (binningVec[binningInd] > phiXy) {
					binningInd--;
					break;
				}
			}
			b_x = b_x + (pow(M, i))*binningInd; //TODO: change this mapping to vector mapping like described in the article (you can use vector and the == operator to in order to check if two vectors are equal
		}
		yBinVec[j] = b_x; //the index of the specific y symbol
		}
	//Transform yBinVec to its "inverse" i.e create a vector that contains for each unique key - a vector that points to the y symbol indices
	//Note: the first entry in each vector contains the key value
	vector< vector <int> > yBinVecInv;
	//yBinVecInv.resize(y_abs); //Maximum number of keys is actually y_abs
	//Go over yBinVec and copy to the right place in yBinVecInv
	for (int i=0; i < y_abs; i++) {
		int key = yBinVec[i];
		for (int j=0; j < y_abs; j++) {
			if (j > (int)yBinVecInv.size()-1) {//we reached a new key
				yBinVecInv.push_back(vector<int>());
				yBinVecInv[j].push_back(key); // first element is the key
				yBinVecInv[j].push_back(i); //add "symbol y" to the "linked list"
				break; //for (j) loop
			}
			else {
				if (yBinVecInv[j][0] == key) {
					yBinVecInv[j].push_back(i); //add "symbol y" to the "linked list"
					break; //for (j) loop
				}
			}
		}
	}
	//Initialize epsilon_x
	vector<myfloat> epsilon_x(x_abs); //0,0,...,0 x_abs times
	vector< vector<myfloat> > Q;
	//resize Q [x_abs, (2*mu_hat)^x_abs]
	int y_abs_Q = yBinVecInv.size();
	Q.resize(x_abs);
	for (int i = 0; i < x_abs; i++)
		Q[i].resize(y_abs_Q);
	//find x_star for each "linked list" and store it in x_star_vec
	vector<int> x_star_vec(yBinVecInv.size());
	for (int i = 0; i < yBinVecInv.size(); i++) {
		myfloat phiXyMax = 0;
		for (int j = 0; j < x_abs; j++) {
			for (int k = 1; k < yBinVecInv[i].size(); k++) { //we start with 1 since the first entry is referring to the key value and it is not relevant

				myfloat P_y = 0;
				for (int l=0; l < x_abs; l++)
					P_y += probMat[l][yBinVecInv[i][k]]*inputDistribution[l];	//reminder: P(y) = sum_x(W(y|x)*P(x))
				myfloat phiXy = probMat[j][yBinVecInv[i][k]]*inputDistribution[j]/P_y;
				if (phiXy > phiXyMax) {
					phiXyMax = phiXy;
					x_star_vec[i] = j;
				}
			}
		}
	}
	for (int i = 0; i < yBinVecInv.size(); i++) {
		//calculate psy_x_z (for specific z i.e for specific yBinVecInv entry
		vector<myfloat> psy_x_z_vec(x_abs);
		myfloat sumTmp = 0;
		for (int j = 0; j < x_abs; j++) {
			myfloat phiMin = 1;
			if (j != x_star_vec[i]) {
				for (int k = 1; k < yBinVecInv[i].size(); k++) { //we start with 1 since the first entry is referring to the key value and it is not relevant
					myfloat P_y = 0;
					for (int l=0; l < x_abs; l++)
						P_y += probMat[l][yBinVecInv[i][k]]*inputDistribution[l];	//reminder: P(y) = sum_x(W(y|x)*P(x))
					myfloat phiXy = probMat[j][yBinVecInv[i][k]]*inputDistribution[j]/P_y;
					if (phiXy < phiMin)
						phiMin = phiXy;
				}
				sumTmp += phiMin;
				psy_x_z_vec[j] = phiMin;
			}
		}
		//updating the x*-th entry:
		psy_x_z_vec[x_star_vec[i]] = 1 - sumTmp;

		for (int j = 0; j < x_abs; j++) {
			for (int k = 1; k < yBinVecInv[i].size(); k++) { //we start with 1 since the first entry is referring to the key value and it is not relevant

				//calculate phiXbar_y
				myfloat P_y = 0;
				for (int l=0; l < x_abs; l++)
					P_y += probMat[l][yBinVecInv[i][k]]*inputDistribution[l];	//reminder: P(y) = sum_x(W(y|x)*P(x))
				myfloat phiXbar_y = probMat[x_star_vec[i]][yBinVecInv[i][k]]*inputDistribution[x_star_vec[i]]/P_y;
				myfloat phiXy = probMat[j][yBinVecInv[i][k]]*inputDistribution[j]/P_y;
				myfloat alpha_x = 0;
				if (phiXy > 0)
					alpha_x = psy_x_z_vec[j]/phiXy*phiXbar_y/psy_x_z_vec[x_star_vec[i]];
				else
					alpha_x = 1;
				epsilon_x[j] += (1-alpha_x)*probMat[j][yBinVecInv[i][k]];
				Q[j][i] += alpha_x*probMat[j][yBinVecInv[i][k]];
			}
		}
	}

	//Produce boost symbols and probabilities
	for (int i = 0; i < x_abs; i++) {
		if (epsilon_x[i]>0)
			for (int j = 0; j < x_abs; j++) {
				if (i == j)
					Q[j].push_back(epsilon_x[i]);
				else
					Q[j].push_back(0);
			}
	}
	//updating probMat and other internal vars:
	probMat = Q;
	y_abs = Q[0].size();
	this->cleanZeroColumns();


}




void Channel::plus() {
	vector< vector<myfloat> > W_plus;
	int y_abs_W_plus = x_abs*pow(y_abs, 2);
	W_plus.resize(x_abs);
	for (int i = 0; i < x_abs; i++)
		W_plus[i].resize(y_abs_W_plus);

	for (int u1_ind = 0; u1_ind < x_abs; u1_ind++)
		for (int u0_ind = 0; u0_ind < x_abs; u0_ind++)
			for (int y1_ind = 0; y1_ind < y_abs; y1_ind++)
				for (int y0_ind = 0; y0_ind < y_abs; y0_ind++)
					W_plus[u1_ind][u0_ind+y1_ind*x_abs+y0_ind*y_abs*x_abs] = 1/((double)x_abs)*probMat[(u0_ind+u1_ind)%x_abs][y0_ind]*probMat[u1_ind][y1_ind];
	probMat = W_plus;
	y_abs = y_abs_W_plus;
	this->cleanZeroColumns();

}

void Channel::minus() {
	vector< vector<myfloat> > W_minus;
	int y_abs_W_minus = pow(y_abs, 2);
	W_minus.resize(x_abs);
	for (int i = 0; i < x_abs; i++)
		W_minus[i].resize(y_abs_W_minus);

	for (int u1_ind = 0; u1_ind < x_abs; u1_ind++)
		for (int u0_ind = 0; u0_ind < x_abs; u0_ind++)
			for (int y1_ind = 0; y1_ind < y_abs; y1_ind++)
				for (int y0_ind = 0; y0_ind < y_abs; y0_ind++) {
					myfloat W_temp = W_minus[u0_ind][y1_ind+y0_ind*y_abs];
					W_minus[u0_ind][y1_ind+y0_ind*y_abs] = W_temp + 1/((double)x_abs)*probMat[(u0_ind+u1_ind)%x_abs][y0_ind]*probMat[u1_ind][y1_ind];
				}
	probMat = W_minus;
	y_abs = y_abs_W_minus;
	this->cleanZeroColumns();

}

void Channel::degradingMerge(double mu, int bin_ind, int bin_ind_length) {
	if (bin_ind_length == 0) { //if bin_ind_length not given, initialize it to length of bin_ind where msb is the last "1"
		int bin_ind_tmp = bin_ind;
		while (bin_ind_tmp >0) {
			bin_ind_length++;
			bin_ind_tmp = bin_ind_tmp >> 1;
		}
	}

	this->degrade(mu);
	for (int i = 0; i < bin_ind_length; i++) {
		if (bin_ind&1)
			this->plus();
		else
			this->minus();
		this->degrade(mu);
		bin_ind = bin_ind>>1;
	}
}

void Channel::upgradingMerge(double mu, int bin_ind, int bin_ind_length) {
	if (bin_ind_length == 0) { //if bin_ind_length not given, initialize it to length of bin_ind where msb is the last "1"
		int bin_ind_tmp = bin_ind;
		while (bin_ind_tmp >0) {
			bin_ind_length++;
			bin_ind_tmp = bin_ind_tmp >> 1;
		}
	}

	this->upgrade(mu);
	for (int i = 0; i < bin_ind_length; i++) {
		if (bin_ind&1)
			this->plus();
		else
			this->minus();
		this->upgrade(mu);
		bin_ind = bin_ind>>1;
	}
}


myfloat Channel::get_p_err() {

	myfloat sumTmp = 0;
	for (int i = 0; i < y_abs; i++) {
		int maxIndTmp = 0;
		for (int j = 1; j < x_abs; j++) //we start with 1 because 0 is maxIndTmp
			if (probMat[j][i] > probMat[maxIndTmp][i])
				maxIndTmp = j;
		//maxIndicesVec[i] = maxIndTmp;
		for (int j = 0; j < x_abs; j++)
			if (j != maxIndTmp)
				sumTmp += inputDistribution[j]*probMat[j][i];
	}
	return sumTmp;
}


myfloat Channel::get_mutual_information() {
	myfloat Ixy = 0;
	for (int i = 0; i < x_abs; i++)
		for (int j=0; j < y_abs; j++)
			if (probMat[i][j] > 0) {
				myfloat P_y = 0;
				for (int k=0; k < x_abs; k++)
					P_y += probMat[k][j]*inputDistribution[k];	//reminder: P(y) = sum_x(W(y|x)*P(x))
				Ixy += inputDistribution[i]*probMat[i][j]*log2(probMat[i][j]/P_y);
			}
	return Ixy;


}


void Channel::cleanZeroColumns() {
	//cleaning columns in which all entries are zeroes (p(y|x)=0 for all y)
	for (int i = 0; i < y_abs; i++)
		for (int j = 0; j < x_abs; j++)
			if (probMat[j][i])
				break;
			else if (j == x_abs-1){
				//all entries in that column are zeroes
				for (int k = 0; k < x_abs; k++)
					probMat[k].erase(probMat[k].begin()+i);
				y_abs--;
				i--;
			}


}


void Channel::printChannel() {
	cout << "      Y ..." << endl;
	cout << "X" << endl;
	for (int i=0; i<x_abs; i++) {
		cout << ".";
		for (int j=0; j<y_abs; j++)
			cout << "    " << setw( -20 ) << probMat[i][j];
		cout << endl;
	}
}


int Channel::get_b_x(myfloat x, myfloat mu) {
	if ((x > 1) || (x < 0))
		throw string("Error in Channel::get_b_x: x should be in [0,1]. This may be caused by overflow in myfloat type");
	myfloat alpha = 1/exp(1);
	myfloat beta = 1/exp(1)/log(2);
	int mu_hat = ceil(beta*mu);
	myfloat etta_x = 0;
	if (x > 0)
		etta_x = -x*log2(x);

	if (x < alpha)
		return (floor(mu*etta_x)+1);
	else if (x==alpha)
		return (mu_hat);
	else {
		int j = floor(mu*etta_x)+1;
		return (2*mu_hat+1-j);
	}
}

vector<myfloat> Channel::get_b_x_vec_upgrade(myfloat mu) { //return an 1xM vector in which the first entry is 0 and the last is b_M
	vector<myfloat> b_x;
	b_x.push_back(0);
	myfloat lastBiElement = 0;
	myfloat alphaSquare = 1/exp(2);
	myfloat delta = 1/mu*0.0000000001; //arbitrary
	myfloat etta_bi = 0;
	while (lastBiElement < 1) {
		if (lastBiElement < alphaSquare) { //calculating the section where the slope is greater than 1
			myfloat binSearchLastStep = 1/mu/2;
			myfloat tmpBinSearchLastRes = lastBiElement+binSearchLastStep;
			while (fabs(-tmpBinSearchLastRes*log(tmpBinSearchLastRes) - etta_bi -1/mu) > delta ) { //since we want to be as close as delta to solve the equation
				binSearchLastStep /= 2;
				if (-tmpBinSearchLastRes*log(tmpBinSearchLastRes) - etta_bi > 1/mu)
					tmpBinSearchLastRes -= binSearchLastStep;
				else
					tmpBinSearchLastRes += binSearchLastStep;
			}
			lastBiElement = tmpBinSearchLastRes;
			etta_bi = -lastBiElement*log(lastBiElement);
			b_x.push_back(lastBiElement);

		}
		else { //slope is strictly less than 1 and strictly greater than -1
			//check the first element after alphaSquare>1 (it has been pushed to the vector in the last iteration) since it might violate the condition of p<=b_i+1/mu
			if (lastBiElement - b_x[b_x.size()-2] > 1/mu)
				b_x[b_x.size()-1] = b_x[b_x.size()-2] + 1/mu; //change the last element to comply the horizontal condition.

			lastBiElement += 1/mu;
			if (lastBiElement < 1)
				b_x.push_back(lastBiElement);
		}
	}
	return b_x;
}


Channel::~Channel() {
	// TODO Auto-generated destructor stub
}

