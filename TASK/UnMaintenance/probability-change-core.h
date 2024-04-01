#pragma once
#include "uncertain-core.h"

class Probability_Core : public Uncertain_Core
{
public:
	//��Ե��������
	void increase_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	vector<vector<double> > increase_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);


	//��Ե���ʼ���
	void decrease_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	vector<vector<double> > decrease_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);



private:

};
