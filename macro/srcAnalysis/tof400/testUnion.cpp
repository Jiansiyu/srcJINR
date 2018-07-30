#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>


using namespace std;

struct TOFHit{
	double t,a,y;
	int st;
	bool hit;
	TOFHit(){
		t = 0;
		a = 0;
		y = 0;
		st = -1;
		hit = false;
	}
};

int main(int argc, char ** argv)
{
	if (argc < 1)
	{
		cerr << "Wrong number of arguments.\n";
		return -1;
	}

	std::vector<int> clusterList[1][48];

	clusterList[0][0] = { 0, 1, 2};
	clusterList[0][1] = { 1, 2};
	clusterList[0][2] = { 2, 3};
	clusterList[0][3] = { 3};
	clusterList[0][4] = { 4, 5};
	clusterList[0][5] = { 5};
	clusterList[0][6] = { 6, 7};
	clusterList[0][7] = { 7, 8};
	clusterList[0][8] = { 8,10};
	clusterList[0][9] = { 9,10};
	clusterList[0][10] = { 10};


	std::vector< std::vector<int> > result;
		// For all strips, take the cluster list
	int pl = 0;
	for( int st = 0 ; st < 48 ; st++){
		if( st > 10 ) break;
		cout << "Working on strip " << st << "\n";
			// For all the groups I already have, find out if any intersection with this cluster list, and if so, add it to group
		bool insert = true;
		cout << "number of current clusters size: " << result.size() << "\n";
		for( int group = 0 ; group < result.size() ; group++){
			std::vector<int> tmp;
			std::vector<int> newGrp;
			std::set_intersection( clusterList[pl][st].begin() , clusterList[pl][st].end(), \
						result.at(group).begin() , result.at(group).end() , back_inserter(tmp) ) ;
			cout << "intersection of how many elements?: " << tmp.size() << "\n";
			if( tmp.size() ){
				// Union two sets into the group and exit
				std::set_union( clusterList[pl][st].begin() , clusterList[pl][st].end(), \
						result.at(group).begin() , result.at(group).end() , back_inserter( newGrp ) );
				insert = false;
				result.at(group).swap(newGrp);
				break;
			}
		}
			// Create new group if no group wiht an intersection
		if( insert ){
			result.push_back(  clusterList[pl][st] );
		}

	}
	
	cout << "number of groups: " << result.size() << "\n";
	for( int group = 0 ; group < result.size() ; group++){
		for( int j = 0 ; j < result.at(group).size() ; j++) cout << result.at(group).at(j) << " ";
		cout << "\n";
	}

	return 0;
}
