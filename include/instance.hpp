#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#define INFINITE 0

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <memory>

class Instance {
private:
	std::string m_name; 
	size_t m_dimension;
	size_t m_num_parameters;

	std::shared_ptr<double*[]> distMatrix;
	std::shared_ptr<double[]> xCoord, yCoord;

	//Computing Distances
	static double CalcDistEuc ( double *, double *, int , int );
	static double CalcDistAtt ( double *, double *, int , int );
	static double CalcDistGeo ( double *, double *, int , int );
	static void CalcLatLong ( double *, double *, int , double *, double* );

	bool explicitCoord;
public:
	Instance( int, char * );

	void read();
	void printMatrixDist();

	inline size_t get_dimension() const { return m_dimension; } 

	inline double get_distance(size_t i, size_t j) const { return distMatrix[i-1][j-1]; }
	
	inline double **getMatrixCost() { return distMatrix.get(); }
	inline double getXCoord(int i){return xCoord[i-1];}
	inline double getYCoord(int i){return yCoord[i-1];}
	inline bool getExplicitCoord(){return explicitCoord; };

	std::string get_name();
};

#endif
