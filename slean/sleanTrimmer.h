/*
Slean Edge Trimmer for Cuckoo Cycle 
Made by Wilke Trei, March 2019
*/

#include <iostream>
#include <cstring>
#include <chrono>
#include <vector>
#include <map>
#include <omp.h>
#include "./crypto/blake2.h"

using namespace std;

 
const uint32_t sleanSegs[] = {16,10,5,3,2,2};
//const uint32_t sleanSegs[] = {128,64,32,16,16,8};

struct sleanEdge {
	uint32_t endpoint;
	uint32_t nonce;
};


struct fullEdge {
	uint32_t endpointU;
	uint32_t endpointV;
	uint32_t nonce;	
	bool root;
};


class sleanTrimmer {
	private:
		uint32_t nSt;
		uint64_t startEdges;

		uint32_t bucketBits;
		uint32_t buckets;

		uint64_t endpointMask;
		uint32_t shiftBits;
		uint32_t addrMask;

		uint32_t bucketSizeA;
		uint32_t bucketSizeASeg;
		uint32_t numThreads;

		vector<uint32_t> nodeDegreeMap;
		vector<uint32_t> edgeLivenessMap;

		vector<sleanEdge> bucketSpaceA;
		vector<uint32_t> countersA;

		/* seeding */
		uint64_t v[4];
		uint64_t dipnode(uint64_t, uint64_t); 			// Parameters: Nonce, uOrv

		/* slean */
		void clearCounters();					
		void seed(uint32_t, uint32_t, uint32_t);		// Parameters: start, end, uOrv
		void mark();
		void trim();			
		void read();

		/* cycle finding */
		uint32_t cycleShiftBits;		

		void seedCycleFinder();	
		uint32_t getElem(uint32_t);
		void setElem(uint32_t, uint32_t);
		uint32_t getPath(uint32_t, vector<uint32_t> &);	
		

	public:
		void setup(uint32_t, vector<uint8_t>, uint32_t, uint32_t, bool); 	//Parameters: log_2 n, header, nonce, threads, appendNonce
		void work(uint32_t);							//Parameters: rounds
	
};	
