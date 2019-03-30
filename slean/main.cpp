/*
Slean Edge Trimmer for Cuckoo Cycle 
Made by Wilke Trei, March 2019
*/

#include <iostream>
#include <unistd.h>

#include "sleanTrimmer.h"

using namespace std;

int main(int argc, char **argv) {
	uint32_t numThreads = 1;
	uint32_t nonce = 0;
	uint32_t range = 1;
	uint32_t ntrims = 0;
	uint32_t edgeBits = 29;

	vector<uint8_t> header;
	header.assign(80,0);

	uint32_t len;
	int32_t c;

	while ((c = getopt (argc, argv, "e:h:m:n:r:t:x:")) != -1) {
		switch (c) {
		case 'e':
			edgeBits = atoi(optarg);
			break;
	        case 'h':
			len = strlen(optarg);
			memcpy(header.data(), optarg, len);
			break;
		case 'x':
			len = strlen(optarg)/2;
			for (uint32_t i=0; i<len; i++) {
				sscanf(optarg+2*i, "%2hhx", header.data()+i);
			}
			break; 
		case 'n':
			nonce = atoi(optarg);
			break;
		case 'r':
			range = atoi(optarg);
			break;
		case 'm':
			ntrims = atoi(optarg);
			break;
		case 't':
			numThreads = atoi(optarg);
			break;
		}
	}

	if (ntrims == 0) {
		ntrims = (edgeBits < 31) ? 20 : 20+numThreads*2;
	}


	for (uint32_t i=nonce; i < nonce+range; i++) {

		cout << endl;
		cout << "---------------------------------------" << endl;
		cout << "     - Cuckoo Slean Edge Trimmer - " << endl;
		cout << "---------------------------------------" << endl << endl;

		cout << "Looking for 42-cycles in cuckoo-" << edgeBits+1 << " with nonce: " << i << endl;
		sleanTrimmer myTrimmer;
		
		myTrimmer.setup(edgeBits, header, i, numThreads, false);
		myTrimmer.work(ntrims);
	}
	
}
