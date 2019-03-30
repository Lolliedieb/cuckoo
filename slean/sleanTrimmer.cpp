/*
Slean Edge Trimmer for Cuckoo Cycle 
Made by Wilke Trei, March 2019
*/

#include "sleanTrimmer.h"

void sleanTrimmer::setup(uint32_t n, vector<uint8_t> header, uint32_t nonce, uint32_t threads, bool attach) {

	/* ---------------------
		Determine number of buckets and assign memory 
	--------------------- */

	omp_set_num_threads(threads);
	numThreads = omp_get_max_threads();

	nSt = n;
	bucketBits = (n-9) / 2;
	buckets = (1 << bucketBits);

	cout << "Number of threads: " << numThreads << endl;

	shiftBits = n-bucketBits;
	startEdges = (1ULL << n);
	endpointMask = startEdges - 1;
	addrMask = (1 << (shiftBits-5)) - 1;

	cout << "Addr Mask: " << addrMask << endl;
	

	bucketSizeA = (4352 << (n-16)) / buckets; 

	cout << "Number of buckets: " << buckets << endl;
	cout << "Elements per bucket: " << bucketSizeA << endl;

	bucketSizeASeg = bucketSizeA/numThreads;

	cout << "Elements per bucket segment: " << bucketSizeASeg << endl;

	countersA.assign(buckets*numThreads, 0);
	bucketSpaceA.resize(bucketSizeA * buckets);

	uint64_t numBytes = bucketSizeA*buckets*8;
	cout << "Using " << numBytes << " bytes (" << numBytes / (1024*1024) << " MBytes) of bucket space" << endl;

	nodeDegreeMap.assign( (2 << (n-5)), 0);
	edgeLivenessMap.assign( (1 << (n-5)), 0xFFFFFFFF);

	numBytes = nodeDegreeMap.size()*4;
	cout << "Using " << numBytes << " bytes (" << numBytes / (1024*1024) << " MBytes) for the node degree map" << endl;
	numBytes = edgeLivenessMap.size()*4;
	cout << "Using " << numBytes << " bytes (" << numBytes / (1024*1024) << " MBytes) for the edge liveness map" << endl;




	/* ---------------------
		Prepare the header and pre-process the seeding 
	--------------------- */
	vector<uint8_t> headerCpy;
	headerCpy.insert(headerCpy.end(), header.begin(), header.end());

	//Normally we should add the nonce to the end, but also allow overwriting for compatiblity 
	if (attach) {
		for (int i=0; i<3; i++) headerCpy.push_back( (uint8_t) 0);
	}

	if (headerCpy.size() < 4) {
		cout << "Initialization failed: too small header \n";
		exit(1);
	}

	// Appending the nonce
	*((uint32_t * ) (&headerCpy.data()[headerCpy.size() - 4])) = htole32(nonce);

	// Blake2B to determine our seeds
	blake2b_state target_state;
	blake2b_init(&target_state, 32);
	blake2b_update(&target_state, headerCpy.data(), headerCpy.size());
	blake2b_final(&target_state, (uint8_t *) &v[0], 32);	
}

#define rotl(x,b) ((x << b) | ( x >> (64 - b)))

#define SIPROUND \
    v0 += v1; v2 += v3; v1 = rotl(v1,13); 	\
    v3 = rotl(v3,16); v1 ^= v0; v3 ^= v2;	\
    v0 = rotl(v0,32); v2 += v1; v0 += v3;	\
    v1 = rotl(v1,17);   v3 = rotl(v3,21);	\
    v1 ^= v2; v3 ^= v0; v2 = rotl(v2,32);


uint64_t sleanTrimmer::dipnode(uint64_t nce, uint64_t uorv) {
	uint64_t nonce = 2 * nce + uorv;
	uint64_t v0 = v[0], v1 = v[1], v2 = v[2], v3 = v[3] ^ nonce;
	SIPROUND; SIPROUND;
	v0 ^= nonce;
	v2 ^= 0xff;
	SIPROUND; SIPROUND; SIPROUND; SIPROUND;
	return (v0 ^ v1 ^ v2  ^ v3) & endpointMask;
}
  

void sleanTrimmer::seed(uint32_t start, uint32_t end, uint32_t uOrv) {
	#pragma omp parallel for schedule(static)
	for (uint32_t blk=start; blk<end; blk++) {

		uint32_t threadId = omp_get_thread_num();
		uint32_t offset = threadId*buckets*bucketSizeASeg;
		uint32_t * counter = &countersA.data()[buckets*threadId]; 
	
		uint32_t live = edgeLivenessMap[blk];
		while (live != 0) {
			uint32_t tZero = __builtin_clz(live);
			uint32_t index = 31-tZero;

			sleanEdge myEdge;
			myEdge.nonce = (uint32_t) ((blk << 5) + index);
			myEdge.endpoint = (uint32_t) dipnode((uint64_t) myEdge.nonce, (uint64_t) uOrv );

			uint32_t bucket = myEdge.endpoint >> shiftBits;
			uint32_t pos = counter[bucket]++;
		
			bucketSpaceA[offset + bucket * bucketSizeASeg + pos] = myEdge;

			live = live & (0x7FFFFFFF >> tZero); 
		}
	} 	
}


void sleanTrimmer::mark() {
	#pragma omp parallel for schedule(static)
	for (uint32_t bucket=0; bucket<buckets; bucket++) {
		
		uint32_t nodesPerBucket = addrMask+1;
		uint32_t * bucketNodes = &nodeDegreeMap[bucket*2*nodesPerBucket];

		for (uint32_t seg = 0; seg<numThreads; seg++) {
			uint32_t elements = countersA[seg*buckets + bucket];
			uint32_t offset = seg*buckets*bucketSizeASeg + bucket*bucketSizeASeg;
			
			for (uint32_t elem=0; elem<elements; elem++) {
				sleanEdge myEdge = bucketSpaceA[offset + elem];

				uint32_t targ = (myEdge.endpoint >> 5) & addrMask;			
				uint32_t one = (1 << (myEdge.endpoint & 0x1F)); 

				// Two bit counter for cuckoo
				if ((bucketNodes[targ] & one) != 0) {
					bucketNodes[nodesPerBucket + targ] |= one;
				}
				
				bucketNodes[targ] |= one;
			}
		}
	}
}


void sleanTrimmer::trim() {

	#pragma omp parallel for schedule(static,1)
	for (uint32_t seg = 0; seg<numThreads; seg++) {
		uint32_t nodesPerBucket = addrMask+1;

		for (uint32_t bucket=0; bucket<buckets; bucket++) {

			uint32_t * bucketNodes = &nodeDegreeMap[bucket*2*nodesPerBucket];
			uint32_t offset = seg*buckets*bucketSizeASeg + bucket*bucketSizeASeg;
			uint32_t elements = countersA[seg*buckets + bucket];

			for (uint32_t elem=0; elem<elements; elem++) {

				sleanEdge myEdge = bucketSpaceA[offset + elem];

				//uint32_t addr = bucketSpaceA[offset + 2*elem];
				//uint32_t nonce = bucketSpaceA[offset + 2*elem + 1];

				uint32_t targ = (myEdge.endpoint >> 5) & addrMask;			
				uint32_t one = (1 << (myEdge.endpoint & 0x1F)); 

				// Two bit counter for cuckoo
				if ((bucketNodes[nodesPerBucket + targ] & one) == 0) {
					uint32_t targ2 = (myEdge.nonce >> 5);			
					uint32_t one2 = ~(1 << (myEdge.nonce & 0x1F));

					edgeLivenessMap[targ2] &= one2;
				}
			}	
		}
	}
}


void sleanTrimmer::clearCounters() {
	memset(countersA.data(), 0, buckets*numThreads*4);
}



uint32_t sleanTrimmer::getElem(uint32_t uIndex) {
	uint32_t start = uIndex >> cycleShiftBits;
	for (uint32_t ui = start; ; ui++) {
		sleanEdge edge = bucketSpaceA[ui];
		
		if ((edge.endpoint == 0) && (edge.nonce == 0)) {
			return 0;
		} 

		if (edge.endpoint == uIndex) {
			return edge.nonce;
		}
	}
}


void sleanTrimmer::setElem(uint32_t uIndex, uint32_t dest) {
	uint32_t start = uIndex >> cycleShiftBits;
	sleanEdge nEdge;
	nEdge.endpoint = uIndex; nEdge.nonce = dest;

	for (uint32_t ui = start; ; ui++) {
		sleanEdge edge = bucketSpaceA[ui];
		
		if (((edge.endpoint == 0) && (edge.nonce == 0)) || (edge.endpoint == uIndex)) {
			bucketSpaceA[ui] = nEdge;
			break;
		} 
	}
}


uint32_t sleanTrimmer::getPath(uint32_t uIndex, vector<uint32_t> &path) {
	uint32_t length;
	uint32_t lastIndex = uIndex;
	for (length = 0; (lastIndex != 0) && (length < 8192); lastIndex=getElem(lastIndex)) {
		path[length++] = lastIndex;
	}

	return length-1;
}


void sleanTrimmer::seedCycleFinder() {
	vector<uint32_t> pathU,pathV;
	pathU.assign(8192,0);
	pathV.assign(8192,0);
	
	for (uint32_t blk=0; blk < (startEdges >> 5); blk++) {
	
		uint32_t live = edgeLivenessMap[blk];
		while (live != 0) {
			uint32_t tZero = __builtin_clz(live);
			uint32_t index = 31-tZero;

			uint32_t nonce = (uint32_t) ((blk << 5) + index);
			uint32_t endpointU = (uint32_t) dipnode((uint64_t) nonce, 0);
			uint32_t endpointV = (uint32_t) dipnode((uint64_t) nonce, 1) | (uint32_t) startEdges;
			
			if (endpointU) {
				int32_t lenU = getPath(endpointU, pathU);
				int32_t lenV = getPath(endpointV, pathV);

				if (pathU[lenU] == pathV[lenV]) {
					while (pathU[lenU] == pathV[lenV]) {
						lenU--; lenV--;
						if ((lenU < 0) || (lenV < 0)) break;
					} 
					lenU++; lenV++;

					uint32_t cycleLength = lenU + lenV + 1;
					cout << cycleLength << " cycle found" << endl;

					if (cycleLength == 42) {
						cout << "Found solution: ";
					
						for (int32_t st = 0; st < lenU; st++) {
							uint32_t first, second;
							uint32_t src = pathU[st];
							uint32_t dst = getElem(src);
							first = (src & (uint32_t) startEdges) ? dst : src; 
							second = (src & (uint32_t) startEdges) ? src : dst; 
							second ^= (uint32_t) startEdges;
							
							cout << "(" << first << "," << second << ") ";
						}

						for (int32_t st = lenV-1; st >= 0; st--) {
							uint32_t first, second;
							uint32_t src = pathV[st];
							uint32_t dst = getElem(src);
							first = (src & (uint32_t) startEdges) ? dst : src; 
							second = (src & (uint32_t) startEdges) ? src : dst; 
							second ^= (uint32_t) startEdges;
							
							cout << "(" << first << "," << second << ") ";
						}

						cout << "(" << endpointU << "," << (endpointV ^ startEdges)  << ")" << endl;
					}
				} else if (lenU < lenV) {
					while (lenU--) {
						setElem(pathU[lenU+1], pathU[lenU]);
					}
					setElem(endpointU, endpointV);
				} else {
					while (lenV--) {
						setElem(pathV[lenV+1], pathV[lenV]);
					}
					setElem(endpointV, endpointU);
				}
				
			}

			live = live & (0x7FFFFFFF >> tZero); 
		}
	} 	
}



void sleanTrimmer::work(uint32_t rounds) {
	
	uint32_t round = 0;

	cout << endl;
	cout << "---------------------------------------" << endl;
	cout << "   - Starting Slean Edge Trimmer - " << endl;
	cout << "---------------------------------------" << endl << endl;

	auto timeS = std::chrono::high_resolution_clock::now();
	while (round < rounds) {
		uint32_t sleanSegments = round > 5 ? 1 : sleanSegs[round];

		memset(nodeDegreeMap.data(), 0, startEdges >> 2);
		
		// Phase 1: Seed / Mark
		for (uint32_t seg = 0; seg < sleanSegments; seg++) {

			uint64_t start = (seg * startEdges) / (sleanSegments);
			uint64_t end = ((seg+1) * startEdges) / (sleanSegments);
			
			clearCounters();
			seed(start >> 5, end >> 5, round & 0x1);
			mark();
		} 


		// Phase 2: Seed / Trim
		trim();
		for (uint32_t seg = 0; seg < sleanSegments-1; seg++) {
			uint64_t start = (seg * startEdges) / sleanSegments;
			uint64_t end = ((seg+1) * startEdges) / sleanSegments;

			clearCounters();
			seed(start >> 5, end >> 5, round & 0x1);
			trim();
		} 
		
		round++;

		if ((round % 4) == 1) {
			cout << "Slean round " << round << " done" << endl;

			auto timeE = std::chrono::high_resolution_clock::now();
			int64_t timeDiff = std::chrono::duration_cast<std::chrono::milliseconds>(timeE - timeS).count();
			cout << "Duration: " << timeDiff << "ms" << endl;
		}
	}

	auto timeE = std::chrono::high_resolution_clock::now();
	int64_t timeDiff = std::chrono::duration_cast<std::chrono::milliseconds>(timeE - timeS).count();

	uint64_t surviving=0;
	for (uint32_t i=0; i<edgeLivenessMap.size(); i++) {
		surviving += __builtin_popcount(edgeLivenessMap[i]);
	} 	

	cout << endl;
	cout << "---------------------------------------" << endl;
	cout << "- Done Slean Edge Trimmer (" << rounds << " rounds) -   " << endl;
	cout << "  Duration: " << timeDiff << "ms" << endl;
	cout << "  Surviving edges: " << surviving << " (" << ((double) surviving * 100.0) / (double) (32*edgeLivenessMap.size()) << " %)" << endl;	
	cout << "---------------------------------------" << endl << endl;


	uint32_t topZ = __builtin_clz(surviving);
	uint32_t bits = (33-topZ);
	cycleShiftBits = nSt + 1 - bits;

	fullEdge rootNode;
	rootNode.root = true;
	memset(bucketSpaceA.data(), 0, bucketSpaceA.size() * sizeof(uint64_t));

	/* fullEdgesU.assign(cycleBucketSize*cycleBuckets,rootNode);
	fullEdgesV.assign(cycleBucketSize*cycleBuckets,rootNode);  */

	seedCycleFinder();
	//findCycles(total);	
}
