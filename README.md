# Cuckoo Cycle Slean Edge Trimmer

A submission to the CPU Speedup Bounty challenge by John Trom,
https://github.com/tromp/cuckoo

The miner is build to use less then one byte per edge but has better 
parallelization scaling then the original lean miner. 


## Bulding

The solver requires cMake 3.9 and newer to handle the parallelization 
with OpenMP correctly.
Then you can build with
`cmake -DCMAKE_BUILD_TYPE=Release . && make `

## Usage

The solver uses the same parameters as the original lean miner.

New parameter:
Use -e  to set the number of edge bits. 

## How it works

Watch my talk at the 2019 Grin Amsterdam Meetup to see how it works.
https://youtu.be/IFtbVUfMYpc

Modifications:
On CPU it is more efficient to store a pair (nonce, edge endpoint), so
when marking and trimming we do not need to recompute the endpoints. 
Also I write directly back to the edge aliveness map instead of bucketing
by the nonces after the trimming rounds.

