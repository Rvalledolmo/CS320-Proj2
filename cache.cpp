#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

class CacheBlock {
public:
    int tag;               
    bool valid;            
    bool dirty;           
    int lru_counter;     
    bool hot_cold_bit;   
    std::vector<char> data;

    CacheBlock(int blockSize)
        : tag(-1), valid(false), dirty(false), lru_counter(0), hot_cold_bit(false), data(blockSize, 0) {}
};

// parse a trace file into a vector of pairs (operation, address)
vector<pair<char, unsigned int> > parseTraceFile(const string& traceFile) {
    vector<pair<char, unsigned int> > traces;
    ifstream infile(traceFile);
    if (!infile) {
        cerr << "Error opening trace file!" << endl;
        return traces;
    }
    char operation;
    unsigned int address;

    while (infile >> operation >> hex >> address) {
        traces.emplace_back(operation, address);
    }

    return traces;
}

bool isSetFull(const vector<CacheBlock> &set) {
    int validCount = 0;
    for (const CacheBlock &block : set) {
        if (block.valid) {
            validCount++;
        }
    }
    return validCount == set.size();
}

int findLRUBlockIndex(const vector<CacheBlock> &set) {
    int lruIndex = 0;
    int maxLRU = set[0].lru_counter;
    for (int i = 1; i < set.size(); i++) {
        if (set[i].lru_counter > maxLRU) {
            lruIndex = i;
            maxLRU = set[i].lru_counter;
        }
    }
    return lruIndex;
}

int findFirstColdBlock(const vector<CacheBlock>& cache) {
    for (int i = 0; i < cache.size(); ++i) {
        if (!cache[i].hot_cold_bit) {
            return i;
        }
    }
    return -1;
}

void DirectMappedCache(const vector<pair<char, unsigned int> > &traces, int cacheSize,  const string& outputFile) {
    int lineSize = 32;
    int numLines = cacheSize / lineSize;
    vector<CacheBlock> cache(numLines, CacheBlock(lineSize));
    ofstream outfile(outputFile, ios::app);

    int hits = 0;
    int misses = 0;

    for (int i = 0; i < traces.size(); ++i) {
        char operation = traces[i].first;
        int address = traces[i].second;
        int blockOffset = address & (lineSize - 1); 
        int index = (address / lineSize) & (numLines - 1);
        int tag = address >> (5 + 5);
        CacheBlock &block = cache[index];
        if (block.valid && block.tag == tag) {
            if (operation == 'S') {
                block.data[blockOffset] = 'X';
            }
                hits++;
            } else {
            misses++;
            block.tag = tag;
            block.valid = true;
            }
        }

    outfile << hits << "," <<  hits + misses << "; ";
    outfile.close();
}

void SetAssociativeCache(const vector<pair<char, unsigned int> > &traces, int associativity, const string& outputFile) {
    int lineSize = 32;
    int cacheSize = 16384;
    int numSets = cacheSize / (lineSize * associativity);
    vector<vector<CacheBlock> > cache(numSets, vector<CacheBlock>(associativity, CacheBlock(lineSize)));
    int numLines = cacheSize / lineSize;
    int hits = 0;
    int misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;
        int blockOffsetBits = log2(32);
        int setIndexBits = log2(numSets);
        int blockOffset = address & (lineSize - 1);
        int setindex = (address / lineSize) & (numSets - 1);
        int tag = address >> (blockOffsetBits + setIndexBits);
        vector<CacheBlock> &set = cache[setindex];
        bool hit = false;

        for (int j = 0; j < associativity; ++j) {
            if (set[j].valid && set[j].tag == tag) { // hit
                set[j].lru_counter = 0;
                for (int k = 0; k < set.size(); ++k) {
                    if (k != j && set[k].valid) {
                        set[k].lru_counter++;
                    }
                }
                if (operation == 'S') {
                    set[j].data[blockOffset] = 'X';
                }
                hit = true;
                hits++;
                break;
            }
        }

        if (!hit) { //miss
            misses++;
            if (isSetFull(set)) {
                int lruIndex = findLRUBlockIndex(set);
                set[lruIndex].tag = tag;
                set[lruIndex].valid = true;
                set[lruIndex].lru_counter = 0;
                if (operation == 'S') {
                    set[lruIndex].data[blockOffset] = 'X';
                }
                for (int k = 0; k < set.size(); ++k) {
                    if (k != lruIndex && set[k].valid) {
                        set[k].lru_counter++;
                    }
                }
            } else {
           for (CacheBlock &block : set) {
                if (!block.valid) {
                    block.tag = tag;
                    block.valid = true;
                    if (operation == 'S') {
                        block.data[blockOffset] = 'X';
                    }
                    break; 
                }
            }
        }
    } 
    }

    ofstream outfile(outputFile, ios::app);
    outfile << hits << "," << (hits + misses) << "; ";
    outfile.close();
}

void FullyAssociativeCacheLRU(const vector<pair<char, unsigned int> > &traces,  const string& outputFile) {
    int lineSize = 32;
    int cacheSize = 16384;
    int numSets = cacheSize / (lineSize);
    vector<CacheBlock> cache(numSets, CacheBlock(lineSize));
    int numLines = cacheSize / lineSize;
    int hits = 0;
    int misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;
        int blockOffset = address & (lineSize - 1);
        int tag = address >> 5;
        bool hit = false;

        for (int j = 0; j < cache.size(); ++j) {
            if (cache[j].valid && cache[j].tag == tag) { // hit
                cache[j].lru_counter = 0;
                for (int k = 0; k < cache.size(); ++k) {
                    if (k != j && cache[k].valid) {
                        cache[k].lru_counter++;
                    }
                }
                if (operation == 'S') {
                    cache[j].data[blockOffset] = 'X';
                }
                hit = true;
                hits++;
                break;
            }
        }

        if (!hit) { //miss
            misses++;
            if (isSetFull(cache)) {
                int lruIndex = findLRUBlockIndex(cache);
                cache[lruIndex].tag = tag;
                cache[lruIndex].valid = true;
                cache[lruIndex].lru_counter = 0;
                if (operation == 'S') {
                    cache[lruIndex].data[blockOffset] = 'X';
                }
                for (int k = 0; k < cache.size(); ++k) {
                    if (k != lruIndex && cache[k].valid) {
                        cache[k].lru_counter++;
                    }
                }
            } else {
                for (CacheBlock &block : cache) {
                    if (!block.valid) {
                        block.tag = tag;
                        block.valid = true;
                            if (operation == 'S') {
                            block.data[blockOffset] = 'X';
                        }
                        break; 
                    }
                }
            }
        } 
    }
    ofstream outfile(outputFile, ios::app);
    outfile << hits << "," << (hits + misses) << "; ";
    outfile.close();
}

void FullyAssociativeCacheHotCold(const vector<pair<char, unsigned int> > &traces, const string& outputFile) {
    int lineSize = 32;
    int cacheSize = 16384;
    int numSets = cacheSize / lineSize;
    vector<CacheBlock> cache(numSets, CacheBlock(lineSize));
    int numLines = cacheSize / lineSize;
    int hits = 0;
    int misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;
        int blockOffset = address & (lineSize - 1);
        int tag = address >> 5;
        bool hit = false;

        for (int j = 0; j < cache.size(); ++j) {
            if (cache[j].valid && cache[j].tag == tag) { // hit
                cache[j].hot_cold_bit = true;
                if (operation == 'S') {
                    cache[j].data[blockOffset] = 'X';
                }
                hit = true;
                hits++;
                break;
            }
        }

        if (!hit) {
            misses++;
            if (isSetFull(cache)) {
                int coldIdx = findFirstColdBlock(cache);
                if(coldIdx == -1){
                    for (int k = 0; k < cache.size(); ++k) {
                        cache[k].hot_cold_bit = false;    
                    }
                    coldIdx = findFirstColdBlock(cache);
                } 
                if (coldIdx != -1) {
                    cache[coldIdx].tag = tag;
                    cache[coldIdx].valid = true;
                    cache[coldIdx].hot_cold_bit = true;
                    if (operation == 'S') {
                        cache[coldIdx].data[blockOffset] = 'X';
                    }
                }
            } else {
                    for (CacheBlock &block : cache) {
                    if (!block.valid) {
                     block.tag = tag;
                     block.valid = true;
                     block.hot_cold_bit = true;
                     if (operation == 'S') {
                         block.data[blockOffset] = 'X';
                     }
                     break; 
                    }
                }   
            }
        } 
    }
    ofstream outfile(outputFile, ios::app);
    outfile << hits << "," << (hits + misses) << "; ";
    outfile.close();
}

void TwoLevelCache(const vector<pair<char, unsigned int> > &traces, const string& outputFile) {
    int l1LineSize = 32;
    int l1CacheSize = 4096;
    int l1Associativity = 4;

    int l2Associativity = 8;
    int l2CacheSize = 65536;
    int l2LineSize = 64;
    
    int l1Numsets = l1CacheSize / (l1LineSize * l1Associativity);
    int l2Numsets = l2CacheSize / (l2LineSize * l2Associativity);

    vector<vector<CacheBlock> > l1cache(l1Numsets, vector<CacheBlock>(l1Associativity, CacheBlock(l1LineSize)));
    vector<vector<CacheBlock> > l2cache(l2Numsets, vector<CacheBlock>(l2Associativity, CacheBlock(l2LineSize)));

    int l1Hits = 0;
    int l2Hits = 0;
    int l1Misses = 0;
    int l2Misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;

        int l1blockOffsetBits = log2(l1LineSize);
        int l1setIndexBits = log2(l1Numsets);
        int l1blockOffset = address & (l1LineSize - 1);
        int l1setindex = (address / l1LineSize) & (l1Numsets - 1);
        int l1tag = address >> (l1blockOffsetBits + l1setIndexBits);
        vector<CacheBlock> &l1set = l1cache[l1setindex];
        bool l1hit = false;

        int l2blockOffsetBits = log2(l2LineSize);
        int l2setIndexBits = log2(l2Numsets);
        int l2blockOffset = address & (l2LineSize - 1);
        int l2setindex = (address / l2LineSize) & (l2Numsets - 1);
        int l2tag = address >> (l2blockOffsetBits + l2setIndexBits);
        vector<CacheBlock> &l2set = l2cache[l2setindex];
        bool l2hit = false;

        for (int j = 0; j < l1Associativity; ++j) {
            if (l1set[j].valid && l1set[j].tag == l1tag) { // hit
                l1set[j].lru_counter = 0;
                for (int k = 0; k < l1set.size(); ++k) {
                    if (k != j && l1set[k].valid) {
                        l1set[k].lru_counter++;
                    }
                }
                    if (operation == 'S') {
                        for (int l = 0; j < l2Associativity; ++l) {
                        if (l2set[l].valid && l2set[l].tag == l2tag) { // hit
                        l2set[l].lru_counter = 0;
                        for (int m = 0; m < l2set.size(); ++m) {
                        if (m != l && l2set[m].valid) {
                            l2set[m].lru_counter++;
                        }
                    }
                    l2hit = true;
                    l2Hits++;
                    break;
                }
            }
                    }
                
                l1hit = true;
                l1Hits++;
                break;
            }
        }

        if (!l1hit) { //miss
            l1Misses++;

            for (int j = 0; j < l2Associativity; ++j) {
                if (l2set[j].valid && l2set[j].tag == l2tag) { // hit
                    l2set[j].lru_counter = 0;
                    for (int k = 0; k < l2set.size(); ++k) {
                        if (k != j && l2set[k].valid) {
                            l2set[k].lru_counter++;
                        }
                    }
                    l2hit = true;
                    l2Hits++;
                    break;
                }
            }
            if (!l2hit) { //miss
            l2Misses++;
            if (isSetFull(l2set)) {
                int lruIndex = findLRUBlockIndex(l2set);
                l2set[lruIndex].tag = l2tag;
                l2set[lruIndex].valid = true;
                l2set[lruIndex].lru_counter = 0;

                for (int k = 0; k < l2set.size(); ++k) {
                    if (k != lruIndex && l2set[k].valid) {
                        l2set[k].lru_counter++;
                    }
                }
            } else {
                for (CacheBlock &block : l2set) {
                    if (!block.valid) {
                        block.tag = l2tag;
                        block.valid = true;
                        break; 
                    }
                }
            }
        } 
        if (isSetFull(l1set)) {
                int lruIndex = findLRUBlockIndex(l1set);
                l1set[lruIndex].tag = l1tag;
                l1set[lruIndex].valid = true;
                l1set[lruIndex].lru_counter = 0;
                for (int k = 0; k < l1set.size(); ++k) {
                    if (k != lruIndex && l1set[k].valid) {
                        l1set[k].lru_counter++;
                    }
                }
            } else {
                for (CacheBlock &block : l1set) {
                    if (!block.valid) {
                        block.tag = l1tag;
                        block.valid = true;
                        break; 
                    }
                }
            }
        } 
    }

    ofstream outfile(outputFile, ios::app);
    outfile << l1Hits << "," << (l1Hits + l1Misses) << ";";
    outfile << l2Hits << "," << (l2Hits + l2Misses) << "; ";
    outfile.close();
}

void TwoLevelCacheWB(const vector<pair<char, unsigned int> > &traces, const string& outputFile) {
    int l1LineSize = 32;
    int l1CacheSize = 4096;
    int l1Associativity = 4;

    int l2Associativity = 8;
    int l2CacheSize = 65536;
    int l2LineSize = 64;
    
    int l1Numsets = l1CacheSize / (l1LineSize * l1Associativity);
    int l2Numsets = l2CacheSize / (l2LineSize * l2Associativity);

    vector<vector<CacheBlock> > l1cache(l1Numsets, vector<CacheBlock>(l1Associativity, CacheBlock(l1LineSize)));
    vector<vector<CacheBlock> > l2cache(l2Numsets, vector<CacheBlock>(l2Associativity, CacheBlock(l2LineSize)));

    int l1Hits = 0;
    int l2Hits = 0;
    int l1Misses = 0;
    int l2Misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;

        int l1blockOffsetBits = log2(l1LineSize);
        int l1setIndexBits = log2(l1Numsets);
        int l1blockOffset = address & (l1LineSize - 1);
        int l1setindex = (address / l1LineSize) & (l1Numsets - 1);
        int l1tag = address >> (l1blockOffsetBits + l1setIndexBits);
        vector<CacheBlock> &l1set = l1cache[l1setindex];
        bool l1hit = false;

        int l2blockOffsetBits = log2(l2LineSize);
        int l2setIndexBits = log2(l2Numsets);
        int l2blockOffset = address & (l2LineSize - 1);
        int l2setindex = (address / l2LineSize) & (l2Numsets - 1);
        int l2tag = address >> (l2blockOffsetBits + l2setIndexBits);
        vector<CacheBlock> &l2set = l2cache[l2setindex];
        bool l2hit = false;

        for (int j = 0; j < l1Associativity; ++j) {
            if (l1set[j].valid && l1set[j].tag == l1tag) { // l1 hit
                l1set[j].lru_counter = 0;            
                for (int k = 0; k < l1Associativity; ++k) {//age rest of blocks in set
                    if (k != j && l1set[k].valid) {
                        l1set[k].lru_counter++;
                    }
                }
                l1hit = true;                
                l1set[j].dirty = (operation == 'S'); 
                l1Hits++;
                break;
            }
        }
       if (!l1hit) { //l1 misses
            l1Misses++;
            for (int j = 0; j < l2Associativity; ++j) {// access
                if (l2set[j].valid && l2set[j].tag == l2tag) { //l2 hit
                        l2set[j].lru_counter = 0;
                        for (int k = 0; k < l2Associativity; ++k) { // update the LRU counters
                        if (k != j && l2set[k].valid) {
                            l2set[k].lru_counter++;
                        }
                    }
                    l2hit = true;
                    l2Hits++;
                    l2set[j].dirty = (operation == 'S');                    
                    break;
                }
            }

            if (!l2hit) { // l2 misses
                l2Misses++;
                if (isSetFull(l2set)) { //every block in set is valid
                    int lruIndex = findLRUBlockIndex(l2set); //write over lru
                    if (l2set[lruIndex].dirty) {
                        int l1lruIndex = findLRUBlockIndex(l1set);
                        // write-back if the block is dirty
                        l1set[l1lruIndex].dirty = true;
                    }
                    l2set[lruIndex].tag = l2tag;
                    l2set[lruIndex].valid = true;
                    l2set[lruIndex].dirty = (operation == 'S');
                } else {
                    for (auto &block : l2set) { 
                        if (!block.valid) {
                            block.tag = l2tag;
                            block.valid = true;
                            block.dirty = (operation == 'S');
                            break;
                        }
                    }
                }
            }
            if (isSetFull(l1set)) {
                int lruIndex = findLRUBlockIndex(l1set);
                if (l1set[lruIndex].dirty) {
                    // write-back dirty block to L2 if needed
                    l2set[lruIndex].dirty = true;
                }
                l1set[lruIndex].tag = l1tag;
                l1set[lruIndex].valid = true;
                l1set[lruIndex].dirty = (operation == 'S');
            } else {
                for (auto &block : l1set) {
                    if (!block.valid) {
                        block.tag = l1tag;
                        block.valid = true;
                        block.dirty = (operation == 'S');
                        break;
                    }
                }
            }
        }
    }           
    ofstream outfile(outputFile, ios::app);
    double utilization = (double)l2Hits / (l2Hits + l2Misses);
    outfile << l1Hits << "," << (l1Hits + l1Misses) << ";";
    outfile << l2Hits << "," << (l2Hits + l2Misses) << ";";
    outfile << fixed << setprecision(5) << utilization << "; ";
    outfile.close();
}

void Utilization(const vector<pair<char, unsigned int> > &traces, const string& outputFile) {
    int l1LineSize = 32;
    int l1CacheSize = 4096;
    int l1Associativity = 4;

    int l2Associativity = 8;
    int l2CacheSize = 65536;
    int l2LineSize = 64;
    
    int l1Numsets = l1CacheSize / (l1LineSize * l1Associativity);
    int l2Numsets = l2CacheSize / (l2LineSize * l2Associativity);

    vector<CacheBlock> cachel1(l1Numsets, CacheBlock(l1LineSize));
    vector<CacheBlock> cachel2(l2Numsets, CacheBlock(l2LineSize));

    int l1Hits = 0;
    int l2Hits = 0;
    int l1Misses = 0;
    int l2Misses = 0;

    for (int i = 0; i < traces.size(); ++i) { //gathers information of current address within trace
        char operation = traces[i].first;
        int address = traces[i].second;

        int l1blockOffsetBits = log2(l1LineSize);
        int l1setIndexBits = log2(l1Numsets);
        int l1blockOffset = address & (l1LineSize - 1);
        int l1setindex = (address / l1LineSize) & (l1Numsets - 1);
        int l1tag = address >> (l1blockOffsetBits + l1setIndexBits);
        bool l1hit = false;

        int l2blockOffsetBits = log2(l2LineSize);
        int l2setIndexBits = log2(l2Numsets);
        int l2blockOffset = address & (l2LineSize - 1);
        int l2setindex = (address / l2LineSize) & (l2Numsets - 1);
        int l2tag = address >> (l2blockOffsetBits + l2setIndexBits);
        bool l2hit = false;

        for (int j = 0; j < l1Associativity; ++j) {
            if (cachel1[j].valid && cachel1[j].tag == l1tag) { // hit
                l1hit = true;
                l1Hits++;
                for (int k = 0; k < cachel1.size(); ++k) {
                    if (k != j && cachel1[k].valid) {
                        cachel1[k].lru_counter++;
                    }
                }
                cachel1[j].lru_counter = 0;
                if (operation == 'L') {
                    cachel1[j].dirty = true;
                }
                break;
            }

       if (!l1hit) {
            l1Misses++;
            for (int j = 0; j < l2Associativity; ++j) {
                if (cachel2[j].valid && cachel2[j].tag == l2tag) {
                    l2hit = true;
                    l2Hits++;
                    cachel2[j].lru_counter = 0;
                    if (operation == 'S') {
                        cachel2[j].dirty = true;
                    }
                        for (auto &block : cachel2) {
                        if (block.valid && block.tag != l2tag) {
                            block.lru_counter++;
                        }
                    }

                    if (isSetFull(cachel1)) {
                        int lruIndex = findLRUBlockIndex(cachel1);

                        if (cachel1[lruIndex].dirty) {
                            int evictedAddress = (cachel1[lruIndex].tag << (l1blockOffsetBits + l1setIndexBits)) | (l1setindex << l1blockOffsetBits);
                            int evictedL2SetIndex = (evictedAddress / l2LineSize) & (l2Numsets - 1);
                            int evictedL2Tag = evictedAddress >> (l2blockOffsetBits + l2setIndexBits);
                            for (auto &block : cachel2) {
                                if (block.valid && block.tag == evictedL2Tag) {
                                    block.dirty = true;
                                    break;
                                }
                            }
                        }

                        cachel1[lruIndex].tag = l1tag;
                        cachel1[lruIndex].valid = true;
                        cachel1[lruIndex].dirty = false;
                    } else {
                        for (auto &block : cachel1) {
                            if (!block.valid) {
                            block.tag = l1tag;
                            block.valid = true;
                            block.dirty = (operation == 'S');
                            break;
                            }
                        }
                    }
                    break;
                }
            }

            if (!l2hit) {
                l2Misses++;

                if (isSetFull(cachel2)) {
                    int lruIndex = findLRUBlockIndex(cachel2);

                    if (cachel2[lruIndex].dirty) {
                    }

                    cachel2[lruIndex].tag = l2tag;
                    cachel2[lruIndex].valid = true;
                    cachel2[lruIndex].dirty = false;
                } else {
                    for (auto &block : cachel2) {
                        if (!block.valid) {
                            block.tag = l2tag;
                            block.valid = true;
                            block.dirty = false;
                            break;
                        }
                    }
                }



            }
        }
    }
    }

    ofstream outfile(outputFile, ios::app);
    double utilization = (double)l2Hits / (l2Hits + l2Misses);
    outfile<< l1Hits << "," << (l1Hits + l1Misses) << ";";
    outfile << l2Hits << "," << (l2Hits + l2Misses) << ";";
    outfile << utilization << "; ";
    outfile.close();
    }

int main(int argc, char *argv[]) {
    vector<pair<char, unsigned int> > traces = parseTraceFile(argv[1]);
    DirectMappedCache(traces, 1024, argv[2]);
    DirectMappedCache(traces, 4096, argv[2]);
    DirectMappedCache(traces, 16384, argv[2]);
    DirectMappedCache(traces, 32768, argv[2]);
    ofstream outfile1(argv[2], ios::app);
    outfile1 << endl;
    outfile1.close();
    SetAssociativeCache(traces, 2, argv[2]);
    SetAssociativeCache(traces, 4, argv[2]);
    SetAssociativeCache(traces, 8, argv[2]);
    SetAssociativeCache(traces, 16, argv[2]);
    ofstream outfile2(argv[2], ios::app);
    outfile2 << endl;
    outfile2.close();
    FullyAssociativeCacheLRU(traces, argv[2]);
    FullyAssociativeCacheHotCold(traces, argv[2]);
    ofstream outfile3(argv[2], ios::app);
    outfile3 << endl;
    outfile3.close();
    TwoLevelCache(traces, argv[2]);
    ofstream outfile4(argv[2], ios::app);
    outfile4 << endl;
    outfile4.close();
    TwoLevelCacheWB(traces, argv[2]);
    ofstream outfile5(argv[2], ios::app);
    outfile5 << endl;
    outfile5.close();
    Utilization(traces, argv[2]);
    return 0;
}