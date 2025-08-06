#include<bits/stdc++.h>
#include "MurmurHash3.h"
#include "Tool.cpp"
using namespace std;
class CellOfBkt_P2S_FC {
public:
    uint32_t flowLabel;
    float est;    
    int priority;
    CellOfBkt_P2S_FC(uint32_t flowLabel, float est, int priority)
        : flowLabel(flowLabel), est(est), priority(priority) {}
    bool operator>(const CellOfBkt_P2S_FC& other) const {
        return priority > other.priority;
    }
    bool operator<(const CellOfBkt_P2S_FC& other) const {
        return priority < other.priority;
    }
};
class Bkt_P2S_FC {
public:
    int sizeOfBkt;
    vector<CellOfBkt_P2S_FC> cellArray;
    Bkt_P2S_FC(int sizeOfBkt) : sizeOfBkt(sizeOfBkt) {}
    void sortBkt() {
        sort(cellArray.begin(), cellArray.end());
    }
    int findMinPriorityInBkt() {
        int min_idx = 0;
        int min_priority = cellArray[0].priority;
        for (size_t idx = 1; idx < sizeOfBkt; ++idx) {
            if (cellArray[idx].priority < min_priority) {
                min_priority = cellArray[idx].priority;
                min_idx = idx;
            }
        }
        return min_idx;
    }
    int posOfFlowInBkt(uint32_t flowLabel) {
        for (size_t idx = 0; idx < cellArray.size(); ++idx) {
            if (cellArray[idx].flowLabel == flowLabel) {
                return static_cast<int>(idx);
            }
        }
        return -1;
    }
    bool isBktFull() {
        return cellArray.size() == static_cast<size_t>(sizeOfBkt);
    }
};
class P2S_FC {
public:
    int gamma;
    vector<int> compressRatio;
    uint32_t compressSeed;
    int bitmapLen;
    vector<bool> bitmap;
    uint32_t bitmapSeed;
    int zeroNumOfBitmap;
    int numOfBkt;
    int sizeOfBkt;
    vector<Bkt_P2S_FC> bktArray;
    uint32_t seed4ChoosingBkt;
    int d;
    int w;
    vector<vector<int>> CS;
    vector<uint32_t> hash_CS_seeds;
    P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step)
        : gamma(gamma), bitmapLen(bitmapLen), numOfBkt(numOfBkt), sizeOfBkt(sizeOfBkt), d(d), w(w) {
        bitmap.resize(bitmapLen, false);
        zeroNumOfBitmap = bitmapLen;
        vector<uint32_t> seeds = generateSeeds(3);
        bitmapSeed = seeds[0];
        seed4ChoosingBkt = seeds[1];
        compressSeed = seeds[2];
        compressRatio.push_back(0);
        for (int i = 1; i <= gamma; ++i) {
            compressRatio.push_back(static_cast<int>(base * pow(static_cast<double>(i), step)));
        }
        cout << "base: " << base << " step: " << step << endl;
        for (int i = 0; i < compressRatio.size(); i++) {
            cout << compressRatio[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < numOfBkt; ++i) {
            bktArray.emplace_back(sizeOfBkt);
        }
        CS.resize(d, vector<int>(w, 0));
        hash_CS_seeds = generateSeeds(2 * d);
    }
    void update(uint32_t src, uint32_t dst, int priority) {
        uint32_t hash_value_dst;
        MurmurHash3_x86_32(reinterpret_cast<const void*>(&dst), sizeof(dst), compressSeed, &hash_value_dst);
        uint32_t elem4Compress = hash_value_dst % static_cast<uint32_t>(compressRatio[priority]);
        uint32_t combined_flow_id_parts = (src ^ elem4Compress); 
        uint32_t hash_value_bitmap;
        MurmurHash3_x86_32(reinterpret_cast<const void*>(&combined_flow_id_parts), 
                           sizeof(combined_flow_id_parts), bitmapSeed, &hash_value_bitmap);
        uint32_t hash_idx = hash_value_bitmap % bitmapLen;
        bool flag = false;
        if (!bitmap[hash_idx]) {
            zeroNumOfBitmap--;
            flag = true;
            bitmap[hash_idx] = true;
        }
        if (!flag) {
            return ;
        }
        float estIncrement = static_cast<float>(bitmapLen * 1.0 / zeroNumOfBitmap);
        if (priority > gamma * 0.6) {
            uint32_t hash_value_bkt;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seed4ChoosingBkt, &hash_value_bkt);
            uint32_t idx4Bkt = hash_value_bkt % numOfBkt;
            int posInTheBkt = bktArray[idx4Bkt].posOfFlowInBkt(src);
            if (posInTheBkt > -1) {
                bktArray[idx4Bkt].cellArray[posInTheBkt].est += estIncrement;
            } else {
                if (bktArray[idx4Bkt].cellArray.size() == static_cast<size_t>(sizeOfBkt)) {
                    int min_idx = bktArray[idx4Bkt].findMinPriorityInBkt();
                    if (priority > bktArray[idx4Bkt].cellArray[min_idx].priority) {
                        CellOfBkt_P2S_FC tempCell(src, estIncrement, priority);
                        coarsePartInsert(bktArray[idx4Bkt].cellArray[min_idx].flowLabel, static_cast<int>(bktArray[idx4Bkt].cellArray[min_idx].est));
                        bktArray[idx4Bkt].cellArray[min_idx] = tempCell;
                    } else {
                        coarsePartInsert(src, static_cast<int>(estIncrement));
                    }
                } else {
                    CellOfBkt_P2S_FC tempCell(src, estIncrement, priority);
                    bktArray[idx4Bkt].cellArray.push_back(tempCell);
                }
            }
        } else {
            coarsePartInsert(src, static_cast<int>(estIncrement));
        }
    }
    int estimate(uint32_t src, int priority) {
        uint32_t hash_value_bkt;
        MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seed4ChoosingBkt, &hash_value_bkt);
        uint32_t idx4Bkt = hash_value_bkt % numOfBkt;
        int posInTheBkt = bktArray[idx4Bkt].posOfFlowInBkt(src);
        int tmpEst = 0;
        if (posInTheBkt > -1) {
            tmpEst = static_cast<int>(bktArray[idx4Bkt].cellArray[posInTheBkt].est);
        } else {
            tmpEst = coarsePartQuery(src);
        }
        int m = compressRatio[priority];
        int num_0 = max(m - tmpEst, 1);
        int est = static_cast<int>(1.0 * -1 * m * log(1.0 * num_0 / m));
        return est;
    }
    void coarsePartInsert(uint32_t src, int estIncrement) {
        for (int i = 0; i < d; ++i) {
            uint32_t hash_value_cs_idx;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), hash_CS_seeds[i], &hash_value_cs_idx);
            uint32_t hash_idx = hash_value_cs_idx % w;
            uint32_t hash_value_cs_val2;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), hash_CS_seeds[d + i], &hash_value_cs_val2);
            uint32_t hash_val2 = hash_value_cs_val2 % 2;
            int coe = (hash_val2 == 1) ? 1 : -1;
            CS[i][hash_idx] += coe * estIncrement;
        }
    }
    int coarsePartQuery(uint32_t src) {
        vector<int> val_lst;
        for (int i = 0; i < d; ++i) {
            uint32_t hash_value_cs_idx;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), hash_CS_seeds[i], &hash_value_cs_idx);
            uint32_t hash_idx = hash_value_cs_idx % w;
            val_lst.push_back(CS[i][hash_idx]);
        }
        sort(val_lst.begin(), val_lst.end());
        size_t middle_idx = val_lst.size() / 2;
        int result;
        if (val_lst.size() % 2 == 0) {
            result = (val_lst[middle_idx] + val_lst[middle_idx - 1]) / 2; 
        } else {
            result = val_lst[middle_idx];
        }
        return max(result, 1); 
    }
};
class Experiment_P2S_FC {
public:
    int gamma;
    P2S_FC sketch;
    unordered_map<uint32_t, set<uint32_t>> real_set_dict;
    unordered_map<uint32_t, uint32_t> real_dict; 
    unordered_map<uint32_t, uint32_t> est_dict;   
    unordered_map<uint32_t, int> flowPriority_dict;
    Experiment_P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step)
        : gamma(gamma), sketch(bitmapLen, numOfBkt, sizeOfBkt, d, w, gamma, base, step) {}
    void start(const string& filename, vector<vector<double>>& results, int idx) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }
        vector<tuple<uint32_t, uint32_t, int>> parsed_data;
        string line;
        while (getline(file, line)) {
            size_t first_space = line.find(' ');
            size_t second_space = line.find(' ', first_space + 1);
            string src_str = line.substr(0, first_space);
            string dst_str = line.substr(first_space + 1, second_space - first_space - 1);
            string priority_str = line.substr(second_space + 1);
            uint32_t src = stoul(src_str);
            uint32_t dst = stoul(dst_str);
            int priority = stoi(priority_str);
            parsed_data.emplace_back(src, dst, priority);
        }
        file.close();
        auto start_time = chrono::high_resolution_clock::now();
        for (const auto& [src, dst, priority] : parsed_data) {
            sketch.update(src, dst, priority);
        }
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        double throughput = static_cast<double>(parsed_data.size()) / elapsed_time.count() / 1e6;
        cout << "throughput : " << throughput << " MPPS" << endl;
        results[4][idx] = throughput;
        for (const auto& [src, dst, priority] : parsed_data) {
            real_set_dict[src].insert(dst);
            flowPriority_dict[src] = priority;
        }
    }
    void estimate() {
        int progress = 0;
        int total_flows = real_set_dict.size();
        int count = 0;
        int update_interval = total_flows / 100 > 0 ? total_flows / 100 : 1;
        for (const auto& pair : real_set_dict) {
            const uint32_t& src = pair.first;
            real_dict[src] = real_set_dict[src].size();
            est_dict[src] = sketch.estimate(src, flowPriority_dict[src]);
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages]: Success: " << progress << "%" << flush;
            }
        }
    }
    void performance(vector<vector<double>>& results, int idx) {
        double highPriorityErrorSUM = 0;
        int highPriorityFlowNum = 0;
        double lowPriorityErrorSUM = 0;
        int lowPriorityFlowNum = 0;
        int progress = 0;
        int total_flows = real_dict.size();
        int count = 0;
        int update_interval = total_flows / 100 > 0 ? total_flows / 100 : 1;
        for (const auto& pair : real_dict) {
            const uint32_t& src = pair.first;
            if (flowPriority_dict[src] > gamma * 0.75) {
                highPriorityFlowNum++;
                highPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - static_cast<double>(est_dict[src])) / real_dict[src];
            } else {
                lowPriorityFlowNum++;
                lowPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - static_cast<double>(est_dict[src])) / real_dict[src];
            }
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages]: Success: " << progress << "%" << flush;
            }
        }
        cout << endl;
        double highPriorityARE = highPriorityErrorSUM / highPriorityFlowNum;
        double lowPriorityARE = lowPriorityErrorSUM / lowPriorityFlowNum;
        double allFlowsARE = (highPriorityErrorSUM + lowPriorityErrorSUM) / (highPriorityFlowNum + lowPriorityFlowNum);
        int countHigh = 0;
        for (size_t i = 0; i < sketch.bktArray.size(); ++i) {
            for (size_t j = 0; j < sketch.bktArray[i].cellArray.size(); ++j) {
                if (sketch.bktArray[i].cellArray[j].priority > gamma * 0.75) {
                    countHigh++;
                }
            }
        }
        double pr = 1.0;
        double rr = static_cast<double>(countHigh) / highPriorityFlowNum;
        double f1 = 2 * pr * rr / (pr + rr);
        cout << "f1-score : " << fixed << setprecision(6) << f1 << endl;
        cout << "highPriorityARE : " << fixed << setprecision(6) << highPriorityARE << endl;
        cout << "lowPriorityARE : " << fixed << setprecision(6) << lowPriorityARE << endl;
        cout << "allFlowsARE : " << fixed << setprecision(6) << allFlowsARE << endl;
        results[0][idx] = highPriorityARE;
        results[1][idx] = lowPriorityARE;
        results[2][idx] = allFlowsARE;
        results[3][idx] = f1;
    }
    void pefmance2() {
        double highPriorityErrorSUM = 0;
        int highPriorityFlowNum = 0;
        int progress = 0;
        int total_flows = real_dict.size();
        int count = 0;
        int update_interval = total_flows / 100 > 0 ? total_flows / 100 : 1;
        for (const auto& pair : real_dict) {
            const uint32_t& src = pair.first;
            if (flowPriority_dict[src] > gamma * 0.75) {
                uint32_t hash_value_bkt;
                MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), sketch.seed4ChoosingBkt, &hash_value_bkt);
                uint32_t idx4Bkt = hash_value_bkt % sketch.numOfBkt;
                int posInTheBkt = sketch.bktArray[idx4Bkt].posOfFlowInBkt(src);
                if (posInTheBkt > -1) {
                    uint32_t estTemp = sketch.bktArray[idx4Bkt].cellArray[posInTheBkt].est;
                    highPriorityFlowNum++;
                    highPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - estTemp) / real_dict[src];
                }
            }
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages]: Success: " << progress << "%" << flush;
            }
        }
    }
};
vector<double> Exp_P2S_FC(double alpha, int gamma, double wholeMemory, string priority_file_name) {
    cout.precision(10);
    vector<vector<double>> results(5, vector<double>(1)); 
    int threeWholeMem = wholeMemory;
    double finePartMem = 0.2 * threeWholeMem;
    double wholeMem = 0.8 * threeWholeMem; 
    double bitmapMemRatio = 0.5;
    int bits = static_cast<int>(wholeMem * bitmapMemRatio * 1024 * 8);
    int sizeOfBkt = 4;
    int numOfBkt = static_cast<int>(finePartMem * 1024 / 9 / sizeOfBkt); 
    int d = 3;
    int w = static_cast<int>(wholeMem * (1 - bitmapMemRatio) * 1024 / (2 * d));
    double base = 100; 
    double step = 2; 
    Experiment_P2S_FC exp(bits, numOfBkt, sizeOfBkt, d, w, gamma, base, step);
    stringstream ss_alpha;
    ss_alpha << fixed << setprecision(6) << alpha;
    string alpha_str = ss_alpha.str();
    string filename = priority_file_name + "_" + to_string(alpha) + "_" + to_string(gamma) + "_.txt"; 
    exp.start(filename, results, 0);
    exp.estimate();
    exp.performance(results, 0);
    exp.pefmance2();
    vector<double> ret(5);
    for (int i = 0; i < 5; ++i) {
        ret[i] = results[i][0];
    }
    return ret;
}
