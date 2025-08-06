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
struct Superspread_P2S_FC {
    uint32_t flowLabel;
    int priority;
    Superspread_P2S_FC(uint32_t f = 0, int p = 0) : flowLabel(f), priority(p) {}
    bool operator<(const Superspread_P2S_FC& other) const {
        return flowLabel > other.flowLabel; 
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
    set<Superspread_P2S_FC> report_Superspread;
    int threshold_Superspread;
    P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step, int threshold_Superspread)
        : gamma(gamma), bitmapLen(bitmapLen), numOfBkt(numOfBkt), sizeOfBkt(sizeOfBkt), d(d), w(w), threshold_Superspread(threshold_Superspread) {
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
        uint64_t combined_flow_id_parts = static_cast<uint64_t>(src);
        combined_flow_id_parts = (combined_flow_id_parts << 32) | elem4Compress;
        uint32_t hash_value_bitmap;
        MurmurHash3_x86_32(reinterpret_cast<const void*>(&combined_flow_id_parts), 
                           sizeof(combined_flow_id_parts), bitmapSeed, &hash_value_bitmap);
        uint32_t hash_idx = hash_value_bitmap % bitmapLen;
        bool flag = false; 
        if (!bitmap[hash_idx]) {
            flag = true;
            if (zeroNumOfBitmap > 1) {
                zeroNumOfBitmap--;
                bitmap[hash_idx] = true;
            }
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
                if (priority > 0.75 * gamma && bktArray[idx4Bkt].cellArray[posInTheBkt].est >= threshold_Superspread) {
                    report_Superspread.insert(Superspread_P2S_FC(src, priority));
                }
            } else {
                if (bktArray[idx4Bkt].cellArray.size() == static_cast<size_t>(sizeOfBkt)) {
                    int min_idx = bktArray[idx4Bkt].findMinPriorityInBkt();
                    if (priority > bktArray[idx4Bkt].cellArray[min_idx].priority) {
                        CellOfBkt_P2S_FC replaced_cell = bktArray[idx4Bkt].cellArray[min_idx]; 
                        CellOfBkt_P2S_FC tempCell(src, estIncrement, priority);
                        bktArray[idx4Bkt].cellArray[min_idx] = tempCell; 
                        coarsePartInsert(replaced_cell.flowLabel, static_cast<int>(replaced_cell.est)); 
                        if (priority > 0.75 * gamma && estIncrement >= threshold_Superspread) {
                            report_Superspread.insert(Superspread_P2S_FC(src, priority));
                        }
                    } else {
                        coarsePartInsert(src, static_cast<int>(estIncrement));
                    }
                } else {
                    CellOfBkt_P2S_FC tempCell(src, estIncrement, priority);
                    bktArray[idx4Bkt].cellArray.push_back(tempCell);
                    if (priority > 0.75 * gamma && estIncrement >= threshold_Superspread) {
                        report_Superspread.insert(Superspread_P2S_FC(src, priority));
                    }
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
    set<Superspread_P2S_FC> getSuperSpread() const {
        return report_Superspread;
    }
};
class Experiment_P2S_FC {
public:
    int gamma;
    int threshold_Superspread;
    P2S_FC sketch; 
    unordered_map<uint32_t, set<uint32_t>> real_set_dict;
    unordered_map<uint32_t, uint32_t> real_dict; 
    unordered_map<uint32_t, uint32_t> est_dict;   
    unordered_map<uint32_t, int> flowPriority_dict;
    Experiment_P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step, int threshold_Superspread)
        : gamma(gamma), threshold_Superspread(threshold_Superspread), sketch(bitmapLen, numOfBkt, sizeOfBkt, d, w, gamma, base, step, threshold_Superspread) {}
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
                cout << "\r[Messages] Success:" << progress << "%" << flush;
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
                highPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - est_dict[src]) / real_dict[src];
            } else {
                lowPriorityFlowNum++;
                lowPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - est_dict[src]) / real_dict[src];
            }
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages] Success:" << progress << "%" << flush;
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
                cout << "\r[Messages] Success:" << progress << "%" << flush;
            }
        }
    }
    void performance_Superspread(vector<vector<double>>& results, int idx) {
        set<Superspread_P2S_FC> true_superspread_flows;
        cout << real_dict.size() << endl;
        cout << threshold_Superspread << endl;
        for (const auto& pair : real_dict) { 
            uint32_t src_flow = pair.first;
            uint32_t true_spread = pair.second;
            if (flowPriority_dict[src_flow] > 0.75 * gamma && true_spread >= static_cast<uint32_t>(threshold_Superspread)) {
                true_superspread_flows.insert(Superspread_P2S_FC(src_flow, flowPriority_dict[src_flow])); 
            }
        }
        set<uint32_t> true_superspread_labels;
        for (const auto& cell : true_superspread_flows) {
            true_superspread_labels.insert(cell.flowLabel);
        }
        set<Superspread_P2S_FC> estimated_superspread_flows = sketch.getSuperSpread();
        set<uint32_t> estimated_superspread_labels;
        for (const auto& cell : estimated_superspread_flows) {
            estimated_superspread_labels.insert(cell.flowLabel);
        }
        vector<uint32_t> common_superspread_labels_vec;
        set_intersection(true_superspread_labels.begin(), true_superspread_labels.end(),
                        estimated_superspread_labels.begin(), estimated_superspread_labels.end(),
                        back_inserter(common_superspread_labels_vec));
        size_t num_common_flows = common_superspread_labels_vec.size();
        size_t num_true_superspreads = true_superspread_labels.size();
        size_t num_estimated_superspreads = estimated_superspread_labels.size();
        double PR = 0.0;
        if (num_estimated_superspreads > 0) {
            PR = static_cast<double>(num_common_flows) / num_estimated_superspreads;
        }
        double RR = 0.0;
        if (num_true_superspreads > 0) {
            RR = static_cast<double>(num_common_flows) / num_true_superspreads;
        }
        double F1 = 0.0;
        if (PR + RR > 0) {
            F1 = 2 * (PR * RR) / (PR + RR);
        }
        double AAE_sum = 0.0;
        double ARE_sum = 0.0;
        for (uint32_t flow_label : common_superspread_labels_vec) {
            uint32_t true_spread = real_dict[flow_label];
            uint32_t estimated_spread = est_dict[flow_label];
            AAE_sum += abs(static_cast<double>(estimated_spread) - true_spread);
            if (true_spread > 0) { 
                ARE_sum += abs(static_cast<double>(estimated_spread) - true_spread) / true_spread;
            }
        }
        double AAE = 0.0;
        double ARE = 0.0;
        if (num_common_flows > 0) {
            AAE = AAE_sum / num_common_flows;
            ARE = ARE_sum / num_common_flows;
        }
        results[5][idx] = PR;
        results[6][idx] = RR;  
        results[7][idx] = F1;  
        results[8][idx] = AAE; 
        results[9][idx] = ARE;
        cout << "Super Spread PR : " << fixed << setprecision(6) << PR << endl;
        cout << "Super Spread RR : " << fixed << setprecision(6) << RR << endl;
        cout << "Super Spread F1-score : " << fixed << setprecision(6) << F1 << endl;
        cout << "Super Spread AAE : " << fixed << setprecision(6) << AAE << endl;
        cout << "Super Spread ARE : " << fixed << setprecision(6) << ARE << endl;
        cout << endl;
    }
};
vector<double> Exp_P2S_FC(double alpha, int gamma, double wholeMemory, string priority_file_name, int threshold_Superspread) {
    cout.precision(10);
    vector<vector<double>> results(10, vector<double>(1)); 
    int threeWholeMem = wholeMemory;
    double finePartMem = 0.2 * threeWholeMem;
    double wholeMem = 0.8 * threeWholeMem; 
    double bitmapMemRatio = 0.5; 
    int bits = static_cast<int>(wholeMem * bitmapMemRatio * 1024 * 8); 
    int sizeOfBkt = 4;
    int numOfBkt = static_cast<int>((finePartMem * 1024) / 9 / sizeOfBkt); 
    int d = 3;
    int w = static_cast<int>(wholeMem * (1 - bitmapMemRatio) * 1024 / (2 * d)); 
    double base = 100; 
    double step = 2; 
    Experiment_P2S_FC exp(bits, numOfBkt, sizeOfBkt, d, w, gamma, base, step, threshold_Superspread);
    stringstream ss_alpha;
    ss_alpha << fixed << setprecision(6) << alpha;
    string alpha_str = ss_alpha.str();
    string filename = priority_file_name + "_" + to_string(alpha) + "_" + to_string(gamma) + "_.txt"; 
    exp.start(filename, results, 0);
    exp.estimate();
    exp.performance(results, 0);
    exp.pefmance2();
    exp.performance_Superspread(results, 0);
    cout << " -------------------- " << endl;
    cout << " -------------------- " << endl;
    vector<double> ret(10);
    for (int i = 0; i < 10; ++i) {
        ret[i] = results[i][0];
    }
    return ret;
}
