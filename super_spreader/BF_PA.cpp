#include <bits/stdc++.h>
#include "MurmurHash3.h" 
#include "Tool.cpp"        
using namespace std;
class CellOfBkt_BF_PA {
public:
    uint32_t flowLabel;
    int est;
    int priority;
    CellOfBkt_BF_PA(uint32_t flowLabel = 0, int est = 0, int priority = 0)
        : flowLabel(flowLabel), est(est), priority(priority) {}
    bool operator>(const CellOfBkt_BF_PA& other) const {
        return priority > other.priority;
    }
    bool operator<(const CellOfBkt_BF_PA& other) const {
        return priority < other.priority;
    }
};
class Bkt_BF_PA {
public:
    int sizeOfBkt;
    vector<CellOfBkt_BF_PA> cellArray;
    Bkt_BF_PA(int sizeOfBkt) : sizeOfBkt(sizeOfBkt) {}
    void sortBkt() {
        sort(cellArray.begin(), cellArray.end());
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
struct Superspread_BF_PA {
    uint32_t flowLabel;
    int priority;
    Superspread_BF_PA(uint32_t f = 0, int p = 0) : flowLabel(f), priority(p) {}
    bool operator<(const Superspread_BF_PA& other) const {
        return flowLabel > other.flowLabel; 
    }
};
class BF_PA {
public:
    int gamma;
    int bitmapLen;
    int k;
    vector<bool> bf;
    vector<uint32_t> bf_seeds; 
    int sizeOfBkt;
    vector<CellOfBkt_BF_PA> cellArray;
    vector<uint32_t> seeds4ChoosingCell; 
    int d;
    int w;
    vector<vector<int>> CM;
    vector<uint32_t> hash_CM_seeds; 
    set<Superspread_BF_PA> report_Superspread;
    int threshold_Superspread;
    BF_PA(int bitmapLen, int k, int sizeOfBkt, int d, int w, int maxPriority, int gamma, int threshold_Superspread)
        : gamma(gamma), bitmapLen(bitmapLen), k(k), sizeOfBkt(sizeOfBkt), d(d), w(w), threshold_Superspread(threshold_Superspread) {
        bf.resize(bitmapLen, 0);
        bf_seeds = generateSeeds(k); 
        cellArray.reserve(sizeOfBkt);
        for (int i = 0; i < sizeOfBkt; ++i) {
            cellArray.emplace_back(0, 0, 0); 
        }
        seeds4ChoosingCell = generateSeeds(gamma); 
        CM.resize(d, vector<int>(w, 0));
        hash_CM_seeds = generateSeeds(d); 
    }
    void update(uint32_t src, uint32_t dst, int priority) {
        bool flag = false;
        uint32_t combined_hash_input = src ^ dst; 
        for (int i = 0; i < k; ++i) {
            uint32_t hash_value;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&combined_hash_input), sizeof(combined_hash_input), bf_seeds[i], &hash_value);
            uint32_t hash_idx = hash_value % bitmapLen;
            if (bf[hash_idx] == 0) {
                flag = true;
                bf[hash_idx] = 1;
            }
        }
        if (!flag) {
            return;
        }
        int estIncrement = 1;
        int minPriority = 99999;
        int minIdx = -1;
        for (int i = 0; i < priority; ++i) {
            uint32_t hash_value_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seeds4ChoosingCell[i], &hash_value_cell);
            uint32_t idx4Cell = hash_value_cell % sizeOfBkt;
            if (cellArray[idx4Cell].flowLabel == src) {
                cellArray[idx4Cell].est += estIncrement;
                if (priority > 0.75 * gamma && cellArray[idx4Cell].est >= threshold_Superspread) {
                    report_Superspread.insert(Superspread_BF_PA(src, priority));
                }
                return;
            }
            if (cellArray[idx4Cell].priority < minPriority) {
                minPriority = cellArray[idx4Cell].priority;
                minIdx = idx4Cell;
            }
        }
        if (priority > minPriority) {
            coarsePartInsert(cellArray[minIdx].flowLabel, cellArray[minIdx].est);
            cellArray[minIdx] = CellOfBkt_BF_PA(src, estIncrement, priority);
            if (priority > 0.75 * gamma && estIncrement >= threshold_Superspread) {
                report_Superspread.insert(Superspread_BF_PA(src, priority));
            }
        } else {
            coarsePartInsert(src, estIncrement);
        }
    }
    int estimate(uint32_t src, int priority) {
        for (int i = 0; i < priority; ++i) {
            uint32_t hash_value_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seeds4ChoosingCell[i], &hash_value_cell);
            uint32_t idx4Cell = hash_value_cell % sizeOfBkt;
            if (cellArray[idx4Cell].flowLabel == src) {
                return cellArray[idx4Cell].est;
            }
        }
        return coarsePartQuery(src);
    }
    void coarsePartInsert(uint32_t src, int estIncrement) {
        for (int i = 0; i < d; ++i) {
            uint32_t hash_value_idx;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), hash_CM_seeds[i], &hash_value_idx);
            uint32_t hash_idx = hash_value_idx % w;
            CM[i][hash_idx] += estIncrement;
        }
    }
    int coarsePartQuery(uint32_t src) {
        int min_val = numeric_limits<int>::max(); 
        for (int i = 0; i < d; ++i) {
            uint32_t hash_val_cm;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), hash_CM_seeds[i], &hash_val_cm);
            uint32_t hash_idx = hash_val_cm % w;
            min_val = min(min_val, CM[i][hash_idx]); 
        }
        return min_val;
    }
    set<Superspread_BF_PA> getSuperSpread() const {
        return report_Superspread;
    }
};
class Experiment_BF_PA {
public:
    int gamma;
    int threshold_Superspread;
    BF_PA sketch;
    unordered_map<uint32_t, set<uint32_t>> real_set_dict;
    unordered_map<uint32_t, int> real_dict;
    unordered_map<uint32_t, int> est_dict;
    unordered_map<uint32_t, int> flowPriority_dict;
    Experiment_BF_PA(int bitmapLen, int k, int sizeOfBkt, int d, int w, int maxPriority, int gamma, int threshold_Superspread)
        : gamma(gamma), threshold_Superspread(threshold_Superspread), sketch(bitmapLen, k, sizeOfBkt, d, w, maxPriority, gamma, threshold_Superspread) {}
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
            string dst_str = line.substr(first_space + 1, second_space - (first_space + 1));
            string priorityStr = line.substr(second_space + 1);
            uint32_t src_flow_label = stoul(src_str);
            uint32_t dst_flow_label = stoul(dst_str);
            int priority = stoi(priorityStr);
            parsed_data.emplace_back(src_flow_label, dst_flow_label, priority);
        }
        file.close();
        auto start_time = chrono::high_resolution_clock::now();
        for (const auto& [src, dst, priority] : parsed_data) {
            sketch.update(src, dst, priority);
        }
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        double throughput = static_cast<double>(parsed_data.size()) / elapsed_time.count() / 1000000;
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
            est_dict[src] = static_cast<int>(sketch.estimate(src, flowPriority_dict[src]));
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "Success: " << progress << "%" << flush;
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
                cout << "\r[Messages] Success: " << progress << "%" << flush;
            }
        }
        cout << endl; 
        double highPriorityARE = highPriorityErrorSUM / highPriorityFlowNum;
        double lowPriorityARE = lowPriorityErrorSUM / lowPriorityFlowNum;
        double allFlowsARE = (highPriorityErrorSUM + lowPriorityErrorSUM) / (highPriorityFlowNum + lowPriorityFlowNum);
        int countHigh = 0;
        for (size_t i = 0; i < sketch.cellArray.size(); ++i) {
            if (sketch.cellArray[i].priority > gamma * 0.75) {
                countHigh++;
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
                for (int i = 0; i < flowPriority_dict[src]; ++i) {
                    uint32_t hash_value_cell;
                    MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), sketch.seeds4ChoosingCell[i], &hash_value_cell);
                    uint32_t idx4Cell = hash_value_cell % sketch.sizeOfBkt;
                    if (sketch.cellArray[idx4Cell].flowLabel == src) {
                        int estTemp = sketch.cellArray[idx4Cell].est;
                        highPriorityFlowNum++;
                        highPriorityErrorSUM += abs(static_cast<double>(real_dict[src]) - estTemp) / real_dict[src];
                        break; 
                    }
                }
            }
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages] Success: " << progress << "%" << flush;
            }
        }
    }
    void performance_Superspread(vector<vector<double>>& results, int idx) {
        set<Superspread_BF_PA> true_superspread_flows;
        cout << real_dict.size() << endl;
        cout << threshold_Superspread << endl;
        for (const auto& pair : real_dict) { 
            uint32_t src_flow = pair.first;
            uint32_t true_spread = pair.second;
            if (flowPriority_dict[src_flow] > 0.75 * gamma && true_spread >= static_cast<uint32_t>(threshold_Superspread)) {
                true_superspread_flows.insert(Superspread_BF_PA(src_flow, flowPriority_dict[src_flow])); 
            }
        }
        set<uint32_t> true_superspread_labels;
        for (const auto& cell : true_superspread_flows) {
            true_superspread_labels.insert(cell.flowLabel);
        }
        set<Superspread_BF_PA> estimated_superspread_flows = sketch.getSuperSpread();
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
vector<double> Exp_BF_PA(double alpha, int gamma, double wholeMemory, string priority_file_name, int threshold_Superspread) {
    cout.precision(10); 
    vector<vector<double>> results(10, vector<double>(1)); 
    int threeWholeMem = wholeMemory;
    double finePartMem = 0.2 * threeWholeMem;
    double wholeMem = 0.8 * threeWholeMem;
    double bitmapMemRatio = 0.5;
    int bits = static_cast<int>(wholeMem * bitmapMemRatio * 1024 * 8);
    int k = 2;
    int sizeOfBkt = static_cast<int>((finePartMem * 1024) / 9); 
    int d = 2;
    int w = static_cast<int>(wholeMem * (1 - bitmapMemRatio) * 1024 / (d * 2)); 
    int maxPriority = gamma;
    Experiment_BF_PA exp(bits, k, sizeOfBkt, d, w, maxPriority, gamma, threshold_Superspread);
    stringstream ss;
    ss << fixed << setprecision(6) << alpha;
    string alpha_str = ss.str();
    string data_filename = priority_file_name + "_" + to_string(alpha) + "_" + to_string(gamma) + "_.txt";
    exp.start(data_filename, results, 0);
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
