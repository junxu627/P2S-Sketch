#include <bits/stdc++.h>
#include "MurmurHash3.h"
#include "Tool.cpp" 
using namespace std;
class CellOfBkt_BF_MC {
public:
    uint32_t flowLabel; 
    int est;            
    int priority;       
    CellOfBkt_BF_MC(uint32_t flowLabel_val, int est_val, int priority_val)
        : flowLabel(flowLabel_val), est(est_val), priority(priority_val) {}
    bool operator<(const CellOfBkt_BF_MC& other) const {
        return priority < other.priority;
    }
    bool operator>(const CellOfBkt_BF_MC& other) const {
        return priority > other.priority;
    }
};
class Bkt_BF_MC {
public:
    int sizeOfBkt;
    vector<CellOfBkt_BF_MC> cellArray;
    Bkt_BF_MC(int sizeOfBkt_val) : sizeOfBkt(sizeOfBkt_val) {}
    void sortBkt() {
        sort(cellArray.begin(), cellArray.end());
    }
    int posOfFlowInBkt(uint32_t flowLabel_to_find) {
        for (size_t idx = 0; idx < cellArray.size(); ++idx) {
            if (cellArray[idx].flowLabel == flowLabel_to_find) {
                return static_cast<int>(idx);
            }
        }
        return -1;
    }
    bool isBktFull() {
        return cellArray.size() == static_cast<size_t>(sizeOfBkt);
    }
};
struct CellOfTopk_BF_MC {
    uint32_t flowLabel;
    int est; 
    int priority;
    CellOfTopk_BF_MC(uint32_t f = 0, int e = 0, int p = 0) : flowLabel(f), est(e), priority(p) {}
    bool operator<(const CellOfTopk_BF_MC& other) const {
        return est > other.est; 
    }
};
class BF_MC {
public:
    int gamma;
    int bitmapLen;
    int k_bf; 
    vector<bool> bf; 
    vector<uint32_t> bf_seeds;
    int numOfLayers;
    int sizeOfBkt_mc; 
    vector<Bkt_BF_MC> bktArray;
    vector<uint32_t> seeds4ChoosingCell;
    int d_cm; 
    int w_cm; 
    vector<vector<int>> CM; 
    vector<uint32_t> hash_CM_seeds;
    vector<CellOfTopk_BF_MC> topkArray;
    int k_topk; 
    BF_MC(int bitmapLen_val, int k_val, int numOfLayers_val, int sizeOfBkt_val, 
            int d_val, int w_val, int gamma, int k_topk_val)
        : gamma(gamma), bitmapLen(bitmapLen_val), k_bf(k_val), 
            numOfLayers(numOfLayers_val), sizeOfBkt_mc(sizeOfBkt_val),
            d_cm(d_val), w_cm(w_val), k_topk(k_topk_val) {
        bf.resize(bitmapLen, false);
        bf_seeds = generateSeeds(k_bf); 
        for (int i = 0; i < numOfLayers; ++i) {
            bktArray.emplace_back(sizeOfBkt_mc); 
            for (int j = 0; j < sizeOfBkt_mc; ++j) {
                bktArray[i].cellArray.emplace_back(0, 0, 0); 
            }
        }
        seeds4ChoosingCell = generateSeeds(numOfLayers_val); 
        CM.resize(d_cm, vector<int>(w_cm, 0));
        hash_CM_seeds = generateSeeds(d_cm); 
    }
    void updateTopk(uint32_t flowLabel, float spread, int priority) {
        if (priority <= 0.75 * gamma) {
            return ;
        }
        int existing_idx = -1;
        for (size_t i = 0; i < topkArray.size(); ++i) {
            if (topkArray[i].flowLabel == flowLabel) {
                existing_idx = static_cast<int>(i);
                break;
            }
        }
        if (existing_idx != -1) {
            topkArray[existing_idx].est = spread;
            topkArray[existing_idx].priority = priority;
            sort(topkArray.begin(), topkArray.end());
        } else {
            if (topkArray.size() < static_cast<size_t>(k_topk)) {
                topkArray.emplace_back(flowLabel, spread, priority);
                sort(topkArray.begin(), topkArray.end()); 
            } else if (spread > topkArray.back().est) { 
                topkArray.back() = CellOfTopk_BF_MC(flowLabel, spread, priority); 
                sort(topkArray.begin(), topkArray.end()); 
            }
        }
    }
    void deleteFromTopk(uint32_t flowLabel) {
        for (auto it = topkArray.begin(); it != topkArray.end(); ++it) {
            if (it->flowLabel == flowLabel) {
                topkArray.erase(it); 
                sort(topkArray.begin(), topkArray.end());
                return; 
            }
        }
    }
    void update(uint32_t src_uint, uint32_t dst_uint, int priority) {
        bool flag = false;
        uint32_t combined_key_bf = src_uint ^ dst_uint; 
        for (int i = 0; i < k_bf; ++i) {
            uint32_t hash_val_bf;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&combined_key_bf), sizeof(combined_key_bf), bf_seeds[i], &hash_val_bf);
            uint32_t hash_idx = hash_val_bf % bitmapLen;
            if (!bf[hash_idx]) {
                flag = true;
                bf[hash_idx] = true;
            }
        }
        if (!flag) {
            return;
        }
        int est = 1; 
        int layerNow = 0;
        mc_update(src_uint, est, priority, layerNow);
    }
    void mc_update(uint32_t f, int est, int priority, int layerNow) {
        if (layerNow == numOfLayers) { 
            coarsePartInsert(f, est);
            deleteFromTopk(f);
            return;
        }
        for (int i = layerNow; i < numOfLayers; ++i) {
            uint32_t hash_val_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&f), sizeof(f), seeds4ChoosingCell[i], &hash_val_cell);
            uint32_t idx4Cell = hash_val_cell % sizeOfBkt_mc;
            if (bktArray[i].cellArray[idx4Cell].flowLabel == f) { 
                bktArray[i].cellArray[idx4Cell].est += est;
                updateTopk(f, bktArray[i].cellArray[idx4Cell].est, priority);
                return;
            } else { 
                if (priority > bktArray[i].cellArray[idx4Cell].priority) { 
                    uint32_t fTemp = bktArray[i].cellArray[idx4Cell].flowLabel;
                    int estTemp = bktArray[i].cellArray[idx4Cell].est;
                    int pTemp = bktArray[i].cellArray[idx4Cell].priority;
                    bktArray[i].cellArray[idx4Cell] = CellOfBkt_BF_MC(f, est, priority);
                    updateTopk(f, est, priority);
                    if (fTemp != 0) { 
                            mc_update(fTemp, estTemp, pTemp, i + 1); 
                    }
                    return;
                }
            }
        }
        coarsePartInsert(f, est);
    }
    int estimate(uint32_t src_uint) {
        for (int i = 0; i < numOfLayers; ++i) {
            uint32_t hash_val_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src_uint), sizeof(src_uint), seeds4ChoosingCell[i], &hash_val_cell);
            uint32_t idx4Cell = hash_val_cell % sizeOfBkt_mc;
            if (bktArray[i].cellArray[idx4Cell].flowLabel == src_uint) {
                return bktArray[i].cellArray[idx4Cell].est; 
            }
        }
        return coarsePartQuery(src_uint); 
    }
    void coarsePartInsert(uint32_t src_uint, int estIncrement) {
        for (int i = 0; i < d_cm; ++i) {
            uint32_t hash_val_cm;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src_uint), sizeof(src_uint), hash_CM_seeds[i], &hash_val_cm);
            uint32_t hash_idx = hash_val_cm % w_cm;
            CM[i][hash_idx] += estIncrement; 
        }
    }
    int coarsePartQuery(uint32_t src_uint) {
        int min_val = numeric_limits<int>::max(); 
        for (int i = 0; i < d_cm; ++i) {
            uint32_t hash_val_cm;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src_uint), sizeof(src_uint), hash_CM_seeds[i], &hash_val_cm);
            uint32_t hash_idx = hash_val_cm % w_cm;
            min_val = min(min_val, CM[i][hash_idx]); 
        }
        return min_val;
    }
    vector<CellOfTopk_BF_MC> getTopkSuperSpread() const {
        return topkArray;
    }
};
class Experiment_BF_MC {
public:
    int gamma;
    BF_MC sketch;
    unordered_map<uint32_t, set<uint32_t>> real_set_dict;
    unordered_map<uint32_t, int> real_dict;
    unordered_map<uint32_t, int> est_dict;
    unordered_map<uint32_t, int> flowPriority_dict;
    Experiment_BF_MC(int bitmapLen, int k, int numOfLayers, int sizeOfBkt, int d, int w, int gamma, int k_topk)
        : gamma(gamma), 
            sketch(bitmapLen, k, numOfLayers, sizeOfBkt, d, w, gamma, k_topk) {}
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
    void estimate_flows() { 
        int total_flows = real_set_dict.size();
        if (total_flows == 0) {
            return;
        }
        int count = 0;
        int update_interval = total_flows / 100 > 0 ? total_flows / 100 : 1; 
        for (auto const& [src, dst_set] : real_set_dict) { 
            real_dict[src] = dst_set.size(); 
            est_dict[src] = sketch.estimate(src); 
            count++;
            if (count % update_interval == 0 || count == total_flows) {
                int progress = static_cast<int>((static_cast<double>(count) / total_flows) * 100);
                cout << "\r[Messages] Success:" << progress << "%" << flush; 
            }
        }
    }
    void performance(vector<vector<double>>& results, int idx) { 
        if (real_dict.empty()) {
            return;
        }
        double highPriorityErrorSUM = 0;
        int highPriorityFlowNum = 0;
        double lowPriorityErrorSUM = 0;
        int lowPriorityFlowNum = 0;
        for (auto const& [src, real_val] : real_dict) {
            if (flowPriority_dict.count(src)) { 
                    if (flowPriority_dict[src] > gamma * 0.75) { 
                    highPriorityFlowNum++;
                    if (real_val > 0) highPriorityErrorSUM += abs(static_cast<double>(real_val) - est_dict[src]) / real_val; 
                } else {
                    lowPriorityFlowNum++;
                    if (real_val > 0) lowPriorityErrorSUM += abs(static_cast<double>(real_val) - est_dict[src]) / real_val; 
                }
            }
        }
        double highPriorityARE = (highPriorityFlowNum > 0) ? highPriorityErrorSUM / highPriorityFlowNum : 0.0; 
        double lowPriorityARE = (lowPriorityFlowNum > 0) ? lowPriorityErrorSUM / lowPriorityFlowNum : 0.0; 
        double allFlowsARE = ((highPriorityFlowNum + lowPriorityFlowNum) > 0) ? 
                                (highPriorityErrorSUM + lowPriorityErrorSUM) / (highPriorityFlowNum + lowPriorityFlowNum) : 0.0; 
        int countHigh = 0;
        for (const auto& layer : sketch.bktArray) { 
            for (const auto& cell : layer.cellArray) {
                if (cell.priority > gamma * 0.75) { 
                    countHigh++;
                }
            }
        }
        double pr = 1.0;  
        double rr = (highPriorityFlowNum > 0) ? static_cast<double>(countHigh) / highPriorityFlowNum : 0.0; 
        double f1_score = (pr + rr > 1e-9) ? (2.0 * pr * rr / (pr + rr)) : 0.0; 
        cout << fixed << setprecision(6);
        cout << "f1-score : " << f1_score << endl; 
        cout << "highPriorityARE : " << highPriorityARE << endl; 
        cout << "lowPriorityARE : " << lowPriorityARE << endl; 
        cout << "allFlowsARE : " << allFlowsARE << endl; 
        results[0][idx] = highPriorityARE;
        results[1][idx] = lowPriorityARE;
        results[2][idx] = allFlowsARE;
        results[3][idx] = f1_score;
    }
    void pefmance2() { 
        if (real_dict.empty()) {
            return;
        }
        double highPriorityErrorSUM = 0;
        int highPriorityFlowNum = 0;
        for (auto const& [src, real_val] : real_dict) {
                if (flowPriority_dict.count(src) && flowPriority_dict[src] > gamma * 0.75) { 
                bool found_in_bkt = false;
                int estTemp = 0;
                for (const auto& layer : sketch.bktArray) {
                    for (const auto& cell : layer.cellArray) {
                        if (cell.flowLabel == src) {
                            estTemp = cell.est;
                            found_in_bkt = true;
                            break;
                        }
                    }
                    if (found_in_bkt) break;
                }
                if (found_in_bkt) {
                    highPriorityFlowNum++;
                    if (real_val > 0) highPriorityErrorSUM += abs(static_cast<double>(real_val) - estTemp) / real_val; 
                }
            }
        }
    }
    void new_performance_topk(vector<vector<double>>& results, int idx) {
        vector<CellOfTopk_BF_MC> true_high_priority_flows;
        for (const auto& pair : real_dict) { 
            uint32_t src_flow = pair.first;
            uint32_t true_spread = pair.second;
            int priority = flowPriority_dict[src_flow];
            if (priority > gamma * 0.75) { 
                true_high_priority_flows.emplace_back(src_flow, static_cast<int>(true_spread), priority);
            }
        }
        sort(true_high_priority_flows.begin(), true_high_priority_flows.end(),
                  [](const CellOfTopk_BF_MC& a, const CellOfTopk_BF_MC& b) {
                      return a.est > b.est; 
                  });
        vector<CellOfTopk_BF_MC> true_topk_flows_extended;
        if (!true_high_priority_flows.empty()) {
            if (static_cast<size_t>(sketch.k_topk) <= true_high_priority_flows.size()) {
                int k_th_spread_value = true_high_priority_flows[sketch.k_topk - 1].est;
                for (const auto& cell : true_high_priority_flows) {
                    if (cell.est >= k_th_spread_value) {
                        true_topk_flows_extended.push_back(cell);
                    } else {
                        break;
                    }
                }
            } else {
                true_topk_flows_extended = true_high_priority_flows;
            }
        }
        set<uint32_t> true_topk_labels_extended;
        for (const auto& cell : true_topk_flows_extended) {
            true_topk_labels_extended.insert(cell.flowLabel);
        }
        vector<CellOfTopk_BF_MC> estimated_topk_flows = sketch.getTopkSuperSpread();
        if (estimated_topk_flows.size() > static_cast<size_t>(sketch.k_topk)) {
            estimated_topk_flows.resize(sketch.k_topk); 
        }
        set<uint32_t> estimated_topk_labels;
        for (const auto& cell : estimated_topk_flows) {
            estimated_topk_labels.insert(cell.flowLabel);
        }
        vector<uint32_t> common_topk_labels_vec;
        set_intersection(true_topk_labels_extended.begin(), true_topk_labels_extended.end(),
                              estimated_topk_labels.begin(), estimated_topk_labels.end(),
                              back_inserter(common_topk_labels_vec));
        size_t num_common_flows = common_topk_labels_vec.size();
        double PR = 0.0;
        if (sketch.k_topk > 0) { 
            PR = static_cast<double>(num_common_flows) / sketch.k_topk;
        }
        double AAE_sum = 0.0;
        double ARE_sum = 0.0;
        for (uint32_t flow_label : common_topk_labels_vec) {
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
        results[6][idx] = AAE;
        results[7][idx] = ARE;
        cout << "Top-K PR : " << fixed << setprecision(6) << PR << endl;
        cout << "Top-K AAE : " << fixed << setprecision(6) << AAE << endl;
        cout << "Top-K ARE : " << fixed << setprecision(6) << ARE << endl;
        cout << endl;
    }
};
vector<double> Exp_BF_MC(double alpha, int gamma, double wholeMemory, string priority_file_name, int k_topk) {
    cout.precision(10); 
    vector<vector<double>> results(8, vector<double>(1)); 
    double threeWholeMem_KB = wholeMemory;
    double finePartMemRatio = 0.2;
    double coarsePartMemRatio = 0.8; 
    double finePartMem_KB = finePartMemRatio * threeWholeMem_KB;
    double coarsePartMem_KB = coarsePartMemRatio * threeWholeMem_KB;
    double bitmapMemRatio_of_coarse = 0.5;
    int bits_main = static_cast<int>(coarsePartMem_KB * bitmapMemRatio_of_coarse * 1024 * 8);
    int k_main = 2;
    int numOfLayers_main = 5;
    int bytes_per_cell_assumption = 9;
    int sizeOfBkt_main = static_cast<int>((finePartMem_KB * 1024  - k_topk * 9) / bytes_per_cell_assumption / numOfLayers_main);
    int d_main = 3;
    int bytes_per_cm_element_assumption = 2;
    int w_main = static_cast<int>(coarsePartMem_KB * (1.0 - bitmapMemRatio_of_coarse) * 1024 / (bytes_per_cm_element_assumption * d_main));
    Experiment_BF_MC exp(bits_main, k_main, numOfLayers_main, sizeOfBkt_main, d_main, w_main, gamma, k_topk);
    string data_filename = priority_file_name + "_" + to_string(alpha) + "_" + to_string(gamma) + "_.txt";
    exp.start(data_filename, results, 0);
    exp.estimate_flows();
    exp.performance(results, 0);
    exp.pefmance2();
    exp.new_performance_topk(results, 0);
    cout << " -------------------- " << endl;
    cout << " -------------------- " << endl;
    vector<double> ret(8);
    for (int i = 0; i < 8; ++i) {
        ret[i] = results[i][0];
    }
    return ret;
}
