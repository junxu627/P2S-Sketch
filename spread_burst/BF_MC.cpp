#include <bits/stdc++.h>
#include "MurmurHash3.h" 
#include "Tool.cpp"      
using namespace std;
struct BurstRecord_BF_MC {
    uint32_t flowLabel;
    int period; 
    BurstRecord_BF_MC(uint32_t f = 0, int p = 0) : flowLabel(f), period(p) {}
    bool operator<(const BurstRecord_BF_MC& other) const {
        if (flowLabel != other.flowLabel) {
            return flowLabel < other.flowLabel;
        }
        return period < other.period;
    }
};
class CellOfBkt_BF_MC {
public:
    uint32_t flowLabel; 
    int pre_est;       
    int cur_est;       
    int priority;      
    CellOfBkt_BF_MC(uint32_t flowLabel_val = 0, int pre_est_val = 0, int cur_est_val = 0, int priority_val = 0)
        : flowLabel(flowLabel_val), pre_est(pre_est_val), cur_est(cur_est_val), priority(priority_val) {}
    bool operator>(const CellOfBkt_BF_MC& other) const {
        return priority > other.priority;
    }
    bool operator<(const CellOfBkt_BF_MC& other) const {
        return priority < other.priority;
    }
};
class Bkt_BF_MC {
public:
    int sizeOfBkt;
    vector<CellOfBkt_BF_MC> cellArray;
    Bkt_BF_MC(int sizeOfBkt_val) : sizeOfBkt(sizeOfBkt_val) {
        cellArray.reserve(sizeOfBkt); 
        for (int i = 0; i < sizeOfBkt_val; ++i) {
            cellArray.emplace_back();
        }
    }
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
class BF_MC {
public:
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
    double ratioThreshold; 
    set<BurstRecord_BF_MC> detected_burst_records; 
    int burstThreshold; 
    int gamma; 
    BF_MC(int bitmapLen_val, int k_val, int numOfLayers_val, int sizeOfBkt_val, 
            int d_val, int w_val, int gamma, double ratioThreshold, int burstThreshold)
        : bitmapLen(bitmapLen_val), k_bf(k_val), 
            numOfLayers(numOfLayers_val), sizeOfBkt_mc(sizeOfBkt_val),
            d_cm(d_val), w_cm(w_val),
            ratioThreshold(ratioThreshold), burstThreshold(burstThreshold),
            gamma(gamma) { 
        bf.resize(bitmapLen, false);
        bf_seeds = generateSeeds(k_bf); 
        for (int i = 0; i < numOfLayers; ++i) {
            bktArray.emplace_back(sizeOfBkt_mc); 
        }
        seeds4ChoosingCell = generateSeeds(numOfLayers_val); 
        CM.resize(d_cm, vector<int>(w_cm, 0));
        hash_CM_seeds = generateSeeds(d_cm); 
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
        int est_increment = 1; 
        int layerNow = 0;
        mc_update(src_uint, est_increment, priority, layerNow);
    }
    void mc_update(uint32_t f, int est_increment, int priority, int layerNow) {
        if (layerNow == numOfLayers) { 
            coarsePartInsert(f, est_increment);
            return;
        }
        for (int i = layerNow; i < numOfLayers; ++i) {
            uint32_t hash_val_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&f), sizeof(f), seeds4ChoosingCell[i], &hash_val_cell);
            uint32_t idx4Cell = hash_val_cell % sizeOfBkt_mc;
            if (bktArray[i].cellArray[idx4Cell].flowLabel == f) { 
                bktArray[i].cellArray[idx4Cell].cur_est += est_increment; 
                return;
            } else { 
                if (priority > bktArray[i].cellArray[idx4Cell].priority) { 
                    uint32_t fTemp = bktArray[i].cellArray[idx4Cell].flowLabel;
                    int curEstTemp = bktArray[i].cellArray[idx4Cell].cur_est; 
                    int pTemp = bktArray[i].cellArray[idx4Cell].priority;
                    bktArray[i].cellArray[idx4Cell] = CellOfBkt_BF_MC(f, 0, est_increment, priority); 
                    if (fTemp != 0) { 
                        mc_update(fTemp, curEstTemp, pTemp, i + 1); 
                    }
                    return;
                }
            }
        }
        coarsePartInsert(f, est_increment);
    }
    void endPeriod(int currentPeriodIdx) {
        for (int i = 0; i < numOfLayers; ++i) { 
            for (size_t j = 0; j < bktArray[i].cellArray.size(); ++j) { 
                CellOfBkt_BF_MC& cell = bktArray[i].cellArray[j];
                if (cell.flowLabel != 0 && cell.priority > 0.75 * gamma) { 
                    if (cell.cur_est > burstThreshold && static_cast<double>(cell.cur_est)  > ratioThreshold * cell.pre_est) {
                        detected_burst_records.insert(BurstRecord_BF_MC(cell.flowLabel, currentPeriodIdx));
                    }
                }
            }
        }
        for (int i = 0; i < numOfLayers; ++i) {
            for (size_t j = 0; j < bktArray[i].cellArray.size(); ++j) {
                CellOfBkt_BF_MC& cell = bktArray[i].cellArray[j];
                if (cell.flowLabel != 0) {
                    cell.pre_est = cell.cur_est;
                    cell.cur_est = 0;
                }
            }
        }
        fill(bf.begin(), bf.end(), false); 
    }
    int estimate(uint32_t src_uint, int priority) { 
        for (int i = 0; i < numOfLayers; ++i) {
            uint32_t hash_val_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src_uint), sizeof(src_uint), seeds4ChoosingCell[i], &hash_val_cell);
            uint32_t idx4Cell = hash_val_cell % sizeOfBkt_mc;
            if (bktArray[i].cellArray[idx4Cell].flowLabel == src_uint) {
                return bktArray[i].cellArray[idx4Cell].cur_est; 
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
    const set<BurstRecord_BF_MC>& getDetectedBurstRecords() const {
        return detected_burst_records;
    }
};
class Experiment_BF_MC {
public:
    int gamma_exp; 
    BF_MC sketch;
    unordered_map<uint32_t, int> flowPriority_dict; 
    set<BurstRecord_BF_MC> true_burst_records; 
    double burst_threshold_exp; 
    int burst_size_threshold_exp; 
    Experiment_BF_MC(int bitmapLen, int k, int numOfLayers, int sizeOfBkt, int d, int w, int gamma_val, double Th_val, int burst_size_val)
        : gamma_exp(gamma_val), burst_threshold_exp(Th_val), burst_size_threshold_exp(burst_size_val),
          sketch(bitmapLen, k, numOfLayers, sizeOfBkt, d, w, gamma_val, Th_val, burst_size_val) {}
    void process_data_and_detect_bursts(const string& filename, int T) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }
        vector<tuple<uint32_t, uint32_t, int>> all_parsed_data;
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
            all_parsed_data.emplace_back(src_flow_label, dst_flow_label, priority);
            flowPriority_dict[src_flow_label] = priority; 
        }
        file.close();
        long long total_packets = all_parsed_data.size();
        long long packets_per_period = (total_packets + T - 1) / T; 
        unordered_map<uint32_t, uint32_t> real_spreads_prev_period; 
        auto start_time = chrono::high_resolution_clock::now();
        for (int period_idx = 0; period_idx < T; ++period_idx) {
            unordered_map<uint32_t, set<uint32_t>> current_period_real_sets; 
            unordered_map<uint32_t, uint32_t> real_spreads_curr_period; 
            long long start_packet_idx = period_idx * packets_per_period;
            long long end_packet_idx = min((period_idx + 1) * packets_per_period, total_packets);
            for (long long i = start_packet_idx; i < end_packet_idx; ++i) {
                const auto& [src, dst, priority] = all_parsed_data[i];
                sketch.update(src, dst, priority); 
                current_period_real_sets[src].insert(dst); 
            }
            for (const auto& pair : current_period_real_sets) {
                real_spreads_curr_period[pair.first] = pair.second.size();
            }
            set<uint32_t> all_active_flows_in_window;
            for(const auto& pair : real_spreads_curr_period) all_active_flows_in_window.insert(pair.first);
            for(const auto& pair : real_spreads_prev_period) all_active_flows_in_window.insert(pair.first);
            for (uint32_t flow_label : all_active_flows_in_window) {
                if (flowPriority_dict[flow_label] <= 0.75 * gamma_exp) { 
                    continue;
                }
                uint32_t current_spread = 0;
                if (real_spreads_curr_period.count(flow_label)) {
                    current_spread = real_spreads_curr_period[flow_label];
                }
                uint32_t previous_spread = 0;
                if (real_spreads_prev_period.count(flow_label)) {
                    previous_spread = real_spreads_prev_period[flow_label];
                }
                if (current_spread > burst_size_threshold_exp && static_cast<double>(current_spread) > burst_threshold_exp * previous_spread) {
                    true_burst_records.insert(BurstRecord_BF_MC(flow_label, period_idx));
                }
            }
            sketch.endPeriod(period_idx);
            real_spreads_prev_period = real_spreads_curr_period; 
        }
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        double throughput = static_cast<double>(total_packets) / elapsed_time.count() / 1000000;
        cout << "throughput : " << throughput << " MPPS" << endl;
    }
    void performance_burst(vector<vector<double>>& results, int idx) {
        const set<BurstRecord_BF_MC>& detected_bursts_set = sketch.getDetectedBurstRecords();
        vector<BurstRecord_BF_MC> true_bursts_vec(true_burst_records.begin(), true_burst_records.end());
        vector<BurstRecord_BF_MC> detected_bursts_vec(detected_bursts_set.begin(), detected_bursts_set.end());
        vector<BurstRecord_BF_MC> common_bursts_vec; 
        set_intersection(true_bursts_vec.begin(), true_bursts_vec.end(),
                                     detected_bursts_vec.begin(), detected_bursts_vec.end(),
                                     back_inserter(common_bursts_vec));
        size_t true_positives = common_bursts_vec.size();
        size_t num_true_bursts = true_burst_records.size();
        size_t num_detected_bursts = detected_bursts_set.size(); 
        double PR = 0.0;
        if (num_detected_bursts > 0) {
            PR = static_cast<double>(true_positives) / num_detected_bursts;
        }
        double RR = 0.0;
        if (num_true_bursts > 0) {
            RR = static_cast<double>(true_positives) / num_true_bursts;
        }
        double F1 = 0.0;
        if (PR + RR > 0) {
            F1 = (2 * PR * RR) / (PR + RR);
        }
        results[0][idx] = PR;
        results[1][idx] = RR;
        results[2][idx] = F1;
        cout << "True Burst Number: " << true_bursts_vec.size() << endl;
        cout << "Report Burst Number: " << detected_bursts_vec.size() << endl;
        cout << "Burst PR : " << fixed << setprecision(6) << PR << endl;
        cout << "Burst RR : " << fixed << setprecision(6) << RR << endl;
        cout << "Burst F1-score : " << fixed << setprecision(6) << F1 << endl;
        cout << endl;
    }
};
vector<double> Exp_BF_MC(double alpha, int gamma, double wholeMemory, string priority_file_name, double burstTh, int T, int burstThreshold, string filename) {
    cout.precision(10); 
    vector<vector<double>> results(5, vector<double>(1)); 
    int threeWholeMem = wholeMemory; 
    double finePartMem = 0.2 * threeWholeMem; 
    double coarsePartMem = 0.8 * threeWholeMem; 
    double bitmapMemRatio = 0.5; 
    int bits = static_cast<int>(coarsePartMem * bitmapMemRatio * 1024 * 8); 
    int k = 2; 
    int numOfLayers = 5; 
    int sizeOfBkt = static_cast<int>(finePartMem * 1024 / 13 / numOfLayers); 
    int d = 3; 
    int w = static_cast<int>(coarsePartMem * (1 - bitmapMemRatio) * 1024 / (d * 2)); 
    int maxPriority = gamma; 
    Experiment_BF_MC exp(bits, k, numOfLayers, sizeOfBkt, d, w, gamma, burstTh, burstThreshold);
    stringstream ss;
    ss << fixed << setprecision(6) << alpha;
    string alpha_str = ss.str();
    exp.process_data_and_detect_bursts(filename, T); 
    exp.performance_burst(results, 0); 
    cout << " -------------------- " << endl;
    cout << " -------------------- " << endl;
    vector<double> ret(3); 
    for (int i = 0; i < 3; ++i) { 
        ret[i] = results[i][0];
    }
    return ret;
}
