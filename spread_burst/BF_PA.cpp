#include<bits/stdc++.h>
#include "MurmurHash3.h" 
#include "Tool.cpp"        
using namespace std;
struct BurstRecord_BF_PA {
    uint32_t flowLabel;
    int period; 
    BurstRecord_BF_PA(uint32_t f = 0, int p = 0) : flowLabel(f), period(p) {}
    bool operator<(const BurstRecord_BF_PA& other) const {
        if (flowLabel != other.flowLabel) {
            return flowLabel < other.flowLabel;
        }
        return period < other.period;
    }
};
class CellOfBkt_BF_PA {
public:
    uint32_t flowLabel;
    int pre_est; 
    int cur_est; 
    int priority;
    CellOfBkt_BF_PA(uint32_t f = 0, int p_est = 0, int c_est = 0, int p = 0)
        : flowLabel(f), pre_est(p_est), cur_est(c_est), priority(p) {}
    bool operator>(const CellOfBkt_BF_PA& other) const {
        return priority > other.priority;
    }
    bool operator<(const CellOfBkt_BF_PA& other) const {
        return priority < other.priority;
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
    double ratioThreshold; 
    set<BurstRecord_BF_PA> detected_burst_records; 
    int burstThreshold;
    BF_PA(int bitmapLen, int k, int sizeOfBkt, int d, int w, int maxPriority, int gamma, double ratioThreshold, int burstThreshold)
        : gamma(gamma), bitmapLen(bitmapLen), k(k), sizeOfBkt(sizeOfBkt), d(d), w(w),
          ratioThreshold(ratioThreshold), burstThreshold(burstThreshold) {
        bf.resize(bitmapLen, false); 
        bf_seeds = generateSeeds(k); 
        cellArray.reserve(sizeOfBkt);
        for (int i = 0; i < sizeOfBkt; ++i) {
            cellArray.emplace_back(); 
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
            if (!bf[hash_idx]) { 
                flag = true;
                bf[hash_idx] = true; 
            }
        }
        if (!flag) {
            return;
        }
        int estIncrement = 1; 
        int minPriority = numeric_limits<int>::max(); 
        int minIdx = -1;
        for (int i = 0; i < priority; ++i) {
            uint32_t hash_value_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seeds4ChoosingCell[i], &hash_value_cell);
            uint32_t idx4Cell = hash_value_cell % sizeOfBkt;
            if (cellArray[idx4Cell].flowLabel == src) {
                cellArray[idx4Cell].cur_est += estIncrement; 
                return; 
            }
            if (cellArray[idx4Cell].flowLabel == 0 || cellArray[idx4Cell].priority < minPriority) {
                minPriority = cellArray[idx4Cell].priority;
                minIdx = idx4Cell;
            }
        }
        if (minIdx == -1) { 
            coarsePartInsert(src, estIncrement);
            return;
        }
        if (cellArray[minIdx].flowLabel == 0 || priority > cellArray[minIdx].priority) {
            if (cellArray[minIdx].flowLabel != 0) { 
                coarsePartInsert(cellArray[minIdx].flowLabel, cellArray[minIdx].cur_est); 
            }
            cellArray[minIdx] = CellOfBkt_BF_PA(src, 0, estIncrement, priority); 
        } else { 
            coarsePartInsert(src, estIncrement);
        }
    }
    void endPeriod(int currentPeriodIdx) {
        for (size_t j = 0; j < cellArray.size(); ++j) {
            CellOfBkt_BF_PA& cell = cellArray[j];
            if (cell.flowLabel != 0 && cell.priority > 0.75 * gamma) { 
                if(cell.cur_est > burstThreshold && static_cast<double>(cell.cur_est) > ratioThreshold * cell.pre_est) {
                    detected_burst_records.insert(BurstRecord_BF_PA(cell.flowLabel, currentPeriodIdx));
                }
            }
        }
        for (size_t j = 0; j < cellArray.size(); ++j) {
            CellOfBkt_BF_PA& cell = cellArray[j];
            if (cell.flowLabel != 0) { 
                cell.pre_est = cell.cur_est; 
                cell.cur_est = 0;           
            }
        }
        fill(bf.begin(), bf.end(), false); 
    }
    int estimate(uint32_t src, int priority) {
        for (int i = 0; i < priority; ++i) {
            uint32_t hash_value_cell;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seeds4ChoosingCell[i], &hash_value_cell);
            uint32_t idx4Cell = hash_value_cell % sizeOfBkt;
            if (cellArray[idx4Cell].flowLabel == src) {
                return cellArray[idx4Cell].cur_est; 
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
    const set<BurstRecord_BF_PA>& getDetectedBurstRecords() const {
        return detected_burst_records;
    }
};
class Experiment_BF_PA {
public:
    int gamma;
    BF_PA sketch;
    int burstThreshold;
    unordered_map<uint32_t, int> flowPriority_dict; 
    set<BurstRecord_BF_PA> true_burst_records; 
    Experiment_BF_PA(int bitmapLen, int k, int sizeOfBkt, int d, int w, int maxPriority, int gamma, double Th, int burstThreshold)
        : gamma(gamma), burstThreshold(burstThreshold), sketch(bitmapLen, k, sizeOfBkt, d, w, maxPriority, gamma, Th, burstThreshold) {}
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
                if (flowPriority_dict[flow_label] <= 0.75 * gamma) {
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
                if (current_spread > burstThreshold && static_cast<double>(current_spread) > previous_spread * sketch.ratioThreshold) {
                    true_burst_records.insert(BurstRecord_BF_PA(flow_label, period_idx));
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
        const set<BurstRecord_BF_PA>& detected_bursts = sketch.getDetectedBurstRecords();
        vector<BurstRecord_BF_PA> true_bursts_vec(true_burst_records.begin(), true_burst_records.end());
        vector<BurstRecord_BF_PA> detected_bursts_vec(detected_bursts.begin(), detected_bursts.end());
        vector<BurstRecord_BF_PA> common_bursts_vec;
        set_intersection(true_bursts_vec.begin(), true_bursts_vec.end(),
                              detected_bursts_vec.begin(), detected_bursts_vec.end(),
                              back_inserter(common_bursts_vec));
        size_t true_positives = common_bursts_vec.size();
        size_t num_true_bursts = true_burst_records.size();
        size_t num_detected_bursts = detected_bursts.size();
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
vector<double> Exp_BF_PA(double alpha, int gamma, double wholeMemory, string priority_file_name, double burstTh, int T, int burstThreshold, string filename) {
    cout.precision(10); 
    vector<vector<double>> results(5, vector<double>(1)); 
    int threeWholeMem = wholeMemory;
    double finePartMem = 0.2 * threeWholeMem;
    double wholeMem = 0.8 * threeWholeMem;
    double bitmapMemRatio = 0.5;
    int bits = static_cast<int>(wholeMem * bitmapMemRatio * 1024 * 8);
    int k = 2;
    int sizeOfBkt = static_cast<int>(finePartMem * 1024 / 13); 
    int d = 2;
    int w = static_cast<int>(wholeMem * (1 - bitmapMemRatio) * 1024 / (d * 2)); 
    int maxPriority = gamma;
    Experiment_BF_PA exp(bits, k, sizeOfBkt, d, w, maxPriority, gamma, burstTh, burstThreshold);
    stringstream ss;
    ss << fixed << setprecision(6) << alpha;
    string alpha_str = ss.str();
    exp.process_data_and_detect_bursts(filename, T);
    exp.performance_burst(results, 0);
    cout << " -------------------- " << endl;
    cout << " -------------------- " << endl;
    vector<double> ret(5);
    for (int i = 0; i < 5; ++i) {
        ret[i] = results[i][0];
    }
    return ret;
}
