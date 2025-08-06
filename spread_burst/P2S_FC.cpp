#include<bits/stdc++.h>
#include "MurmurHash3.h" 
#include "Tool.cpp"        
using namespace std;
struct BurstRecord_P2S_FC {
    uint32_t flowLabel;
    int period; 
    BurstRecord_P2S_FC(uint32_t f = 0, int p = 0) : flowLabel(f), period(p) {}
    bool operator<(const BurstRecord_P2S_FC& other) const {
        if (flowLabel != other.flowLabel) {
            return flowLabel < other.flowLabel;
        }
        return period < other.period;
    }
};
class CellOfBkt_P2S_FC {
public:
    uint32_t flowLabel;
    int pre_est; 
    int cur_est; 
    int priority;
    CellOfBkt_P2S_FC(uint32_t f = 0, int p_est = 0, int c_est = 0, int p = 0)
        : flowLabel(f), pre_est(p_est), cur_est(c_est), priority(p) {}
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
    double ratioThreshold; 
    set<BurstRecord_P2S_FC> detected_burst_records; 
    int burstThreshold;
    P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step, double ratioThreshold, int burstThreshold)
        : gamma(gamma), bitmapLen(bitmapLen), numOfBkt(numOfBkt), sizeOfBkt(sizeOfBkt), d(d), w(w), ratioThreshold(ratioThreshold), burstThreshold(burstThreshold) {
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
            zeroNumOfBitmap--;
            flag = true;
            bitmap[hash_idx] = true;
        }
        if (!flag) {
            return ;
        }
        int estIncrement = static_cast<int>(static_cast<double>(bitmapLen) / zeroNumOfBitmap);
        if (priority > gamma * 0.6) { 
            uint32_t hash_value_bkt;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seed4ChoosingBkt, &hash_value_bkt);
            uint32_t idx4Bkt = hash_value_bkt % numOfBkt;
            int posInTheBkt = bktArray[idx4Bkt].posOfFlowInBkt(src);
            if (posInTheBkt > -1) {
                bktArray[idx4Bkt].cellArray[posInTheBkt].cur_est += estIncrement;
            } else {
                if (bktArray[idx4Bkt].cellArray.size() == static_cast<size_t>(sizeOfBkt)) {
                    int min_idx = bktArray[idx4Bkt].findMinPriorityInBkt();
                    if (priority > bktArray[idx4Bkt].cellArray[min_idx].priority) {
                        CellOfBkt_P2S_FC tempCell(src, 0, estIncrement, priority);
                        coarsePartInsert(bktArray[idx4Bkt].cellArray[min_idx].flowLabel, bktArray[idx4Bkt].cellArray[min_idx].cur_est);
                        bktArray[idx4Bkt].cellArray[min_idx] = tempCell;
                    } else {
                        coarsePartInsert(src, estIncrement);
                    }
                } else {
                    CellOfBkt_P2S_FC tempCell(src, 0, estIncrement, priority);
                    bktArray[idx4Bkt].cellArray.push_back(tempCell);
                }
            }
        } else { 
            coarsePartInsert(src, estIncrement);
        }
    }
    void endPeriod(int currentPeriodIdx) {
        for (size_t i = 0; i < numOfBkt; ++i) {
            for (size_t j = 0; j < bktArray[i].cellArray.size(); ++j) {
                CellOfBkt_P2S_FC& cell = bktArray[i].cellArray[j];
                if (cell.priority > 0.75 * gamma) {
                    if (cell.cur_est > burstThreshold && static_cast<double>(cell.cur_est)  > ratioThreshold * cell.pre_est) {
                        detected_burst_records.insert(BurstRecord_P2S_FC(cell.flowLabel, currentPeriodIdx));
                    }
                }
            }
        }
        for (size_t i = 0; i < numOfBkt; ++i) {
            for (size_t j = 0; j < bktArray[i].cellArray.size(); ++j) {
                CellOfBkt_P2S_FC& cell = bktArray[i].cellArray[j];
                cell.pre_est = cell.cur_est; 
                cell.cur_est = 0;           
            }
        }
        fill(bitmap.begin(), bitmap.end(), false);
        zeroNumOfBitmap = bitmapLen;
    }
    int estimate(uint32_t src, int priority) {
        uint32_t hash_value_bkt;
        MurmurHash3_x86_32(reinterpret_cast<const void*>(&src), sizeof(src), seed4ChoosingBkt, &hash_value_bkt);
        uint32_t idx4Bkt = hash_value_bkt % numOfBkt;
        int posInTheBkt = bktArray[idx4Bkt].posOfFlowInBkt(src);
        int tmpEst = 0;
        if (posInTheBkt > -1) {
            tmpEst = bktArray[idx4Bkt].cellArray[posInTheBkt].cur_est;
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
    const set<BurstRecord_P2S_FC>& getDetectedBurstRecords() const {
        return detected_burst_records;
    }
};
class Experiment_P2S_FC {
public:
    int gamma;
    P2S_FC sketch; 
    int burstThreshold;
    unordered_map<uint32_t, int> flowPriority_dict; 
    set<BurstRecord_P2S_FC> true_burst_records; 
    Experiment_P2S_FC(int bitmapLen, int numOfBkt, int sizeOfBkt, int d, int w, int gamma, double base, double step, double Th, int burstThreshold)
        : gamma(gamma), burstThreshold(burstThreshold), sketch(bitmapLen, numOfBkt, sizeOfBkt, d, w, gamma, base, step, Th, burstThreshold) {}
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
                if (current_spread > burstThreshold && static_cast<double>(current_spread) > sketch.ratioThreshold * previous_spread) {
                    true_burst_records.insert(BurstRecord_P2S_FC(flow_label, period_idx));
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
        cout << "True Burst: " << true_burst_records.size() << endl;
        cout << "Record Burst: " << sketch.detected_burst_records.size() << endl;
        const set<BurstRecord_P2S_FC>& detected_bursts = sketch.getDetectedBurstRecords();
        vector<BurstRecord_P2S_FC> true_bursts_vec(true_burst_records.begin(), true_burst_records.end());
        vector<BurstRecord_P2S_FC> detected_bursts_vec(detected_bursts.begin(), detected_bursts.end());
        vector<BurstRecord_P2S_FC> common_bursts_vec;
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
        cout << "Burst PR : " << fixed << setprecision(6) << PR << endl;
        cout << "Burst RR : " << fixed << setprecision(6) << RR << endl;
        cout << "Burst F1-score : " << fixed << setprecision(6) << F1 << endl;
        cout << endl;
    }
};
vector<double> Exp_P2S_FC(double alpha, int gamma, double wholeMemory, string priority_file_name, double burstTh, int T, int burstThreshold, string filename) {
    cout.precision(10);
    vector<vector<double>> results(5, vector<double>(1)); 
    int threeWholeMem = wholeMemory;
    double finePartMem = 0.2 * threeWholeMem;
    double wholeMem = 0.8 * threeWholeMem; 
    double bitmapMemRatio = 0.5; 
    int bits = static_cast<int>(wholeMem * bitmapMemRatio * 1024 * 8); 
    int sizeOfBkt = 4;
    int numOfBkt = static_cast<int>(finePartMem * 1024 / 13 / sizeOfBkt); 
    int d = 3;
    int w = static_cast<int>(wholeMem * (1 - bitmapMemRatio) * 1024 / (2 * d)); 
    double base = 100; 
    double step = 2; 
    Experiment_P2S_FC exp(bits, numOfBkt, sizeOfBkt, d, w, gamma, base, step, burstTh, burstThreshold);
    stringstream ss_alpha;
    ss_alpha << fixed << setprecision(6) << alpha;
    string alpha_str = ss_alpha.str();
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
