#include <bits/stdc++.h>
#include "MurmurHash3.h"
#include "Tool.cpp"
using namespace std;
uint32_t hash1(uint32_t val, uint32_t seed) {
    uint32_t hash_value;
    MurmurHash3_x86_32(reinterpret_cast<const void*>(&val), sizeof(val), seed, &hash_value);
    return hash_value;
}
class HyperLogLog {
public:
    int m;
    vector<int> registers;
    HyperLogLog(vector<int> reg) : registers(reg), m(reg.size()) {}
    double estimate() {
        double alpha_m = 0.673;
        double sum = 0.0;
        for (int r : registers)
            sum += 1.0 / (1 << r);
        double raw_est = alpha_m * m * m / sum;
        int V = count(registers.begin(), registers.end(), 0);
        if (raw_est <= 2.5 * m && V > 0)
            return m * log(static_cast<double>(m) / V);
        else if (raw_est > (1.0 / 30.0) * (1LL << 32))
            return -(1LL << 32) * log(1.0 - raw_est / (1LL << 32));
        return raw_est;
    }
};
class VHLL {
public:
    int num_phy_registers, num_virt_registers;
    vector<int> pyh_registers;
    vector<uint32_t> seedArray;
    uint32_t seed_fe;
    int count_rank_10 = 0;
    int count_src_rank_10 = 0;
    uint32_t seed_tmp = generateSeeds(1)[0];
    VHLL(int phy, int virt) : num_phy_registers(phy), num_virt_registers(virt) {
        pyh_registers.resize(phy, 0);
        seedArray = generateSeeds(virt);
        seed_fe = generateSeeds(1)[0];
    }
    void insert(uint32_t src, uint32_t dst) {
        uint32_t fe = src ^ dst;
        uint32_t hashfe = hash1(fe, seed_fe);
        int b = log2(num_virt_registers);
        int v_pos = hashfe & (num_virt_registers - 1);
        uint32_t remain = hashfe >> b;
        int rank = (hashfe == 0) ? 32 : __builtin_clz(hashfe) + 1;
        int phy_pos = hash1(src, seedArray[v_pos]) % num_phy_registers;
        pyh_registers[phy_pos] = max(pyh_registers[phy_pos], rank);
    }
    double query(uint32_t src) {
        vector<int> regs;
        for (int i = 0; i < num_virt_registers; ++i) {
            int pos = hash1(src, seedArray[i]) % num_phy_registers;
            regs.push_back(pyh_registers[pos]);
        }
        return HyperLogLog(regs).estimate();
    }
};
vector<double> Exp_VHLL(double alpha, int gamma, double wholeMemory, string priority_file_name) {
    string input_file = priority_file_name + "_" + to_string(alpha) + "_" + to_string(gamma) + "_.txt";
    ifstream file(input_file);
    if (!file.is_open()) {
        return {};
    }
    VHLL vhll(wholeMemory * 1024 * 8 / 4, 16);
    unordered_map<uint32_t, set<uint32_t>> real_map;
    unordered_map<uint32_t, int> flowPriority_dict;
    vector<tuple<uint32_t, uint32_t, int>> parsed_data;
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        uint32_t src, dst;
        int priority;
        ss >> src >> dst >> priority;
        parsed_data.push_back({src, dst, priority}); 
    }
    file.close();
    auto start_time = chrono::high_resolution_clock::now();
    for (const auto& [src, dst, priority] : parsed_data) {
        vhll.insert(src, dst); 
    }
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time = end_time - start_time;
    double throughput = static_cast<double>(parsed_data.size()) / elapsed_time.count() / 1000000;
    cout << "VHLL Insertion Throughput: " << throughput << " MPPS" << endl;
    for (const auto& [src, dst, priority] : parsed_data) {
        real_map[src].insert(dst);
        flowPriority_dict[src] = priority;
    }
    file.close();
    double highPriorityErrorSUM = 0;
    int highPriorityFlowNum = 0;
    double lowPriorityErrorSUM = 0;
    int lowPriorityFlowNum = 0;
    for (const auto& [src, dst_set] : real_map) {
        int priority = flowPriority_dict[src];
        double real = dst_set.size();
        double est = vhll.query(src);
        if (priority > gamma * 0.75) {
            highPriorityFlowNum++;
            highPriorityErrorSUM += abs(est - real) / real;
        } else {
            lowPriorityFlowNum++;
            lowPriorityErrorSUM += abs(est - real) / real;;
        }
    }
    double highPriorityARE = highPriorityErrorSUM / highPriorityFlowNum;
    double lowPriorityARE = lowPriorityErrorSUM / lowPriorityFlowNum;
    double allFlowsARE = (highPriorityErrorSUM + lowPriorityErrorSUM) / (highPriorityFlowNum + lowPriorityFlowNum);
    cout << "highPriorityARE : " << fixed << setprecision(6) << highPriorityARE << endl;
    cout << "lowPriorityARE : " << fixed << setprecision(6) << lowPriorityARE << endl;
    cout << "allFlowsARE : " << fixed << setprecision(6) << allFlowsARE << endl;
    vector<double> ret(5);
    ret[0] = highPriorityARE;
    ret[1] = lowPriorityARE;
    ret[2] = allFlowsARE;
    ret[3] = 0;
    ret[4] = throughput;
    return ret;
}
