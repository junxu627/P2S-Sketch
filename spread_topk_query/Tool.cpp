#ifndef tool_H
#define tool_H
#include<bits/stdc++.h>
#include "MurmurHash3.h" 
#include <immintrin.h>     
#include <x86intrin.h>     
using namespace std;
vector<uint32_t> generateSeeds(int num) {
    vector<uint32_t> seeds;
    set<uint32_t> unique_seeds;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint32_t> dist(67297, UINT32_MAX);
    while (unique_seeds.size() < num) {
        uint32_t seed = dist(gen);
        unique_seeds.insert(seed);
    }
    seeds.assign(unique_seeds.begin(), unique_seeds.end());
    return seeds;
}
uint32_t ipv4_to_int(const string& ip_str) {
    uint32_t ip_int = 0;
    stringstream ss(ip_str);
    string segment;
    for (int i = 0; i < 4; ++i) {
        getline(ss, segment, '.');
        ip_int = (ip_int << 8) | stoul(segment);
    }
    return ip_int;
}
void convert_and_save_ips(const string& InputFilename, const string& OutputFilename) {
    ifstream infile(InputFilename);
    ofstream outfile(OutputFilename);
    string line;
    if (!infile.is_open()) {
        return;
    }
    if (!outfile.is_open()) {
        return;
    }
    int line_count = 0;
    while (getline(infile, line)) {
        if (line.size() >= 32) {
            continue;
        }
        size_t first_space = line.find(' ');
        if (first_space == string::npos) {
            continue; 
        }
        string src_ip_str = line.substr(0, first_space);
        string dst_ip_str = line.substr(first_space + 1); 
        if (dst_ip_str.empty()) {
            continue;
        }
        uint32_t src_ip_int = ipv4_to_int(src_ip_str);
        uint32_t dst_ip_int = ipv4_to_int(dst_ip_str);
        outfile << src_ip_int << " " << dst_ip_int << "\n";
        line_count++;
    }
    infile.close();
    outfile.close();
}
void generate_prioritized_data(const string& InputConvertedFilename, double alpha, int gamma, const string& OutputPrioritizedFilename) {
    vector<uint32_t> seeds = generateSeeds(1);
    uint32_t seed = seeds[0];
    uint32_t max_uint32 = numeric_limits<uint32_t>::max();
    unordered_map<uint32_t, int> flowPriorityMap; 
    double bPart = 0;
    for (int i = 1; i < gamma + 1; ++i) {
        bPart += pow((1.0 / i), alpha);
    }
    double sum = 0;
    vector<double> sumVec;
    vector<double> priorityRankVec; 
    sumVec.push_back(0.0); 
    priorityRankVec.push_back(0.0); 
    map<int, double> theoretical_probabilities;
    for (int x = 1; x < gamma + 1; ++x) {
        double aPart = pow(static_cast<double>(x), alpha);
        double res = 1 / (aPart * bPart);
        sum += res;
        sumVec.push_back(sum);
        priorityRankVec.push_back(sum * max_uint32);
        theoretical_probabilities[x] = res;
    }
    ifstream infile(InputConvertedFilename);
    ofstream outfile(OutputPrioritizedFilename);
    string line;
    if (!infile.is_open()) {
        return;
    }
    if (!outfile.is_open()) {
        return;
    }
    int line_count = 0;
    map<int, long long> priorityCounts;
    for (int i = 1; i <= gamma; ++i) {
        priorityCounts[i] = 0; 
    }
    while (getline(infile, line)) {
        stringstream ss(line);
        uint32_t src_ip, dst_ip;
        ss >> src_ip >> dst_ip;
        int priority = 0;
        auto it = flowPriorityMap.find(src_ip);
        if (it != flowPriorityMap.end()) {
            priority = it->second; 
        } else {
            uint32_t hashValue;
            MurmurHash3_x86_32(reinterpret_cast<const void*>(&src_ip), sizeof(src_ip), 627, &hashValue);
            for (int j = 1; j < priorityRankVec.size(); ++j) {
                if (hashValue <= priorityRankVec[j]) {
                    priority = j;
                    flowPriorityMap[src_ip] = j; 
                    break;
                }
            }
        }
        outfile << src_ip << " " << dst_ip << " " << priority << "\n";
        priorityCounts[priority]++; 
        line_count++;
    }
    infile.close();
    outfile.close();
}
void calculate_flow_cardinalities(const string& InputConvertedFilename, const string& OutputFlowCardinalitiesFilename) {
    ifstream infile(InputConvertedFilename);
    string line;
    unordered_map<uint32_t, set<uint32_t>> flow_destination_ips;
    if (!infile.is_open()) {
        return;
    }
    int line_count = 0;
    while (getline(infile, line)) {
        stringstream ss(line);
        uint32_t src_ip, dst_ip;
        if (!(ss >> src_ip >> dst_ip)) {
            continue;
        }
        flow_destination_ips[src_ip].insert(dst_ip);
        line_count++;
    }
    infile.close();
    vector<pair<uint32_t, size_t>> flow_cardinalities;
    for (const auto& pair : flow_destination_ips) {
        flow_cardinalities.push_back({pair.first, pair.second.size()});
    }
    sort(flow_cardinalities.begin(), flow_cardinalities.end(), 
              [](const pair<uint32_t, size_t>& a, const pair<uint32_t, size_t>& b) {
                  return a.second > b.second; 
              });
    ofstream outfile(OutputFlowCardinalitiesFilename);
    if (!outfile.is_open()) {
        return;
    }
    for (const auto& pair : flow_cardinalities) {
        outfile << pair.first << " " << pair.second << "\n";
    }
    outfile.close();
}
void generate_volume_dependent_prioritized_data(
    const string& InputOriginalConvertedFilename,
    const string& InputCardinalitiesFilename,
    double alpha,
    int gamma,
    const string& OutputPrioritizedFilename,
    const string& OutputUniqueFlowsPrioritiesFilename) {
    unordered_map<uint32_t, int> flowPriorityMap; 
    double bPart = 0;
    for (int k = 1; k < gamma + 1; ++k) {
        bPart += pow((1.0 / k), alpha);
    }
    vector<double> sumVec; 
    sumVec.push_back(0.0); 
    double current_sum_prob = 0.0;
    map<int, double> theoretical_probabilities;
    for (int k = gamma; k > 0; --k) {
        double prob_rank_k = 1 / (pow(static_cast<double>(k), alpha) * bPart);
        current_sum_prob += prob_rank_k;
        sumVec.push_back(current_sum_prob);
        theoretical_probabilities[k] = prob_rank_k;
    }
    vector<uint32_t> sorted_flows_by_cardinality; 
    ifstream card_infile(InputCardinalitiesFilename);
    if (!card_infile.is_open()) {
        return;
    }
    string line;
    uint32_t src_ip_temp;
    size_t cardinality_temp;
    while (getline(card_infile, line)) {
        stringstream ss(line);
        if (!(ss >> src_ip_temp >> cardinality_temp)) {
             continue;
        }
        sorted_flows_by_cardinality.push_back(src_ip_temp);
    }
    card_infile.close();
    int flowVecSize = sorted_flows_by_cardinality.size();
    if (flowVecSize == 0) {
        return;
    }
    vector<double> priorityRankThresholds; 
    priorityRankThresholds.push_back(0.0); 
    for (int k = 1; k < gamma + 1; ++k) {
        priorityRankThresholds.push_back(sumVec[k] * flowVecSize);
    }
    for (int i = 0; i < flowVecSize; ++i) { 
        for (int j = 1; j <= gamma; ++j) { 
            if (i >= priorityRankThresholds[j-1] && i < priorityRankThresholds[j]) {
                flowPriorityMap[sorted_flows_by_cardinality[i]] = gamma - j + 1;
                break;
            }
        }
    }
    ifstream original_infile(InputOriginalConvertedFilename);
    ofstream outfile(OutputPrioritizedFilename);
    if (!original_infile.is_open()) {
        return;
    }
    if (!outfile.is_open()) {
        return;
    }
    int line_count_output = 0;
    outfile << "--- Theoretical Zipf Probability Distribution (alpha=" << alpha << ", gamma=" << gamma << ") ---\n";
    for (int i = 1; i <= gamma; ++i) {
        outfile << "Priority " << i << " (Zipf Rank " << i << "): P = " 
                << fixed << setprecision(10) << theoretical_probabilities[i] << "\n";
    }
    map<int, long long> priorityCounts;
    for (int i = 1; i <= gamma; ++i) {
        priorityCounts[i] = 0; 
    }
    while (getline(original_infile, line)) {
        stringstream ss(line);
        uint32_t src_ip, dst_ip;
        if (!(ss >> src_ip >> dst_ip)) {
            continue;
        }
        auto it = flowPriorityMap.find(src_ip);
        if (it != flowPriorityMap.end()) {
            outfile << src_ip << " " << dst_ip << " " << it->second << "\n";
        } else {
            continue; 
        }
        line_count_output++;
        priorityCounts[it->second]++; 
    }
    original_infile.close();
    outfile << "\n--- Priority Distribution Statistics ---\n";
    for (int i = 1; i <= gamma; ++i) {
        outfile << "Priority " << i << ": " << priorityCounts[i] << " occurrences\n";
    }
    outfile << "Total lines processed: " << line_count_output << "\n";
    outfile << "--------------------------------------\n";
    outfile.close();
}
void generate_volume_dependent_prioritized_data_new(
    const string& InputOriginalConvertedFilename,
    const string& InputCardinalitiesFilename,
    double alpha,
    int gamma,
    const string& OutputPrioritizedFilename,
    const string& OutputUniqueFlowsPrioritiesFilename) { 
    unordered_map<uint32_t, int> flowPriorityMap; 
    unordered_map<uint32_t, size_t> flowCardinalityMap; 
    double bPart = 0;
    for (int k = 1; k < gamma + 1; ++k) {
        bPart += pow((1.0 / k), alpha);
    }
    vector<double> sumVec; 
    sumVec.push_back(0.0); 
    double current_sum_prob = 0.0;
    map<int, double> theoretical_probabilities; 
    for (int k = gamma; k > 0; --k) { 
        double prob_rank_k = 1 / (pow(static_cast<double>(k), alpha) * bPart);
        current_sum_prob += prob_rank_k;
        sumVec.push_back(current_sum_prob); 
        theoretical_probabilities[k] = prob_rank_k; 
    }
    vector<pair<uint32_t, size_t>> raw_flows_with_cardinality; 
    ifstream card_infile(InputCardinalitiesFilename);
    if (!card_infile.is_open()) {
        return;
    }
    string line;
    uint32_t src_ip_temp;
    size_t cardinality_temp;
    while (getline(card_infile, line)) {
        stringstream ss(line);
        if (!(ss >> src_ip_temp >> cardinality_temp)) {
             continue;
        }
        raw_flows_with_cardinality.push_back({src_ip_temp, cardinality_temp});
        flowCardinalityMap[src_ip_temp] = cardinality_temp; 
    }
    card_infile.close();
    sort(raw_flows_with_cardinality.begin(), raw_flows_with_cardinality.end(), 
              [](const pair<uint32_t, size_t>& a, const pair<uint32_t, size_t>& b) {
                  return a.second > b.second; 
              });
    vector<uint32_t> sorted_flows_by_cardinality;
    for(const auto& p : raw_flows_with_cardinality) {
        sorted_flows_by_cardinality.push_back(p.first);
    }
    int flowVecSize = sorted_flows_by_cardinality.size();
    if (flowVecSize == 0) {
        return;
    }
    vector<size_t> priorityRankThresholds; 
    priorityRankThresholds.push_back(0); 
    for (int k = 1; k < gamma + 1; ++k) { 
        priorityRankThresholds.push_back(static_cast<size_t>(sumVec[k] * flowVecSize));
    }
    priorityRankThresholds[gamma] = flowVecSize; 
    double prev_sum_prob_debug = 0.0;
    size_t prev_threshold_debug = 0;
    for (int j = 1; j < sumVec.size(); ++j) {
        double current_rank_prob_debug = sumVec[j] - prev_sum_prob_debug;
        size_t current_flows_in_rank = priorityRankThresholds[j] - prev_threshold_debug;
        cout << "  Priority " << (gamma - j + 1) << " (Zipf Rank " << j << "): "
                  << "Prob = " << fixed << setprecision(10) << current_rank_prob_debug
                  << ", CumProb = " << sumVec[j]
                  << ", Threshold (Flow Index) = " << priorityRankThresholds[j]
                  << ", Flows in this segment: " << current_flows_in_rank << endl;
        prev_sum_prob_debug = sumVec[j];
        prev_threshold_debug = priorityRankThresholds[j];
    }
    cout << endl;
    for (int i = 0; i < flowVecSize; ++i) { 
        for (int j = 1; j <= gamma; ++j) { 
            if (i >= priorityRankThresholds[j-1] && i < priorityRankThresholds[j]) {
                flowPriorityMap[sorted_flows_by_cardinality[i]] = gamma - j + 1;
                break;
            }
        }
    }
    ofstream unique_flows_outfile(OutputUniqueFlowsPrioritiesFilename);
    if (!unique_flows_outfile.is_open()) {
        return;
    }
    map<int, long long> uniquePriorityCounts;
    for (int i = 1; i <= gamma; ++i) {
        uniquePriorityCounts[i] = 0;
    }
    for (const auto& flow_ip : sorted_flows_by_cardinality) {
        int priority = flowPriorityMap[flow_ip];
        size_t cardinality = flowCardinalityMap[flow_ip];
        unique_flows_outfile << flow_ip << " " << cardinality << " " << priority << "\n";
        uniquePriorityCounts[priority]++; 
    }
    unique_flows_outfile.close();
    ifstream original_infile(InputOriginalConvertedFilename);
    ofstream outfile(OutputPrioritizedFilename);
    if (!original_infile.is_open()) {
        return;
    }
    if (!outfile.is_open()) {
        return;
    }
    int line_count_output = 0;
    map<int, long long> totalPriorityCounts;
    for (int i = 1; i <= gamma; ++i) {
        totalPriorityCounts[i] = 0; 
    }
    while (getline(original_infile, line)) {
        stringstream ss(line);
        uint32_t src_ip, dst_ip;
        if (!(ss >> src_ip >> dst_ip)) {
            continue;
        }
        auto it = flowPriorityMap.find(src_ip);
        if (it != flowPriorityMap.end()) {
            outfile << src_ip << " " << dst_ip << " " << it->second << "\n";
            totalPriorityCounts[it->second]++; 
        } else {
            continue; 
        }
        line_count_output++;
    }
    original_infile.close();
    outfile.close();
}
inline int SIMD_MIN_4(uint32_t* LBT) {
    int matched;
    const uint32_t mask_base = 0x7FFFFFFF;
    const __m128i* counters = reinterpret_cast<const __m128i*>(LBT);
    __m128i masks = _mm_set1_epi32(mask_base);
    __m128i results = _mm_and_si128(counters[0], masks); 
    __m128i min1 = _mm_shuffle_epi32(results, _MM_SHUFFLE(0, 0, 3, 2));
    __m128i min2 = _mm_min_epi32(results, min1);
    __m128i min3 = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128i min4 = _mm_min_epi32(min2, min3);
    uint32_t min_counter_val = _mm_cvtsi128_si32(min4); 
    const __m128i ct_item = _mm_set1_epi32(min_counter_val);
    __m128i ct_a_comp = _mm_cmpeq_epi32(ct_item, results); 
    matched = _mm_movemask_epi8(ct_a_comp); 
    return _tzcnt_u32((uint32_t)matched) >> 2;
}
#endif
