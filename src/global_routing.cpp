#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <algorithm>

#include <cstdio>
#include <cmath>
#include <ctime>

using namespace std;

typedef struct net {
    int sx, sy;
    int tx, ty;
} Net;

typedef struct edge {
    int idx;
    int direction;
} Edge;

// # grids in x and y
int grid_x, grid_y;
// horizontal and verical capacity
int h_cap, v_cap;
int num_nets;
vector<Net> nets;

// routing grids
vector<double> grids;
// routing directions, 0: top, 1: bottom, 2: left, 3: right
vector<int> grids_direction;

int max_grids;
// grids capacity usage
vector<int> grids_h_cap, grids_v_cap;
// nets pass through
vector<unordered_set<int>> grids_h_nets, grids_v_nets;
// overflow
vector<int> h_overflow, v_overflow;
// historical terms
vector<int> h_historical, v_historical;

// routing results of each net
vector<vector<int>> nets_routing;

// adaptive base cost function
double b = 1.0;

void ReadInputFile(string input_file)
{
    ifstream file;
    file.open(input_file);

    string str1, str2;
    file >> str1 >> grid_x >> grid_y;
    file >> str1 >> str2 >> v_cap;
    file >> str1 >> str2 >> h_cap;
    file >> str1 >> str2 >> num_nets;

    nets = vector<Net>(num_nets);
    max_grids = grid_x * grid_y;
    grids_h_cap = vector<int>((grid_x - 1) * grid_y, 0);
    grids_v_cap = vector<int>(grid_x * (grid_y - 1), 0);
    grids_h_nets = vector<unordered_set<int>>((grid_x - 1) * grid_y);
    grids_v_nets = vector<unordered_set<int>>(grid_x * (grid_y - 1));
    h_overflow = vector<int>((grid_x - 1) * grid_y, 0);
    v_overflow = vector<int>(grid_x * (grid_y - 1), 0);
    h_historical = vector<int>((grid_x - 1) * grid_y, 1);
    v_historical = vector<int>(grid_x * (grid_y - 1), 1);
    nets_routing = vector<vector<int>>(num_nets);

    int net_id, num_pins;
    for (int i = 0; i < num_nets; i++) {
        file >> str1 >> net_id >> num_pins;
        file >> nets[i].sx >> nets[i].sy >> nets[i].tx >> nets[i].ty;
    }

    file.close();
}

double GetRoutingCost(const int edge, const int direction)
{
    if (direction == 0)
        return b + h_historical[edge] * pow((double)(grids_h_cap[edge] + 1) / h_cap, 5);
    else
        return b + v_historical[edge] * pow((double)(grids_v_cap[edge] + 1) / v_cap, 5);
}

void WavePropagation(const int source, const int target, const int net_id, const bool monotonic)
{
    queue<int> wavefront;
    wavefront.push(source);

    int top_boundary = grid_y - 1;
    int bottom_boundary = 0;
    int left_boundary = 0;
    int right_boundary = grid_x - 1;
    if (monotonic) {
        int sx = source % grid_x;
        int sy = source / grid_x;
        int tx = target % grid_x;
        int ty = target / grid_x;
        if (sx < tx) {
            left_boundary = sx;
            right_boundary = tx;
        }
        else {
            left_boundary = tx;
            right_boundary = sx;
        }
        if (sy < ty) {
            top_boundary = ty;
            bottom_boundary = sy;
        }
        else {
            top_boundary = sy;
            bottom_boundary = ty;
        }
    }

    while (1) {
        queue<int> wavefront_temp;
        while (!wavefront.empty()) {
            int wf = wavefront.front();
            wavefront.pop();
            int wf_x = wf % grid_x;
            int wf_y = wf / grid_x;

            int top_wf = wf + grid_x;
            if (wf_y != top_boundary) {
                int edge_idx = wf_x + wf_y * grid_x;
                double cost = GetRoutingCost(edge_idx, 1);
                if (grids[top_wf] < 0 || grids[top_wf] > grids[wf] + cost) {
                    grids[top_wf] = grids[wf] + cost;
                    grids_direction[top_wf] = 1;
                    wavefront_temp.push(top_wf);
                }
            }

            int bottom_wf = wf - grid_x;
            if (wf_y != bottom_boundary) {
                int edge_idx = wf_x + (wf_y - 1) * grid_x;
                double cost = GetRoutingCost(edge_idx, 1);
                if (grids[bottom_wf] < 0 || grids[bottom_wf] > grids[wf] + cost) {
                    grids[bottom_wf] = grids[wf] + cost;
                    grids_direction[bottom_wf] = 0;
                    wavefront_temp.push(bottom_wf);
                }
            }

            int left_wf = wf - 1;
            if (wf_x != left_boundary) {
                int edge_idx = (wf_x - 1) + wf_y * (grid_x - 1);
                double cost = GetRoutingCost(edge_idx, 0);
                if (grids[left_wf] < 0 || grids[left_wf] > grids[wf] + cost) {
                    grids[left_wf] = grids[wf] + cost;
                    grids_direction[left_wf] = 3;
                    wavefront_temp.push(left_wf);
                }
            }

            int right_wf = wf + 1;
            if (wf_x != right_boundary) {
                int edge_idx = wf_x + wf_y * (grid_x - 1);
                double cost = GetRoutingCost(edge_idx, 0);
                if (grids[right_wf] < 0 || grids[right_wf] > grids[wf] + cost) {
                    grids[right_wf] = grids[wf] + cost;
                    grids_direction[right_wf] = 2;
                    wavefront_temp.push(right_wf);
                }
            }
        }

        if (wavefront_temp.size() == 0)
            break;

        wavefront = wavefront_temp;
    }
}

void Retrace(const int source, const int target, const int net_id)
{
    nets_routing[net_id].emplace_back(target);

    int wf = target;

    int prev_dir = -1; // 0: horizontal, 1: vertical
    while (wf != source) {
        int wf_x = wf % grid_x;
        int wf_y = wf / grid_x;

        int direction = grids_direction[wf];
        if (direction == 0) {
            int edge_idx = wf_x + wf_y * grid_x;
            grids_v_cap[edge_idx]++;
            grids_v_nets[edge_idx].insert(net_id);
            if (prev_dir == 0)
                nets_routing[net_id].emplace_back(wf);
            prev_dir = 1;

            wf += grid_x;
        }
        else if (direction == 1) {
            int edge_idx = wf_x + (wf_y - 1) * grid_x;
            grids_v_cap[edge_idx]++;
            grids_v_nets[edge_idx].insert(net_id);
            if (prev_dir == 0)
                nets_routing[net_id].emplace_back(wf);
            prev_dir = 1;

            wf -= grid_x;
        }
        else if (direction == 2) {
            int edge_idx = (wf_x - 1) + wf_y * (grid_x - 1);
            grids_h_cap[edge_idx]++;
            grids_h_nets[edge_idx].insert(net_id);
            if (prev_dir == 1)
                nets_routing[net_id].emplace_back(wf);
            prev_dir = 0;

            wf -= 1;
        }
        else if (direction == 3) {
            int edge_idx = wf_x + wf_y * (grid_x - 1);
            grids_h_cap[edge_idx]++;
            grids_h_nets[edge_idx].insert(net_id);
            if (prev_dir == 1)
                nets_routing[net_id].emplace_back(wf);
            prev_dir = 0;

            wf += 1;
        }
    }

    nets_routing[net_id].emplace_back(source);
}

void MazeRouting(const int net, const bool monotonic)
{
    grids = vector<double>(grid_x * grid_y, -1);
    grids_direction = vector<int>(grid_x * grid_y, -1);

    int source = nets[net].sx + nets[net].sy * grid_x;
    int target = nets[net].tx + nets[net].ty * grid_x;
    grids[source] = 0;

    WavePropagation(source, target, net, monotonic);
    Retrace(source, target, net);
}

void InitialRouting()
{
    for (int i = 0; i < num_nets; i++)
        MazeRouting(i, true);
}

int CalculateOverflow()
{
    int total_overflow = 0;

    for (int y = 0; y < grid_y; y++) {
        for (int x = 0; x < grid_x - 1; x++) {
             int idx = x + y * (grid_x - 1);
             h_overflow[idx] = max(grids_h_cap[idx] - h_cap, 0);
             total_overflow += h_overflow[idx];
        }
    }

    for (int y = 0; y < grid_y - 1; y++) {
        for (int x = 0; x < grid_x; x++) {
            int idx = x + y * grid_x;
            v_overflow[idx] = max(grids_v_cap[idx] - v_cap, 0);
            total_overflow += v_overflow[idx];
        }
    }

    return total_overflow;
}

struct CompareOverflow {
    bool operator() (const Edge e1, const Edge e2) {
        int overflow_1 = e1.direction == 0 ? h_overflow[e1.idx] : v_overflow[e1.idx];
        int overflow_2 = e2.direction == 0 ? h_overflow[e2.idx] : v_overflow[e2.idx];
        return overflow_1 > overflow_2;
    }
};

void GetRipupEdges(priority_queue<Edge, vector<Edge>, CompareOverflow> &ripup_queue)
{
    Edge e;

    for (int y = 0; y < grid_y; y++) {
        for (int x = 0; x < grid_x - 1; x++) {
             int idx = x + y * (grid_x - 1);
             if (h_overflow[idx] > 0) {
                 h_historical[idx]++;
                 e.idx = idx;
                 e.direction = 0;
                 ripup_queue.push(e);
             }
        }
    }

    for (int y = 0; y < grid_y - 1; y++) {
        for (int x = 0; x < grid_x; x++) {
            int idx = x + y * grid_x;
            if (v_overflow[idx] > 0) {
                v_historical[idx]++;
                e.idx = idx;
                e.direction = 1;
                ripup_queue.push(e);
            }
        }
    }
}

struct CompareBoundingBox {
    bool operator() (const int net1, const int net2) {
        int wl1 = abs(nets[net1].sx - nets[net1].tx) + abs(nets[net1].sy - nets[net1].ty);
        int wl2 = abs(nets[net2].sx - nets[net2].tx) + abs(nets[net2].sy - nets[net2].ty);
        return wl1 < wl2;
    }
};

void RipupReroute(const int ripup_edge, const int ripup_direction)
{
    // rip up
    unordered_set<int> ripup_nets;
    if (ripup_direction == 0)
        ripup_nets = grids_h_nets[ripup_edge];
    else
        ripup_nets = grids_v_nets[ripup_edge];

    priority_queue<int, vector<int>, CompareBoundingBox> reroute_order;
    for (const int net : ripup_nets) {
        int size = nets_routing[net].size() - 1;
        for (int i = 0; i < size; i++) {
            int x1 = nets_routing[net][i] % grid_x;
            int y1 = nets_routing[net][i] / grid_x;
            int x2 = nets_routing[net][i + 1] % grid_x;
            int y2 = nets_routing[net][i + 1] / grid_x;

            int s, t, direction;
            if (x1 == x2) {
                s = y1 < y2 ? y1 : y2;
                t = y1 < y2 ? y2 : y1;
                direction = 1;
            }
            else {
                s = x1 < x2 ? x1 : x2;
                t = x1 < x2 ? x2 : x1;
                direction = 0;
            }

            if (direction == 0) {
                for (int j = s; j < t; j++) {
                    int idx = j + y1 * (grid_x - 1);
                    grids_h_nets[idx].erase(net);
                    grids_h_cap[idx]--;
                }
            }
            else {
                for (int j = s; j < t; j++) {
                    int idx = x1 + j * grid_x;
                    grids_v_nets[idx].erase(net);
                    grids_v_cap[idx]--;
                }
            }
        }
/*
        for (int y = 0; y < grid_y; y++) {
            for (int x = 0; x < grid_x - 1; x++) {
                int idx = x + y * (grid_x - 1);
                if (grids_h_nets[idx].find(net) != grids_h_nets[idx].end()) {
                    grids_h_nets[idx].erase(net);
                    grids_h_cap[idx]--;
                }
            }
        }
        for (int y = 0; y < grid_y - 1; y++) {
            for (int x = 0; x < grid_x; x++) {
                int idx = x + y * grid_x;
                if (grids_v_nets[idx].find(net) != grids_v_nets[idx].end()) {
                    grids_v_nets[idx].erase(net);
                    grids_v_cap[idx]--;
                }
            }
        }
*/

        nets_routing[net] = vector<int>();
        reroute_order.push(net);
    }

    // reroute
    while (!reroute_order.empty()) {
        int net = reroute_order.top();
        reroute_order.pop();
        MazeRouting(net, false);
    }
}

void GlobalRouting()
{
    clock_t start = clock();
    float time_elapsed;

    InitialRouting();

    int total_overflow = CalculateOverflow();
    int iter = 0;
    while (total_overflow > 0) {
        iter++;
        b = 1.0 - exp(-5 * exp(-0.1 * iter));

        priority_queue<Edge, vector<Edge>, CompareOverflow> ripup_queue;
        GetRipupEdges(ripup_queue);

        while (total_overflow > 0 && !ripup_queue.empty()) {
            int edge_idx = ripup_queue.top().idx;
            int direction = ripup_queue.top().direction;
            ripup_queue.pop();
            if (direction == 0 && grids_h_cap[edge_idx] > h_cap || direction == 1 && grids_v_cap[edge_idx] > v_cap)
                RipupReroute(edge_idx, direction);
            total_overflow = CalculateOverflow();

            time_elapsed = (clock() - start) / CLOCKS_PER_SEC;
            if (time_elapsed >= 590)
                total_overflow = 0;
        }
    }
}

void OutputFile(string output_file)
{
    ofstream file;
    file.open(output_file);

    for (int i = 0; i < num_nets; i++) {
        file << "net" << i << " " << i << '\n';

        int size = nets_routing[i].size() - 1;
        for (int j = size; j > 0; j--) {
            file << "(" << nets_routing[i][j] % grid_x << ", " << nets_routing[i][j] / grid_x << ", 1)";
            file << "-";
            file << "(" << nets_routing[i][j - 1] % grid_x << ", " << nets_routing[i][j - 1] / grid_x << ", 1)";
            file << '\n';
        }

        file << "!\n";
    }

    file.close();
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        cout << "Usage: \n";
        cout << "    ./route <path/to/input_file> <path/to/output_file>\n";
        exit(0);
    }
    string input_file = argv[1];
    string output_file = argv[2];

    ReadInputFile(input_file);

    GlobalRouting();

    OutputFile(output_file);

    return 0;
}
