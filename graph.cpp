#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm> 

using namespace std;
const int dirX[] = {-1, 1, 0, 0};
const int dirY[] = {0, 0, -1, 1};

class graph {
public:
  int n, m, sx, sy, gx, gy;
  double *deg;
  vector<string> map;
  vector<vector<int>> adj;

  graph(int n, int m, string map[]) : n(n), m(m){
    deg = new double[n*m];
    adj.resize(n*m);
    fill(deg, deg + n, 0.0);

    for (int i = 0; i < n; ++i) {
      this->map.push_back(map[i]);
    }
    
    bfs();
  }

  ~graph() {
    delete[] deg;
  }
  
  void bfs() {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        if (map[i][j] == '#') continue;
        else if (map[i][j] == 'T') {
          sx = i; sy = j;
        } else if (map[i][j] == 'W') {
          gx = i; gy = j;
        }
        int u = i * m + j;
        for (int d = 0; d < 4; ++d) {
          int ni = i + dirX[d];
          int nj = j + dirY[d];
          if (ni >= 0 && ni < n && nj >= 0 && nj < m && map[ni][nj] != '#') {
            int v = ni * m + nj;
            adj[u].push_back(v);
            deg[u]++;
          }
        }
      }
    }
  }
};

int main() {
  int t, n, m;
  scanf("%d", &t);
  string map[t];
  for (int i = 0; i < t; ++i) {
    cin >> n >> m;
    for (int j = 0; j < n; ++j) {
      cin >> map[j];
    }
    graph g(n, m, map);
  }


  return 0;
}
