#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;

mt19937 gen;
uniform_int_distribution<int> uniform(0, 1);

// Вычисление значения на ребре
int edge_num(int layer, int cur_v, vector<vector<int>> &p, vector<vector<int>> &activity, int removed_v) {
    int res = 0;
    for (int i = 0; i < activity[layer].size(); ++i) {
        int cur_bit = (cur_v >> i) & 1;
        res ^= cur_bit * p[activity[layer][i]][layer - 1];
    }
    if (removed_v >= 0) {
        res ^= removed_v;
    }
    return res;
}

// Кодирование перемножением матрицы на вектор
vector<int> encode(int n, vector<int> &input, const vector<vector<int>> &g) {
    vector<int> result(n);
    for (int i = 0; i < g.size(); ++i) {
        for (int j = 0; j < g[i].size(); ++j) {
            result[j] ^= input[i] * g[i][j];
        }
    }
    return result;
}

// Декодирование
vector<int> decode(int n, vector<double> &input, const vector<vector<vector<pair<int, int>>>> &graph) {
    // Поиск лучшего пути
    vector<vector<pair<double, int>>> dynamic(n + 1);
    dynamic[0].push_back({0, -1});
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j < graph[i].size(); ++j) {
            double best = 0;
            int best_i = 0;
            for (int l = 0; l < graph[i][j].size(); ++l) {
                pair<int, int> in_e = graph[i][j][l];
                int edge_val = 1 - 2 * in_e.second;
                double sum = edge_val * input[i - 1] + dynamic[i - 1][in_e.first].first;
                if (best < sum || l == 0) {
                    best = sum;
                    best_i = l;
                }
            }
            dynamic[i].push_back({best, best_i});
        }
    }
    // Восстановление ответа
    int cur_v = 0;
    int cur_l = n;
    vector<int> result(n);
    while (cur_l > 0) {
        pair<int, int> edge = graph[cur_l][cur_v][dynamic[cur_l][cur_v].second];
        cur_v = edge.first;
        result[cur_l - 1] = edge.second;
        --cur_l;
    }
    return result;
}

// Симуляция
double simulate(int n, int k, int it, int max_err, double div, const vector<vector<int>> &g, const vector<vector<vector<pair<int, int>>>> &graph) {
    // Дисперсия
    double d = 0.5 * pow(10, -div / 10.0) * n / k;
    normal_distribution<double> distribution(0, sqrt(d));
    int cnt_err = 0;
    for (int i = 0; i < it; ++i) {
        // Генерация вектора
        vector<int> original(k);
        for (int j = 0; j < k; ++j) {
            original[j] = uniform(gen);
        }

        // Кодирование
        vector<int> encoded = encode(n, original, g);

        // Внесение шума
        vector<double> with_noise(n);
        for (int j = 0; j < n; ++j) {
            with_noise[j] = (1 - 2 * encoded[j]) + distribution(gen);
        }

        // Декодирование
        vector<int> res = decode(n, with_noise, graph);
        for (int j = 0; j < n; ++j) {
            if (encoded[j] != res[j]) {
                ++cnt_err;
                break;
            }
        }
        if (cnt_err == max_err) {
            return cnt_err * 1.0 / (i + 1);
        }
    }
    return cnt_err * 1.0 / it;
}

int main() {
    ifstream cin("input.txt");
    ofstream cout("output.txt");

    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    int n, k;
    cin >> n >> k;

    vector<vector<int>> g(k, vector<int>(n));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> g[i][j];
        }
    }

    // Приведение к МСФ
    vector<vector<int>> p = g;
    // Приведение начал
    for (int i = 0; i < k; ++i) {
        int best = n + 1;
        int best_i = -1;
        for (int l = i; l < k; ++l) {
            for (int t = i; t < n; ++t) {
                if (p[l][t] > 0) {
                    if (t < best) {
                        best = t;
                        best_i = l;
                        break;
                    }
                }
            }
        }
        for (int t = best; t < n; ++t) {
            swap(p[best_i][t], p[i][t]);
        }
        for (int l = i + 1; l < k; ++l) {
            if (p[l][best] > 0) {
                for (int t = 0; t < n; ++t) {
                    p[l][t] ^= p[i][t];
                }
            }
        }
    }
    // Приведение концов
    for (int i = k - 1; i >= 0; --i) {
        int end;
        for (int j = n - 1; j >= 0; --j) {
            if (p[i][j] > 0) {
                end = j;
                break;
            }
        }
        for (int l = i - 1; l >= 0; --l) {
            if (p[l][end] > 0) {
                for (int t = 0; t < n; ++t) {
                    p[l][t] ^= p[i][t];
                }
            }
        }
    }

    // Предподсчет активных элементов
    vector<vector<int>> activity(n + 1);
    for (int i = 0; i < k; ++i) {
        int begin, end;
        for (int j = 0; j < n; ++j) {
            if (p[i][j] > 0) {
                begin = j;
                break;
            }
        }
        for (int j = n - 1; j >= 0; --j) {
            if (p[i][j] > 0) {
                end = j;
                break;
            }
        }
        for (int j = begin; j < end; ++j) {
            activity[j + 1].push_back(i);
        }
    }

    // Подсчет количества элементов в решетке
    vector<int> vertexes(n + 1);
    for (int j = 0; j <= n; ++j) {
        vertexes[j] = 1 << (activity[j].size());
    }

    // Построение решетки
    vector<vector<vector<pair<int, int>>>> graph(n + 1);
    for (int i = 1; i <= n; ++i) {
        int inserted_v = -1;
        int removed_v = -1;

        // Детекция изменения активных строк
        int k1 = 0;
        int k2 = 0;
        while (k1 < activity[i - 1].size() || k2 < activity[i].size()) {
            if (k1 == activity[i - 1].size()) {
                inserted_v = k2++;
            } else if (k2 == activity[i].size()) {
                removed_v = k1++;
            } else {
                if (activity[i - 1][k1] == activity[i][k2]) {
                    ++k1;
                    ++k2;
                } else if (k1 + 1 < activity[i - 1].size() && activity[i - 1][k1 + 1] == activity[i][k2]) {
                    removed_v = k1++;
                } else if (k2 + 1 < activity[i].size() && activity[i - 1][k1] == activity[i][k2 + 1]) {
                    inserted_v = k2++;
                } else {
                    inserted_v = k2++;
                    removed_v = k1++;
                }
            }
        }

        // Добавление необходимых ребер слоя
        graph[i] = vector<vector<pair<int, int>>>(vertexes[i]);
        for (int j = 0; j < vertexes[i - 1]; ++j) {
            int with_removed = j;
            if (removed_v >= 0) {
                int before = (1 << removed_v) - 1;
                int before_inc = (1 << (removed_v + 1)) - 1;
                int after = ~before_inc;
                with_removed = (with_removed & before) + ((with_removed & after) >> 1);
            }

            int nxt_1 = with_removed;
            int nxt_2 = with_removed;
            if (inserted_v >= 0) {
                int before = (1 << inserted_v) - 1;
                int after = ~before;
                nxt_1 = (with_removed & before) + ((with_removed & after) << 1);
                nxt_2 = (with_removed & before) + ((with_removed & after) << 1) + (1 << inserted_v);
            }
            graph[i][nxt_1].push_back({j, edge_num(i, nxt_1, p, activity, removed_v >= 0 ? (j >> removed_v) & 1 : -1)});
            if (nxt_2 != nxt_1) {
                graph[i][nxt_2].push_back({j, edge_num(i, nxt_2, p, activity, removed_v >= 0 ? (j >> removed_v) & 1 : -1)});
            }
        }
    }

    for (int vertex : vertexes) {
        cout << vertex << " ";
    }
    cout << endl;

    string s;
    while (cin >> s) {
        if (s == "Encode") {
            vector<int> input(k);
            for (int i = 0; i < k; ++i) {
                cin >> input[i];
            }
            vector<int> result = encode(n, input, g);
            for (int r : result) {
                cout << r << " ";
            }
        } else if (s == "Decode") {
            vector<double> input(n);
            for (int i = 0; i < n; ++i) {
                cin >> input[i];
            }
            vector<int> result = decode(n, input, graph);
            for (int r : result) {
                cout << r << " ";
            }
        } else if (s == "Simulate") {
            int it, max_err;
            double div;
            cin >> div >> it >> max_err;
            double result = simulate(n, k, it, max_err, div, g, graph);
            cout << result;
        }
        cout << endl;
    }

    return 0;
}
