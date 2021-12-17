#include <bits/stdc++.h>
#define endl "\n"

using namespace std;

typedef long long ll;
typedef long double ld;

mt19937 gen;
uniform_int_distribution<int> uniform(0, 1);
vector<int> direct_;
vector<int> reverse_;
int n, m, k, d, bs_s;
const int BS_SIZE = 256;

// Сумма полиномов
inline int add(int a, int b) {
    return a ^ b;
}

// Обратный полином по модулю
inline int modulo_inverse(int a) {
    return direct_[(n - reverse_[a]) % n];
}

// Возведение полинома в степень
int pw(int a, int b) {
    if (a == 0) {
        return 0;
    } else {
        return direct_[(reverse_[a] * b) % n];
    }
}

// Перемножение полиномов
int multiply(int a, int b) {
    if (a == 0 || b == 0) {
        return 0;
    } else {
        return direct_[(reverse_[a] + reverse_[b]) % n];
    }
}

// Сложение многочленов
vector<int> polypolynomial_add(const vector<int> &a, const vector<int> &b) {
    vector<int> result = a;
    if (b.size() > a.size()) {
        result.resize(b.size(), 0);
    }
    for (int i = 0; i < b.size(); ++i) {
        result[i] = add(result[i], b[i]);
    }
    return result;
}

// Умножение многочлена на x^b (сдвиг вектора коэффициентов)
vector<int> polypolynomial_shift(const vector<int> &a, int b) {
    if (b > 0) {
        vector<int> result(a.size() + b);
        for (int i = 0; i < a.size(); ++i) {
            result[i + b] = a[i];
        }
        return result;
    } else {
        return a;
    }
}

// Умножение многочлена на константу
void _polypolynomial_const_mul(vector<int> &a, int b) {
    for (int & i : a) {
        i = multiply(b, i);
    }
}

// Умножение многочлена на константу
vector<int> polypolynomial_const_mul(const vector<int> &a, int b) {
    vector<int> result = a;
    _polypolynomial_const_mul(result, b);
    return result;
}

// Перемножение многочленов
vector<int> polypolynomial_multiply(const vector<int> &a, const vector<int> &b) {
    vector<int> result(a.size() + b.size() - 1);
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] > 0) {
            for (int j = 0; j < b.size(); ++j) {
                result[i + j] = add(result[i + j], multiply(a[i], b[j]));
            }
        }
    }
    return result;
}

// Вывод значимых битов многочлена
void bit_print(ostream &os, const bitset<BS_SIZE> &a, int b) {
    for (int i = 0; i < b; ++i) {
        os << a[i] << " ";
    }
    os << endl;
}

// Умножение многочлена на x^b (сдвиг вектора коэффициентов)
bitset<BS_SIZE> bit_shift(const bitset<BS_SIZE> &a, int b) {
    return a << b;
}

// Сумма многочленов
inline void _bit_add(bitset<BS_SIZE> &a, const bitset<BS_SIZE> &b) {
    a ^= b;
}

// Поиск старшего бита многочлена
int highest_bit(const bitset<BS_SIZE> &a) {
    for (int i = BS_SIZE - 1; i >= 0; --i) {
        if (a.test(i)) {
            return i;
        }
    }
    return -1;
}

// Подсчет остатка от деления многочленов
bitset<BS_SIZE> bit_mod(bitset<BS_SIZE> a, const bitset<BS_SIZE> &b) {
    int b_h = highest_bit(b);
    int a_h = BS_SIZE - 1;
    while (true) {
        while (a_h >= 0 && !a[a_h]) --a_h;
        if (a_h < b_h) break;
        int s = a_h - b_h;
        _bit_add(a, bit_shift(b, s));
    }
    return a;
}

// Кодирование
void encode(bitset<BS_SIZE> &input, const bitset<BS_SIZE> &bs) {
    input <<= n - k;
    _bit_add(input, bit_mod(input, bs));
}

// Декодирование
void decode(bitset<BS_SIZE> &input) {
    // Вычисление синдрома
    vector<int> tmm(d - 1);
    for (int i = 1; i < d; ++i) {
        if (i % 2 == 0) {
            tmm[i - 1] = multiply(tmm[i / 2 - 1], tmm[i / 2 - 1]);
        } else {
            tmm[i - 1] = 0;
            for (int j = 0; j < n; ++j) {
                tmm[i - 1] = add(tmm[i - 1], multiply(input[j], pw(direct_[i], j)));
            }
        }
    }

    // Вычисление многочлена локаторов ошибок
    vector<int> A(1, 1);
    vector<int> B(1, 1);
    int L = 0;
    int m = 0;
    for (int r = 1; r < d; ++r) {
        int delta = 0;
        for (int j = 0; j <= L; ++j) {
            delta = add(delta, multiply(A[j], tmm[r - 1 - j]));
        }
        if (delta != 0) {
            vector<int> shft = polypolynomial_shift(B, r - m);
            _polypolynomial_const_mul(shft, delta);
            vector<int> T = polypolynomial_add(A, shft);
            if (2 * L <= r - 1) {
                B = polypolynomial_const_mul(A, modulo_inverse(delta));
                A = T;
                L = r - L;
                m = r;
            } else {
                A = T;
            }
        }
    }

    for (int j = 1; j <= n; ++j) {
        // Проверка потенциального корня
        int res = 0;
        for (int it = 0; it < A.size(); ++it) {
            res = add(res, multiply(A[it], pw(j, it)));
        }
        if (res == 0) {
            input.flip((n - reverse_[j]) % n);
        }
    }
}

// Симуляция
double simulate(int it, int max_err, double div, const bitset<BS_SIZE> &bs) {
    uniform_real_distribution<double> distribution(0, 1);
    int cnt_err = 0;
    bitset<BS_SIZE> original;
    for (int i = 0; i < it; ++i) {
        // Генерация вектора
        original.reset();
        for (int j = 0; j < k; ++j) {
            original[j] = uniform(gen) == 1;
        }

        // Кодирование
        encode(original, bs);

        // Внесение шума
        bitset<BS_SIZE> with_noise = original;
        for (int j = 0; j < n; ++j) {
            if (distribution(gen) < div) {
                with_noise.flip(j);
            }
        }

        // Декодирование
        decode(with_noise);
        if (original != with_noise) {
            ++cnt_err;
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

    int pp;
    cin >> n >> pp >> d;

    m = 0;
    while (n > (1 << m) - 1) {
        ++m;
    }

    // Вычисление отображений степень <-> многочлен
    direct_.resize(n);
    reverse_.resize(n + 1);
    for (int i = 0; i < n; ++i) {
        int v = bit_mod(bit_shift(1, i), pp).to_ulong();
        direct_[i] = v;
        reverse_[v] = i;
    }

    // Вычисление пораждающего многочлена
    vector<bool> used(n);
    vector<int> mm(1, 1);
    for (int i = 1; i < d; ++i) {
        // Используем только уникальные минимальные многочлены
        if (!used[i]) {
            int cur = i;
            vector<int> pows;
            while (!(!pows.empty() && pows.front() == cur)) {
                used[cur] = true;
                pows.push_back(cur);
                cur = (cur << 1) % n;
            }
            for (int p : pows) {
                mm = polypolynomial_multiply(mm, vector<int>{direct_[p], 1});
            }
        }
    }
    bs_s = mm.size() - 1;
    bitset<BS_SIZE> bs;
    for (int i = 0; i < mm.size(); ++i) {
        bs.set(i, mm[i] == 1);
    }

    k = n - bs_s;
    cout << k << endl;

    bit_print(cout, bs, mm.size());

    bitset<BS_SIZE> input_enc;
    bitset<BS_SIZE> input_dec;
    string s;
    while (cin >> s) {
        if (s == "Encode") {
            input_enc.reset();
            for (int i = 0; i < k; ++i) {
                int a;
                cin >> a;
                input_enc.set(i, a == 1);
            }
            encode(input_enc, bs);
            bit_print(cout, input_enc, n);
        } else if (s == "Decode") {
            input_dec.reset();
            for (int i = 0; i < n; ++i) {
                int a;
                cin >> a;
                input_dec.set(i, a == 1);
            }
            decode(input_dec);
            bit_print(cout, input_dec, n);
        } else if (s == "Simulate") {
            int it, max_err;
            double div;
            cin >> div >> it >> max_err;
            double result = simulate(it, max_err, div, bs);
            cout << result << endl;
        }
    }

    return 0;
}
