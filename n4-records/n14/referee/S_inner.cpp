// Claim S, part 2: branch and bound for the positive roots of F in logarithmic coordinates,
// on the cube |y|_inf <= T0, with interval arithmetic (outward widening after every operation).
//
//   G(y) = sum_i C_i exp(<a_i - b, y>),   y = log x.
// On a box B the function is divided by the positive constant exp(M) (M = largest centre exponent),
// which changes neither the roots nor the sign of the Jacobian determinant.
//
// A box is discarded only if  (1) some linear combination lam.G has an interval enclosure on B that
// excludes 0, or (2) the Krawczyk image K(B) misses B.  A box is declared to hold exactly one root
// if K(B) lies in the interior of B (Krawczyk/Moore: existence and uniqueness).  Everything else is
// bisected.  Approximate (non-rigorous) quantities - the preconditioner Y and the vectors lam - are
// only ever used as *choices*; all enclosures are computed with the interval type below.
//
// Rounding: doubles in round-to-nearest; after each +,-,* the result is widened by one ulp on each
// side (covers the <=1/2 ulp rounding error).  exp() is widened by 4 ulps on each side (glibc's exp
// is documented to be accurate to < 1 ulp).  This is the one non-exact ingredient of the review.
//
// Build: g++ -O2 -std=c++17 -o S_inner S_inner.cpp      (do NOT use -ffast-math)
// Run:   ./S_inner T0 < system.txt    where system.txt = N, then N rows "a1 a2 a3 a4 c1 c2 c3 c4"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <limits>
using namespace std;
static const double INF = numeric_limits<double>::infinity();
static inline double dn(double x) { return nextafter(x, -INF); }
static inline double up(double x) { return nextafter(x, INF); }
struct I {
    double lo, hi;
    I() : lo(0), hi(0) {}
    I(double a) : lo(a), hi(a) {}
    I(double a, double b) : lo(a), hi(b) {}
};
// Any NaN (inf-inf or 0*inf, possible only when an enclosure has overflowed) is turned into the
// whole line, so an overflow can never produce a spuriously narrow enclosure.
static inline I fix(double l, double h) { return (std::isnan(l) || std::isnan(h)) ? I(-INF, INF) : I(l, h); }
static inline I operator+(I a, I b) { return fix(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline I operator-(I a, I b) { return fix(dn(a.lo - b.hi), up(a.hi - b.lo)); }
static inline I operator*(I a, I b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    if (std::isnan(p1) || std::isnan(p2) || std::isnan(p3) || std::isnan(p4)) return I(-INF, INF);
    return fix(dn(min(min(p1, p2), min(p3, p4))), up(max(max(p1, p2), max(p3, p4))));
}
static inline I iexp(I a) {
    double l = exp(a.lo), h = exp(a.hi);
    for (int k = 0; k < 4; k++) { l = dn(l); h = up(h); }
    if (l < 0) l = 0;
    if (std::isnan(l) || std::isnan(h)) return I(0.0, INF);
    return I(l, h);
}
static inline bool has0(I a) { return !(a.lo > 0 || a.hi < 0); }   // NaN-safe: "contains 0" unless proved otherwise

int N;
double A[64][4], C[64][4];   // A[i] = shifted exponent (integers), C[i] = coefficient column
struct Box { double lo[4], hi[4]; };
struct Root { Box b; int sign; };
vector<Root> roots;
vector<Box> unresolved;
long long nboxes = 0, nexcl0 = 0, nexclK = 0;

static bool inv4(const double M[4][4], double Y[4][4]) {   // plain floating point, a choice only
    double a[4][8];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { a[i][j] = M[i][j]; a[i][4 + j] = (i == j); }
    for (int c = 0; c < 4; c++) {
        int p = c;
        for (int r = c + 1; r < 4; r++) if (fabs(a[r][c]) > fabs(a[p][c])) p = r;
        if (fabs(a[p][c]) < 1e-280) return false;
        for (int j = 0; j < 8; j++) swap(a[c][j], a[p][j]);
        double d = a[c][c];
        for (int j = 0; j < 8; j++) a[c][j] /= d;
        for (int r = 0; r < 4; r++) if (r != c) { double f = a[r][c]; for (int j = 0; j < 8; j++) a[r][j] -= f * a[c][j]; }
    }
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) Y[i][j] = a[i][4 + j];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (!isfinite(Y[i][j])) return false;
    return true;
}
static I idet(I M[4][4], int n, int *rows, int col) {      // Laplace expansion, interval
    if (n == 1) return M[rows[0]][col];
    I s(0.0);
    for (int k = 0; k < n; k++) {
        int sub[4], t = 0;
        for (int j = 0; j < n; j++) if (j != k) sub[t++] = rows[j];
        I m = M[rows[k]][col] * idet(M, n - 1, sub, col + 1);
        s = (k % 2 == 0) ? s + m : s - m;
    }
    return s;
}

// returns 0 = excluded, 1 = unique root (sign stored), 2 = undecided
static int examine(Box &B, int &sign, bool &sign_known) {
    double c[4], r[4];
    I Y_[4];
    for (int k = 0; k < 4; k++) { c[k] = 0.5 * (B.lo[k] + B.hi[k]); Y_[k] = I(B.lo[k], B.hi[k]); }
    // exponent enclosures on the box and at the centre
    I E[64], Ec[64];
    double M = -INF;
    for (int i = 0; i < N; i++) {
        I e(0.0), ec(0.0);
        for (int k = 0; k < 4; k++) { e = e + I(A[i][k]) * Y_[k]; ec = ec + I(A[i][k]) * I(c[k]); }
        E[i] = e; Ec[i] = ec;
        M = max(M, ec.hi);
    }
    I m[64], mc[64];
    for (int i = 0; i < N; i++) { m[i] = iexp(E[i] - I(M)); mc[i] = iexp(Ec[i] - I(M)); }
    // (1) zero-order tests
    auto excl = [&](const double lam[4]) {
        I s(0.0);
        for (int i = 0; i < N; i++) {
            I w(0.0);
            for (int k = 0; k < 4; k++) w = w + I(lam[k]) * I(C[i][k]);
            s = s + w * m[i];
        }
        return !has0(s);
    };
    for (int k = 0; k < 4; k++) { double lam[4] = {0, 0, 0, 0}; lam[k] = 1; if (excl(lam)) return 0; }
    {   // dominance: lam with lam.C_i = 1 on the k largest monomials (least-norm), k = 1..4
        int ord[64];
        for (int i = 0; i < N; i++) ord[i] = i;
        sort(ord, ord + N, [&](int a, int b) { return E[a].hi > E[b].hi; });
        for (int k = 1; k <= 4; k++) {
            double Gm[4][4] = {{0}}, Gi[4][4];
            for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) Gm[a][b] = (a == b && a >= k) ? 1 : 0;
            for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) { double s = 0; for (int t = 0; t < 4; t++) s += C[ord[a]][t] * C[ord[b]][t]; Gm[a][b] = s; }
            if (!inv4(Gm, Gi)) continue;
            double lam[4] = {0, 0, 0, 0};
            for (int a = 0; a < k; a++) { double w = 0; for (int b = 0; b < k; b++) w += Gi[a][b]; for (int t = 0; t < 4; t++) lam[t] += w * C[ord[a]][t]; }
            if (excl(lam)) return 0;
        }
    }
    // (2) Krawczyk with preconditioner Y ~ J(c)^-1
    double Jc[4][4] = {{0}}, Y[4][4];
    for (int i = 0; i < N; i++) { double mm = 0.5 * (mc[i].lo + mc[i].hi); for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) Jc[a][b] += C[i][a] * mm * A[i][b]; }
    if (!inv4(Jc, Y)) return 2;
    for (int a = 0; a < 4; a++) if (excl(Y[a])) return 0;
    I Gc[4], JB[4][4];
    for (int a = 0; a < 4; a++) { Gc[a] = I(0.0); for (int b = 0; b < 4; b++) JB[a][b] = I(0.0); }
    for (int i = 0; i < N; i++) for (int a = 0; a < 4; a++) {
        Gc[a] = Gc[a] + I(C[i][a]) * mc[i];
        for (int b = 0; b < 4; b++) JB[a][b] = JB[a][b] + I(C[i][a] * A[i][b]) * m[i];   // C*A product of small integers: exact
    }
    bool inside = true;
    for (int a = 0; a < 4; a++) {
        I yg(0.0);
        for (int b = 0; b < 4; b++) yg = yg + I(Y[a][b]) * Gc[b];
        I K = I(c[a]) - yg;
        for (int b = 0; b < 4; b++) {
            I yj(0.0);
            for (int t = 0; t < 4; t++) yj = yj + I(Y[a][t]) * JB[t][b];
            I e = I(a == b ? 1.0 : 0.0) - yj;
            K = K + e * (Y_[b] - I(c[b]));
        }
        if (K.hi < B.lo[a] || K.lo > B.hi[a]) return 0;      // false if K is NaN/infinite
        if (!(K.lo > B.lo[a] && K.hi < B.hi[a])) inside = false;
    }
    if (!inside) return 2;
    int rows[4] = {0, 1, 2, 3};
    I d = idet(JB, 4, rows, 0);
    sign_known = !has0(d);
    sign = d.lo > 0 ? 1 : -1;
    return 1;
}

static void solve(Box B, int depth) {
    nboxes++;
    int sign = 0; bool sk = false;
    int res = examine(B, sign, sk);
    if (res == 0) return;
    if (res == 1 && sk) { roots.push_back({B, sign}); return; }
    // bisect the widest coordinate (slightly off-centre so that roots do not sit on dyadic faces)
    int k = 0;
    for (int j = 1; j < 4; j++) if (B.hi[j] - B.lo[j] > B.hi[k] - B.lo[k]) k = j;
    if (B.hi[k] - B.lo[k] < 1e-10 || depth > 400) { unresolved.push_back(B); return; }
    double mid = B.lo[k] + 0.4993 * (B.hi[k] - B.lo[k]);
    Box L = B, R = B;
    L.hi[k] = mid; R.lo[k] = mid;
    solve(L, depth + 1); solve(R, depth + 1);
}

int main(int argc, char **argv) {
    double T0 = atof(argv[1]);
    if (scanf("%d", &N) != 1) return 1;
    for (int i = 0; i < N; i++) for (int k = 0; k < 8; k++) { double v; if (scanf("%lf", &v) != 1) return 1; if (k < 4) A[i][k] = v; else C[i][k - 4] = v; }
    Box B;
    for (int k = 0; k < 4; k++) { B.lo[k] = -T0; B.hi[k] = T0; }
    solve(B, 0);
    printf("boxes examined: %lld\n", nboxes);
    printf("unresolved boxes: %zu\n", unresolved.size());
    for (auto &u : unresolved) printf("  U %.17g %.17g %.17g %.17g  w=%.3g\n", u.lo[0], u.lo[1], u.lo[2], u.lo[3], u.hi[0] - u.lo[0]);
    printf("boxes with exactly one root: %zu\n", roots.size());
    int tot = 0;
    for (auto &r : roots) {
        tot += r.sign;
        printf("ROOT %+d", r.sign);
        for (int k = 0; k < 4; k++) printf(" %.17g %.17g", r.b.lo[k], r.b.hi[k]);
        printf("\n");
    }
    int pos = 0; for (auto &r : roots) pos += r.sign > 0;
    printf("positive Jacobian: %d, negative Jacobian: %d, signed sum: %d\n", pos, (int)roots.size() - pos, tot);
    return 0;
}
