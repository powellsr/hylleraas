#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>

int binomial(int c, int n) {
//    if (c < 0)
//        std::cout << "c: " << c << std::endl;
//    if (n < 0)
//        std::cout << "n: " << n << std::endl;
    if (n == 0 || n == c)
        return 1;
    return binomial(c-1, n-1) + binomial(c-1, n);
}

int factorial(int n) {
    assert(n >= 0);
    return (n == 1 || n == 0) ? 1 : factorial(n-1) * n;
}

double k_pfac(double pfac, int n, int l, int m, double alpha,
         double beta, double gamma) {
    if (pfac == 0.0)
        return 0;
    double sum_term = 0.0;
    for (std::size_t a = 0; a != n+2; ++a) {
        for (std::size_t b = 0; b != l+2; ++b) {
            for (std::size_t c = 0; c != m+2; ++c) {
                sum_term += binomial(l + 1 - b + a, a) *
                            binomial(m + 1 - c + b, b) * binomial(n + 1 - a + c, c)
                            / (pow(alpha + beta, l - b + a + 2)
                               * pow(alpha + gamma, n - a + c + 2)
                               * pow(beta + gamma, m - c + b + 2));
            }
        }
    }

    return pfac * 16.0 * M_PI * M_PI * factorial(n+1) * factorial(l+1) * factorial(m+1) * sum_term;
}

double k(int n, int l, int m, double alpha,
        double beta, double gamma) {
    
    double sum_term = 0.0;
    for (std::size_t a = 0; a != n+2; ++a) {
        for (std::size_t b = 0; b != l+2; ++b) {
            for (std::size_t c = 0; c != m+2; ++c) {
                sum_term += binomial(l + 1 - b + a, a) * 
                    binomial(m + 1 - c + b, b) * binomial(n + 1 - a + c, c) 
                    / (pow(alpha + beta, l - b + a + 2)
                    * pow(alpha + gamma, n - a + c + 2) 
                    * pow(beta + gamma, m - c + b + 2));
            }
        }
    }

    return 16.0 * M_PI * M_PI * factorial(n+1) * factorial(l+1) * factorial(m+1) * sum_term;
}

double s(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma) {
    std::cout << "ni: " << ni << " nj: " << nj << " li: " << li << " lj: " << lj << " mi: " << mi << " mj: " << mj << std::endl;
    std::cout << "alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl;
    return k(ni + nj, li + lj, mi + mj, alpha, beta, gamma);
}

double v_ne(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma,
        int z) {
    return -1 * z * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma));
}

double v_ee(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma) {
    return k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);
}

double t_e(int ni, int li, int mi, int nj,
           int lj, int mj, double alpha, double beta, double gamma) {
    return (alpha * alpha + beta * beta + gamma * gamma) * s(ni, li, mi, nj, lj, mj, alpha, beta, gamma) / 8.0
        + (nj * alpha / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + (lj * beta / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma) 
        + (mj * gamma) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma) 
//        - (nj * (nj - 1) / 2) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
        - k_pfac((nj * (nj - 1) / 2), ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
//        - (lj * (lj - 1) / 2) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
        - k_pfac((lj * (lj - 1) / 2), ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
//
        - k_pfac((mj * (mj - 1)), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        + (alpha / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + (beta / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + (gamma) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
//        - (nj) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
        - k_pfac((nj), ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
//        - (lj) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
        - k_pfac((lj), ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
//        - (2 * mj) * k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        - k_pfac((2 * mj), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        - (alpha * gamma / 8) * (k(ni + nj - 1, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj + 1, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj - 1, li + lj + 2, mi + mj - 1, alpha, beta, gamma))
        - (beta * gamma / 8) * (k(ni + nj, li + lj - 1, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj, li + lj + 2, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 1, mi + mj - 1, alpha, beta, gamma))
//        + (nj * gamma / 4) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma))
        + (nj * gamma / 4) * (k_pfac((nj * gamma / 4), ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((nj * gamma / 4), ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k_pfac((nj * gamma / 4), ni + nj - 2, li + lj + 2, mi + mj - 1, alpha, beta, gamma))
//        + (mj * alpha / 4) * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
        + (k_pfac((mj * alpha / 4), ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
            + k_pfac((mj * alpha / 4), ni + nj + 1, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((mj * alpha / 4), ni + nj - 1, li + lj + 2, mi + mj - 2, alpha, beta, gamma))
//        - (nj * mj / 2) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
        - (k_pfac((nj * mj / 2), ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((nj * mj / 2), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((nj * mj / 2), ni + nj - 2, li + lj + 2, mi + mj - 2, alpha, beta, gamma))
//        + (lj * gamma / 4) * (k(ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma)
          + (k_pfac((lj * gamma / 4), ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((lj * gamma / 4), ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k_pfac((lj * gamma / 4), ni + nj + 2, li + lj - 2, mi + mj - 1, alpha, beta, gamma))
//        + (mj * beta / 4) * (k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + (k_pfac((mj * beta / 4), ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
            + k_pfac((mj * beta / 4), ni + nj, li + lj + 1, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((mj * beta / 4), ni + nj + 2, li + lj - 1, mi + mj - 2, alpha, beta, gamma))
//        - (lj * mj / 2) * (k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
          - (k_pfac((lj * mj / 2), ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
            + k_pfac((lj * mj / 2), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((lj * mj / 2), ni + nj + 2, li + lj - 2, mi + mj - 2, alpha, beta, gamma));
}

Eigen::MatrixXd s_matrix_builder(std::vector<std::pair<std::vector<int>, std::vector<double> > > basis) {
    Eigen::MatrixXd result;
    result.setZero(basis.size(), basis.size());
    for (std::size_t i = 0; i != basis.size(); ++i) {
        for (std::size_t j = 0; j != basis.size(); ++j) {
            //std::vector<int> q_nos_i, q_nos_j;
            //std::vector<double> exp_i, exp_j;
            //auto [q_nos_i, exp_i] = basis[i];
            //auto [q_nos_j, exp_j] = basis[j];
            std::vector<int> q_nos_i = basis[i].first;
            std::vector<double> exp_i = basis[i].second;
            std::vector<int> q_nos_j = basis[j].first;
            result(i,j) = s(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2]);
        }
    }
    return result;
}

Eigen::MatrixXd h_matrix_builder(std::vector<std::pair<std::vector<int>, std::vector<double> > > basis, int z) {
    Eigen::MatrixXd result;
    result.setZero(basis.size(), basis.size());
    for (std::size_t i = 0; i != basis.size(); ++i) {
        for (std::size_t j = 0; j != basis.size(); ++j) {
            //std::vector<int> q_nos_i, q_nos_j;
            //std::vector<double> exp_i, exp_j;
            //auto [q_nos_i, exp_i] = basis[i];
            //auto [q_nos_j, exp_j] = basis[j];
            std::vector<int> q_nos_i = basis[i].first;
            std::vector<double> exp_i = basis[i].second;
            std::vector<int> q_nos_j = basis[j].first;
            result(i,j) = v_ne(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2], z)
                + v_ee(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2])
                + t_e(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2]);
        }
    }
    return result;
}

int main() {
    std::vector<std::pair<std::vector<int>, std::vector<double> > > basis =
            {std::pair<std::vector<int>, std::vector<double> >(std::vector<int>{0,0,0}, std::vector<double>{1.6875, 1.6875, 0.0})};
    Eigen::MatrixXd S = s_matrix_builder(basis);
    std::cout << S << std::endl;

    Eigen::MatrixXd H = h_matrix_builder(basis, 2);
    std::cout << H << std::endl;

    return 0;
}
