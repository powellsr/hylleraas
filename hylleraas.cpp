#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>

int binomial(std::size_t c, std::size_t n) {
    if (n == 0 || n == c)
        return 1;
    return binomial(c-1, n-1) + binomial(c-1, n);
}

int factorial(std::size_t n) {
    return (n == 1 || n == 0) ? 1 : factorial(n-1) * n;
}

double k(std::size_t n, std::size_t l, std::size_t m, double alpha, 
        double beta, double gamma) {
    
    double sum_term = 0.0;
    for (std::size_t a = 0; a != n+2; ++a) {
        for (std::size_t b = 0; b != l+2; ++b) {
            for (std::size_t c = 0; c != m+2; ++c) {
                sum_term += binomial(l + 1 - b + a, a) * 
                    binomial(m + 1 - c + b, b) * binomial(n + 1 - a + c, c) 
                    / (pow(alpha + beta, l - b + a +2) 
                    * pow(alpha + gamma, n - a + c + 2) 
                    * pow(beta + gamma, m - c + b + 2));
            }
        }
    }

    return 16.0 * M_PI * M_PI * factorial(n+1) * factorial(l+1) * factorial(m+1) * sum_term;
}

double s(std::size_t ni, std::size_t li, std::size_t mi, std::size_t nj, 
        std::size_t lj, std::size_t mj, double alpha, double beta, double gamma) {
    return k(ni + nj, li + lj, mi + mj, alpha, beta, gamma);
}

double v_ne(std::size_t ni, std::size_t li, std::size_t mi, std::size_t nj, 
        std::size_t lj, std::size_t mj, double alpha, double beta, double gamma, 
        std::size_t z) {
    return -1 * z * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma));
}

double v_ee(std::size_t ni, std::size_t li, std::size_t mi, std::size_t nj, 
        std::size_t lj, std::size_t mj, double alpha, double beta, double gamma) {
    return k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);
}

double t_e(std::size_t ni, std::size_t li, std::size_t mi, std::size_t nj, 
        std::size_t lj, std::size_t mj, double alpha, double beta, double gamma) {
    return (alpha * alpha + beta * beta + gamma * gamma) * s(ni, li, mi, nj, lj, mj, alpha, beta, gamma) / 8.0
        + (nj * alpha / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + (lj * beta / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma) 
        + (mj * gamma) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma) 
        - (nj * (nj - 1) / 2) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma) 
        - (lj * (lj - 1) / 2) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma) 
        - (mj * (mj - 1)) * k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma) 
        + (alpha / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) 
        + (beta / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + (gamma) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
        - (nj) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
        - (lj) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
        - (2 * mj) * k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        - (alpha * gamma / 8) * (k(ni + nj - 1, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj + 1, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj - 1, li + lj + 2, mi + mj - 1, alpha, beta, gamma))
        - (beta * gamma / 8) * (k(ni + nj, li + lj - 1, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj, li + lj + 2, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 1, mi + mj - 1, alpha, beta, gamma))
        + (nj * gamma / 4) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma))
            + k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj - 2, li + lj + 2, mi + mj - 1, alpha, beta, gamma)
        + (mj * alpha / 4) * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma))
            + k(ni + nj + 1, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k(ni + nj - 1, li + lj + 2, mi + mj - 2, alpha, beta, gamma)
        - (nj * mj / 2) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma))
            + k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k(ni + nj - 2, li + lj + 2, mi + mj - 2, alpha, beta, gamma)
        + (lj * gamma / 4) * (k(ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma))
            + k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 2, mi + mj - 1, alpha, beta, gamma)
        + (mj * beta / 4) * (k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma))
            + k(ni + nj, li + lj + 1, mi + mj - 2, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 1, mi + mj - 2, alpha, beta, gamma)
        - (lj * mj / 2) * (k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma))
            + k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 2, mi + mj - 2, alpha, beta, gamma);
}

Eigen::MatrixXd s_matrix_builder(std::vector<std::pair<std::vector<int>, std::vector<double>>> basis) {
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
}

Eigen::MatrixXd h_matrix_builder(std::vector<std::pair<std::vector<int>, std::vector<double>>> basis, int z) {
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
}

int main() {
    std::vector<std::pair<std::vector<int>, std::vector<double>>> basis = {std::pair<std::vector<int>, std::vector<double>>{std::vector<int>{0,0,0}, std::vector<double>{1.6875, 1.6875, 0.0}}};
    Eigen::MatrixXd S = s_matrix_builder(basis);
    std::cout << S << std::endl;

    return 0;
}