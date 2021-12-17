#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>

int binomial(int n, int k) {
//    if (c < 0)
//        std::cout << "c: " << c << std::endl;
//    if (n < 0)
//        std::cout << "n: " << n << std::endl;
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
    return binomial(n-1, k-1) + binomial(n-1, k);
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
    //std::cout << "n,l,m: " << n << "," << l << "," << m << std::endl;
    //std::cout << "alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl;
    for (std::size_t a = 0; a != (n+2); ++a) {
        for (std::size_t b = 0; b != (l+2); ++b) {
            for (std::size_t c = 0; c != (m+2); ++c) {
                double sum_term_part = binomial(l + 1 - b + a, a) *
                            binomial(m + 1 - c + b, b) * binomial(n + 1 - a + c, c)
                            / (pow((alpha + beta)*2, l - b + a + 2)
                               * pow((alpha + gamma)*2, n - a + c + 2)
                               * pow((beta + gamma)*2, m - c + b + 2));
                //std::cout << "a,b,c: " << a << "," << b << "," << c << " sum term contrib: " << sum_term_part << std::endl;
                sum_term += sum_term_part;
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
                    / (pow((alpha + beta)*2, l - b + a + 2)
                    * pow((alpha + gamma)*2, n - a + c + 2)
                    * pow((beta + gamma)*2, m - c + b + 2));
            }
        }
    }

    return 16.0 * M_PI * M_PI * factorial(n+1) * factorial(l+1) * factorial(m+1) * sum_term;
}

double s(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma) {
    //std::cout << "ni: " << ni << " nj: " << nj << " li: " << li << " lj: " << lj << " mi: " <<
    //mi << " mj: " << mj << std::endl;
    //std::cout << "alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl;
    return k_pfac(1, ni + nj, li + lj, mi + mj, alpha, beta, gamma);
}

double v_ne(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma,
        int z) {
    return -z * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
        + k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma));
}

double v_ee(int ni, int li, int mi, int nj,
        int lj, int mj, double alpha, double beta, double gamma) {
    return k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);
}

double t_e(int ni, int li, int mi, int nj,
           int lj, int mj, double alpha, double beta, double gamma) {
    return -((alpha*2) * (alpha*2) + (beta*2) * (beta*2) + (gamma*2) * (gamma*2)) * s(ni, li, mi, nj, lj, mj, alpha, beta, gamma) / 8.0
        + (nj * (alpha*2) / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
        + (lj * (beta*2) / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + (mj * (gamma*2)) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
//        - (nj * (nj - 1) / 2) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
        - k_pfac((nj * (nj - 1) / 2), ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
//        - (lj * (lj - 1) / 2) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
        - k_pfac((lj * (lj - 1) / 2), ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
//
        - k_pfac((mj * (mj - 1)), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        + ((alpha*2) / 2) * k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
        + ((beta*2) / 2) * k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + ((gamma*2)) * k(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
//        - (nj) * k(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
        - k_pfac((nj), ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma)
//        - (lj) * k(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
        - k_pfac((lj), ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma)
//        - (2 * mj) * k(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        - k_pfac((2 * mj), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
        - ((alpha*2) * (gamma*2) / 8) * (k(ni + nj - 1, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj + 1, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj - 1, li + lj + 2, mi + mj - 1, alpha, beta, gamma))
        - ((beta*2) * (gamma*2) / 8) * (k(ni + nj, li + lj - 1, mi + mj + 1, alpha, beta, gamma)
            + k(ni + nj, li + lj + 2, mi + mj - 1, alpha, beta, gamma)
            - k(ni + nj + 2, li + lj - 1, mi + mj - 1, alpha, beta, gamma))
//        + (nj * gamma / 4) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma))
        + (k_pfac((nj * (gamma*2) / 4), ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((nj * (gamma*2) / 4), ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k_pfac((nj * (gamma*2) / 4), ni + nj - 2, li + lj + 2, mi + mj - 1, alpha, beta, gamma))
//        + (mj * alpha / 4) * (k(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
        + (k_pfac((mj * (alpha*2) / 4), ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma)
            + k_pfac((mj * (alpha*2) / 4), ni + nj + 1, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((mj * (alpha*2) / 4), ni + nj - 1, li + lj + 2, mi + mj - 2, alpha, beta, gamma))
//        - (nj * mj / 2) * (k(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
        - (k_pfac((nj * mj / 2), ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((nj * mj / 2), ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((nj * mj / 2), ni + nj - 2, li + lj + 2, mi + mj - 2, alpha, beta, gamma))
//        + (lj * gamma / 4) * (k(ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma)
          + (k_pfac((lj * (gamma*2) / 4), ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma)
            + k_pfac((lj * (gamma*2) / 4), ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma)
            - k_pfac((lj * (gamma*2) / 4), ni + nj + 2, li + lj - 2, mi + mj - 1, alpha, beta, gamma))
//        + (mj * beta / 4) * (k(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
        + (k_pfac((mj * (beta*2) / 4), ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma)
            + k_pfac((mj * (beta*2) / 4), ni + nj, li + lj + 1, mi + mj - 2, alpha, beta, gamma)
            - k_pfac((mj * (beta*2) / 4), ni + nj + 2, li + lj - 1, mi + mj - 2, alpha, beta, gamma))
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
            //std::cout << "S(" << i << "," << j << ") = s(" << q_nos_i[0] << "," <<  q_nos_i[1] << "," << q_nos_i[2] << "," << q_nos_j[0] << "," <<  q_nos_j[1] << "," << q_nos_j[2] << "," << exp_i[0] << "," << exp_i[1] << "," << exp_i[2] << std::endl;
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
            result(i,j) = (v_ne(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2], z)
                + v_ee(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2])
                + t_e(q_nos_i[0], q_nos_i[1], q_nos_i[2], q_nos_j[0], q_nos_j[1], q_nos_j[2], exp_i[0], exp_i[1], exp_i[2]));
            if (result(i,j) == nan(0)) {
                result(i,j) = 0;
            }
        }
    }
    return result;
}

//template<typename T,typename U>
//std::pair<T, U>
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> hylleraas(std::vector<std::pair<std::vector<int>, std::vector<double> > > basis, int z) {
    auto H = h_matrix_builder(basis, z);
    auto S = s_matrix_builder(basis);

    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
    ges.compute(H,S);

    Eigen::MatrixXd evecs = ges.eigenvectors().real();
    Eigen::MatrixXd evals = ges.eigenvalues().real();

    return std::pair(evecs, evals);
}

std::vector<std::pair<std::vector<int>, std::vector<double> > > generate_basis(std::size_t N, double alpha, double gamma) {
    std::vector<std::pair<std::vector<int>, std::vector<double> > > result;
    for (int n = 0; n != N; ++n) {
        for (int l = 0; l != (N-n); ++l) {
            for (int m = 0; m != (N-n-l); ++m) {
                if ((n + l + m) <= N)
                    result.push_back(std::pair<std::vector<int>, std::vector<double>>{std::vector<int>{n,l,m}, std::vector<double>{alpha, alpha, gamma}});
            }
        }
    }
    return result;
}

template<typename T,typename U>
std::pair<T, U> hylleraas_generic(int N, int z, double alpha, double gamma) {
    auto basis = generate_basis(N, alpha, gamma);
    auto [evecs, evals] = hylleraas(basis, z); //<Eigen::MatrixXd, Eigen::MatrixXd>
    return std::pair(evecs, evals);
}

int main() {
    std::vector<std::pair<std::vector<int>, std::vector<double> > > basis =
            {std::pair(std::vector<int>{0,0,0}, std::vector<double>{1.6875, 1.6875, 0.0})};
    Eigen::MatrixXd S = s_matrix_builder(basis);
    std::cout << "S: \n" << S << std::endl;

    Eigen::MatrixXd H = h_matrix_builder(basis, 2);
    std::cout << "H: \n" << H << std::endl;

//    double k_test = k(0, 0, 0, 0.5, 0.5, 0.5);
//    std::cout << "k(-1,-1,-1,0.5,0.5,0.5): " << k_test << std::endl;

//    double k_test_pfac = k_pfac(1,0, 0, 0, 1, 1, 1);
//    std::cout << "k(1,0,0,0,1,1,1): " << k_test_pfac << std::endl;
//    double k_test_pfac2 = k_pfac(1,0, 0, 0, 1.6875, 1.6875, 0);
//    std::cout << "k(1,0,0,0,1.6875,1.6875,0): " << k_test_pfac2 << std::endl;

//    int binom_test = binomial(2,1);
//    std::cout << "binom test C(2,1)\n" << binom_test;

    std::vector<std::pair<std::vector<int>, std::vector<double> > > func3_basis =
            {std::pair(std::vector<int>{0,0,0}, std::vector<double>{1.8, 1.8, 0.0}), std::pair(std::vector{1,1,0}, std::vector{1.8, 1.8, 0.0}), std::pair(std::vector{0,0,1}, std::vector{1.8, 1.8, 0.0})};
    Eigen::MatrixXd S_3func = s_matrix_builder(func3_basis);
    std::cout << "S_3func: \n" << S_3func << std::endl;

    Eigen::MatrixXd H_3func = h_matrix_builder(func3_basis, 2);
    std::cout << "H_3func: \n" << H_3func << std::endl;

    auto [evecs_3func, evals_3func] = hylleraas(func3_basis, 2);
    std::cout << "3basis evecs: \n" << evecs_3func << "\n\n3basis evals: \n" << evals_3func << std::endl;

    //He+
//    std::vector<std::pair<std::vector<int>, std::vector<double> > > func3_basis_he_plus =
//            {std::pair(std::vector<int>{0,0,0}, std::vector<double>{1.8, 0.0, 0.0}), std::pair(std::vector{1,1,0}, std::vector{1.8, 1.8, 0.0}), std::pair(std::vector{0,0,1}, std::vector{1.8, 1.8, 0.0})};
//    Eigen::MatrixXd S_3func_he_plus = s_matrix_builder(func3_basis);
//    std::cout << "S_3func he+: \n" << S_3func_he_plus << std::endl;
//
//    Eigen::MatrixXd H_3func_he_plus = h_matrix_builder(func3_basis_he_plus, 2);
//    std::cout << "H_3func he+: \n" << H_3func_he_plus << std::endl;
//
//    auto [evecs_3func_he_plus, evals_3func_he_plus] = hylleraas(func3_basis_he_plus, 2);
//    std::cout << "3basis evecs, he+: \n" << evecs_3func_he_plus << "\n\n3basis evals: \n" << evals_3func << std::endl;

    //hylleraas generic
    auto [evecs_gen_func, evals_gen_func] = hylleraas_generic<Eigen::MatrixXd, Eigen::MatrixXd>(3, 2, 1.8, 0.0);
    std::cout << "gen basis evecs: \n" << evecs_gen_func << "\n\ngen basis evals: \n" << evals_gen_func << std::endl;

    std::vector<double> evals; //(15);
    for (size_t i = 1; i != 16; ++i) {
        auto [evecs_gen_func, evals_gen_func] = hylleraas_generic<Eigen::MatrixXd, Eigen::MatrixXd>(i, 2, 1.8, 0.0);
        //std::cout << "gen basis evecs: \n" << evecs_gen_func << "\n\ngen basis evals: \n" << evals_gen_func << std::endl;
        //evals[i] = evals_gen_func.minCoeff();
        evals.push_back(evals_gen_func.minCoeff());
        std::cout << "Completed " << i << "th computation: " << evals_gen_func.minCoeff() << "\n";
    }

    for (std::size_t i = 1; i != 16; ++i) {
        std::cout << "He energy computed with i=" << i << ": " << evals[i] << "\n";
    }

    return 0;
}
