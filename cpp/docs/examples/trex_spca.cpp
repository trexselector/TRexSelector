// Example program of the guide page "T-Rex SPCA". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_spca/trex_spca.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

using trex::trex_selector_methods::ad::trex_spca::SPCAMode;
using trex::trex_selector_methods::ad::trex_spca::TRexSPCA;
using trex::trex_selector_methods::ad::trex_spca::TRexSPCAControlParameter;
using trex::utils::memmap::AccessMode;
using trex::utils::memmap::MemoryMappedMatrix;
//! [includes]

//! [helpers]
void print_active_sets(const std::vector<std::vector<Eigen::Index>>& sets) {
    for (std::size_t m = 0; m < sets.size(); ++m) {
        std::cout << "component " << m + 1 << ":";
        for (Eigen::Index j : sets[m]) std::cout << ' ' << j;
        std::cout << '\n';
    }
}
//! [helpers]

int main() {
    DocSections doc("trex_spca", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(77);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 200, p = 60;

    Eigen::MatrixXd Z_true = Eigen::MatrixXd::NullaryExpr(n, 2, [&] { return normal(gen); });
    Eigen::MatrixXd V_true = Eigen::MatrixXd::Zero(p, 2);
    V_true.col(0).segment(0, 10).setConstant(1.0 / std::sqrt(10.0));
    V_true.col(1).segment(30, 10).setConstant(1.0 / std::sqrt(10.0));
    Eigen::MatrixXd noise = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::MatrixXd X = 6.0 * Z_true.col(0) * V_true.col(0).transpose()
                      + 3.0 * Z_true.col(1) * V_true.col(1).transpose() + noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    //! [data]

    doc.begin("run");
    //! [run]
    TRexSPCAControlParameter spca_ctrl;
    spca_ctrl.gvs_ctrl.trex_ctrl.K = 10;

    TRexSPCA spca(X_map, /*M=*/2, /*tFDR=*/0.1, spca_ctrl, /*seed=*/9, /*verbose=*/false);
    const auto& res = spca.select();
    print_active_sets(res.active_sets);
    //! [run]

    doc.begin("fields");
    //! [fields]
    std::cout << "Z: " << res.Z.rows() << " x " << res.Z.cols() << '\n';
    std::cout << "V: " << res.V.rows() << " x " << res.V.cols() << '\n';
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "adjusted_ev:   " << res.adjusted_ev.transpose() << '\n';
    std::cout << "cumulative_ev: " << res.cumulative_ev.transpose() << '\n';
    //! [fields]

    doc.begin("loadings");
    //! [loadings]
    for (Eigen::Index m = 0; m < 2; ++m) {
        std::cout << "loading " << m + 1 << ", non-zero at:";
        for (Eigen::Index j = 0; j < p; ++j) {
            if (res.V(j, m) != 0.0) std::cout << ' ' << j;
        }
        const double on_support = res.V.col(m).segment(m == 0 ? 0 : 30, 10).squaredNorm();
        std::cout << "\n   squared norm on the true support: " << on_support << '\n';
    }
    std::ofstream v_csv(std::string(DOC_OUTPUT_DIR) + "/trex_spca.figure_loadings.csv");
    for (Eigen::Index j = 0; j < p; ++j) {   // the data of the figure: variable, loading 1, loading 2
        v_csv << j << ',' << res.V(j, 0) << ',' << res.V(j, 1) << '\n';
    }
    //! [loadings]

    doc.begin("thresholded");
    //! [thresholded]
    TRexSPCAControlParameter thr_ctrl;
    thr_ctrl.mode = SPCAMode::Thresholded;
    thr_ctrl.gvs_ctrl.trex_ctrl.K = 10;

    TRexSPCA thr(X_map, 2, 0.1, thr_ctrl, 9, false);
    const auto& thr_res = thr.select();
    print_active_sets(thr_res.active_sets);
    std::cout << "cumulative_ev: " << thr_res.cumulative_ev.transpose() << '\n';
    //! [thresholded]

    doc.begin("mmap");
    //! [mmap]
    namespace fs = std::filesystem;
    const fs::path workdir = fs::temp_directory_path() / "trex_guide_trex_spca";
    fs::create_directories(workdir);
    const std::string path = (workdir / "X.bin").string();
    {
        MemoryMappedMatrix<double> file(path, n, p, AccessMode::ReadWrite);
        file.getMap() = X;                                 // writes the design into the file
        Eigen::Map<Eigen::MatrixXd> X_file = file.getMap();

        TRexSPCA mapped(X_file, 2, 0.1, spca_ctrl, 9, false);
        const auto& mapped_res = mapped.select();
        const bool same_sets = mapped_res.active_sets == res.active_sets;
        std::cout << std::boolalpha;
        std::cout << "same active sets: " << same_sets << '\n';
        std::cout << "file restored:    " << file.getMap().isApprox(X) << '\n';
    }
    fs::remove_all(workdir);
    //! [mmap]

    return 0;
}
