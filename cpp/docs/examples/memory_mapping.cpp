// Example program of the guide page "Memory mapping". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <ml_methods/scaler_methods/z_score_scaler.hpp>
#include <ml_methods/model_selection/enet_cv_ccd.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::tsolvers::ad::linear_model::lars_based::TLARS_Solver;
using trex::ml_methods::scaler_methods::ZScoreScaler;
using trex::ml_methods::model_selection::enet_gaussian;
using trex::utils::memmap::AccessMode;
using MappedMatrix = trex::utils::memmap::MemoryMappedMatrix<double>;
namespace fs = std::filesystem;
//! [includes]

namespace {

// The mapped objects live inside this function, so they are unmapped before main() removes
// the folder.
void run(DocSections& doc, const fs::path& workdir) {
    doc.begin("create");
    //! [create]
    std::mt19937_64 gen(5);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, p = 200;
    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    const std::string path = (workdir / "design.bin").string();

    // The constructor creates the file of n * p doubles, the map writes into it
    MappedMatrix mm(path, n, p, AccessMode::ReadWrite);
    mm.getMap() = X;
    std::cout << mm.rows() << " x " << mm.cols() << '\n';
    std::cout << std::boolalpha
              << (fs::file_size(path) == static_cast<std::uintmax_t>(8 * n * p)) << '\n';
    //! [create]

    doc.begin("access");
    //! [access]
    auto view = mm.getMap();                        // an Eigen::Map of the file, no copy
    std::cout << std::fixed << std::setprecision(3) << view.block(0, 0, 2, 3) << '\n';
    std::cout << (view.block(0, 0, 2, 3) == X.block(0, 0, 2, 3)) << '\n';
    mm(0, 0) = 42.0;
    std::cout << mm(0, 0) << '\n';
    mm(0, 0) = X(0, 0);
    std::cout << (view - X).cwiseAbs().maxCoeff() << '\n';

    MappedMatrix ro(path, n, p, AccessMode::ReadOnly);
    try {
        ro(0, 0) = 0.0;
    } catch (const std::runtime_error& e) {
        std::cout << "error: " << e.what() << '\n';
    }
    //! [access]

    doc.begin("stream");
    //! [stream]
    const std::string path2 = (workdir / "streamed.bin").string();
    {
        std::ofstream out(path2, std::ios::binary);
        Eigen::VectorXd column(n);
        const auto bytes = static_cast<std::streamsize>(sizeof(double)) * n;
        for (Eigen::Index j = 0; j < p; ++j) {
            std::mt19937_64 gen_j(1000 + static_cast<unsigned>(j));
            std::normal_distribution<double> normal_j(0.0, 1.0);
            for (Eigen::Index i = 0; i < n; ++i) column[i] = normal_j(gen_j);
            out.write(reinterpret_cast<const char*>(column.data()), bytes);
        }
    }
    MappedMatrix big(path2, n, p, AccessMode::ReadWrite);

    std::mt19937_64 gen_2(1002);
    std::normal_distribution<double> normal_2(0.0, 1.0);
    const auto draw_2 = [&] { return normal_2(gen_2); };
    Eigen::VectorXd column_2 = Eigen::VectorXd::NullaryExpr(n, draw_2);
    std::cout << (big.getMap().col(2) - column_2).cwiseAbs().maxCoeff() << '\n';
    //! [stream]

    doc.begin("select");
    //! [select]
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.head(4).setConstant(1.5);
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    TRexControlParameter ctrl;
    ctrl.K = 10;

    TRexSelector sel_mem(X_map, y_map, 0.1, ctrl, /*seed=*/3, /*verbose=*/false);
    sel_mem.select();
    TRexSelector sel_map(view, y_map, 0.1, ctrl, 3, false);   // the map of the file
    sel_map.select();
    std::cout << (sel_mem.getSelectedIndices() == sel_map.getSelectedIndices()) << '\n';
    std::cout << (view - X).cwiseAbs().maxCoeff() << '\n';   // the original values are back
    //! [select]

    doc.begin("solver");
    //! [solver]
    Eigen::MatrixXd D = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), n, p);

    TLARS_Solver s(view, D_map, y_map);
    s.executeStep(/*T_stop=*/2);
    for (std::size_t j : s.getActivePredictorIndices()) std::cout << j << ' ';
    std::cout << '\n';
    s.restore(view, D_map, y_map);

    ZScoreScaler z;
    z.fit_transform_inplace(view);
    std::cout << (std::abs(view.col(0).mean()) < 1e-12) << '\n';   // the file is centered
    z.inverse_transform_inplace(view);

    enet_gaussian en(/*low_memory=*/true);
    en.fit(view, y, /*alpha=*/1.0, /*n_lambda=*/20);
    std::cout << en.coef().rows() << " x " << en.coef().cols() << '\n';
    //! [solver]
}

}  // namespace

int main() {
    DocSections doc("memory_mapping", DOC_OUTPUT_DIR);
    const fs::path workdir = fs::temp_directory_path() / "trex_guide_memory_mapping";
    fs::create_directories(workdir);
    run(doc, workdir);
    fs::remove_all(workdir);
    return 0;
}
