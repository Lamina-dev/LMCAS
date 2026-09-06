#pragma once

#include <memory>
#include <variant>

namespace LMCAS {

class SymbolicExpr;



using ProofExprPtr = std::shared_ptr<SymbolicExpr>;

struct ByConstructionProof {};
struct ExactNormalizationProof {
    ProofExprPtr normalized_residual;
};
struct ExactResidualProof {
    ProofExprPtr normalized_residual;
};
struct ExactRoundTripProof {
    ProofExprPtr normalized_residual;
};

using ProofCertificate = std::variant<
    ByConstructionProof,
    ExactNormalizationProof,
    ExactResidualProof,
    ExactRoundTripProof>;

template <class T>
struct Verified {
    T value;
    ProofCertificate certificate;
};

} // namespace LMCAS
